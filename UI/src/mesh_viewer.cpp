// Mesh viewer (written with Claude): draws the O-grid from mesh::generateOGrid,
// and on top of it the exact NACA airfoil from the CAD equations, so you can see the wall sits on it.
// A panel shows whether the mesh is valid: flipped cells (filled in red), wall deviation, seam gap.
//   buttons            Computed / Profil NACA: show or hide the grid and the exact airfoil
//                      8 ... 128: number of cells nc (regenerates the grid)
//                      NACA M P TT: the airfoil digits, with - / + or click the number to type it
//   mouse wheel        zoom toward the mouse
//   right-click drag   move the view
//   R                  back to the airfoil view
//   F                  show the whole grid, out to the far field

#include <algorithm>
#include <cmath>
#include <vector>
#include "raylib.h"
#define RAYGUI_IMPLEMENTATION               // raygui is a single header: this line puts its code in this file
#include "raygui.h"
#include "OGrid.h"
#include "nacawing/NACA.hpp"

int main()
{
    SetConfigFlags(FLAG_WINDOW_RESIZABLE);  // must come before InitWindow: lets you drag the window edges
    InitWindow(800, 450, "O-grid viewer");  // a second program, so a second window next to naca_viewer
    SetWindowMinSize(640, 360);             // small enough to shrink, big enough for the buttons
    SetTargetFPS(60);
    int screenW = GetScreenWidth(), screenH = GetScreenHeight();   // last known size, to measure a resize

    int digitM = 0, digitP = 0, digitTT = 12;   // start with NACA 0012
    int m = digitM, p = digitP, tt = digitTT;   // what the spinners show: can be invalid for a moment (2012 on the way to 2412)

    const int levels[] = {8, 16, 32, 64, 128};  // the nc choices, one per button
    int level = 3;                          // start at nc = 64
    int nc = 0;
    mesh::StructuredMesh2D grid;            // computed by generateOGrid
    std::vector<Vector2> airfoil;           // exact NACA curve, TE -> lower -> LE -> upper -> TE, screen y = -y

    // Mesh validity, recomputed by load()
    struct Cell { int i, j; };
    std::vector<Cell> badCells;             // cells whose signed area is zero or has the minority sign (flipped)
    double wallError = 0.0;                 // max | |y| - y_t(x) | on the wall j = 0 (symmetric airfoils only)
    double seamGap = 0.0;                   // max distance between at(0, j) and at(nc, j): the two sides of the seam

    // Camera2D: we draw in chord units and the camera turns them into pixels.
    // offset = the screen pixel where target sits, zoom = pixels per chord.
    Camera2D airfoilView{};
    airfoilView.offset = {400.0f, 225.0f};
    airfoilView.target = {0.5f, 0.0f};   // mid-chord in the middle of the window
    airfoilView.zoom = 600.0f;
    Camera2D farView = airfoilView;
    float farRadius = 1.0f;                 // distance from mid-chord to the far field, in chords

    // Rebuilds everything that depends on the airfoil or on nc.
    auto load = [&]() {
        nc = levels[level];
        nacawing::NacaParams naca = nacawing::fromDigits(digitM, digitP, digitTT);
        grid = mesh::generateOGrid(naca, nc);

        // Exact airfoil from the CAD equations, out to the trailing edge of the grid (x_TE = 1.00893...,
        // where y_t = 0): node (0, 0) is the trailing edge, so we take x_TE from there, no second copy.
        const double xTE = grid.at(0, 0).x;
        const int n = 400;
        std::vector<nacawing::Point2D> up = nacawing::upperSurface(naca.m, naca.p, naca.t, n, xTE);
        std::vector<nacawing::Point2D> lo = nacawing::lowerSurface(naca.m, naca.p, naca.t, n, xTE);
        airfoil.clear();
        for (int k = n - 1; k >= 0; --k) airfoil.push_back({static_cast<float>(lo[k].x), static_cast<float>(-lo[k].y)});
        for (int k = 1; k < n; ++k)      airfoil.push_back({static_cast<float>(up[k].x), static_cast<float>(-up[k].y)});

        // Validity 1: signed area of each cell, half the cross product of its diagonals (same as test_ogrid).
        // In a valid grid every cell has the same sign; a cell with the other sign (or zero) is flipped.
        std::vector<double> area;
        int positives = 0;
        for (int j = 0; j + 1 < grid.nj; ++j)
            for (int i = 0; i + 1 < grid.ni; ++i) {
                const auto& a = grid.at(i, j);
                const auto& b = grid.at(i + 1, j);
                const auto& c = grid.at(i + 1, j + 1);
                const auto& d = grid.at(i, j + 1);
                area.push_back(0.5 * ((c.x - a.x) * (d.y - b.y) - (d.x - b.x) * (c.y - a.y)));
                if (area.back() > 0) positives++;
            }
        double sign = (2 * positives >= static_cast<int>(area.size())) ? 1.0 : -1.0;   // the majority sign
        badCells.clear();
        for (int j = 0; j + 1 < grid.nj; ++j)
            for (int i = 0; i + 1 < grid.ni; ++i)
                if (sign * area[i + j * (grid.ni - 1)] <= 0.0) badCells.push_back({i, j});

        // Validity 2: the wall is the airfoil, |y| = y_t(x). Leading edge (i = nc/2) skipped: sqrt(x) is
        // vertical at x = 0, so 1e-17 of rounding on x gives 1e-9 on y_t. x clamped to [0, x_TE] (sqrt(<0) = NaN).
        wallError = 0.0;
        for (int i = 0; i < grid.ni; ++i) {
            if (i == nc / 2) continue;
            const auto& w = grid.at(i, 0);
            double x = std::clamp(w.x, 0.0, xTE);
            wallError = std::max(wallError, std::abs(std::abs(w.y) - nacawing::y_t(naca.t, x)));
        }

        // Validity 3: the seam, the first and last spokes (i = 0 and i = nc) must be the same points.
        seamGap = 0.0;
        for (int j = 0; j < grid.nj; ++j)
            seamGap = std::max(seamGap, std::hypot(grid.at(0, j).x - grid.at(nc, j).x, grid.at(0, j).y - grid.at(nc, j).y));

        // the far field is ~150 chords away (more for thick or cambered airfoils): fit its radius in the window
        nacawing::Point2D far = grid.at(0, grid.nj - 1);
        farRadius = static_cast<float>(std::hypot(far.x - 0.5, far.y));
        farView.zoom = 0.45f * static_cast<float>(GetScreenHeight()) / farRadius;
    };
    load();

    bool showComputed = true;
    bool showAirfoil = true;
    bool editM = false, editP = false, editTT = false;   // true while typing in that number box

    Camera2D cam = airfoilView;

    // Draws every line of a grid. -y because screen y points down.
    auto drawGrid = [](mesh::StructuredMesh2D& g, Color rings, Color spokes, Color wall, float width) {
        auto node = [&](int i, int j) {
            nacawing::Point2D p = g.at(i, j);
            return Vector2{static_cast<float>(p.x), static_cast<float>(-p.y)};
        };
        for (int j = 0; j < g.nj; ++j)              // rings: j fixed, i goes around
            for (int i = 0; i + 1 < g.ni; ++i)
                DrawLineEx(node(i, j), node(i + 1, j), width, rings);
        for (int i = 0; i < g.ni; ++i)              // spokes: i fixed, j goes out
            for (int j = 0; j + 1 < g.nj; ++j)
                DrawLineEx(node(i, j), node(i, j + 1), width, spokes);
        for (int i = 0; i + 1 < g.ni; ++i)          // the wall (j = 0) on top, thicker
            DrawLineEx(node(i, 0), node(i + 1, 0), 2.0f * width, wall);
    };

    while (!WindowShouldClose())            // true when you click X or press Esc
    {
        // Window resized: move every camera's offset by half the size change, so the view stays centred.
        if (IsWindowResized()) {
            float dx = 0.5f * static_cast<float>(GetScreenWidth() - screenW);
            float dy = 0.5f * static_cast<float>(GetScreenHeight() - screenH);
            screenW = GetScreenWidth();
            screenH = GetScreenHeight();
            for (Camera2D* c : {&cam, &airfoilView, &farView}) {
                c->offset.x += dx;
                c->offset.y += dy;
            }
            farView.zoom = 0.45f * static_cast<float>(screenH) / farRadius;   // the far field still fits
        }
        float right = static_cast<float>(GetScreenWidth());   // buttons are placed from the right edge
        float bottom = static_cast<float>(GetScreenHeight());

        // UI scale: 1 for the 800 x 450 start window, 2 at 1600 x 900. Buttons and text are sized for 800 x 450
        // and multiplied by s. S() is for rectangles (float pixels), F() for DrawText font sizes (int pixels).
        float s = std::fmin(right / 800.0f, bottom / 450.0f);
        auto S = [s](float v) { return v * s; };
        auto F = [s](int size) { return static_cast<int>(static_cast<float>(size) * s); };
        GuiSetStyle(DEFAULT, TEXT_SIZE, F(10));            // raygui defaults: text 10, spinner buttons 24 wide,
        GuiSetStyle(SPINNER, SPIN_BUTTON_WIDTH, F(24));    // 2 px between spinner buttons and between toggles
        GuiSetStyle(SPINNER, SPIN_BUTTON_SPACING, F(2));
        GuiSetStyle(TOGGLE, GROUP_PADDING, F(2));

        // Zoom toward the mouse: pin the point under the mouse, then scale around it.
        float wheel = GetMouseWheelMove();
        if (wheel != 0.0f) {
            cam.target = GetScreenToWorld2D(GetMousePosition(), cam);
            cam.offset = GetMousePosition();
            cam.zoom *= std::exp(0.2f * wheel);
        }
        if (IsMouseButtonDown(MOUSE_BUTTON_RIGHT)) {
            Vector2 d = GetMouseDelta();    // pixels moved since last frame
            cam.target.x -= d.x / cam.zoom;
            cam.target.y -= d.y / cam.zoom;
        }
        bool typing = editM || editP || editTT;
        if (!typing && IsKeyPressed(KEY_R)) cam = airfoilView;
        if (!typing && IsKeyPressed(KEY_F)) cam = farView;

        float px = 1.0f / cam.zoom;         // one pixel, in chord units, so lines stay 1 px at any zoom

        BeginDrawing();
        ClearBackground(RAYWHITE);
        BeginMode2D(cam);

        if (showComputed) {
            // flipped cells filled in red, under the lines. Both triangle orders are drawn because raylib
            // only fills counter-clockwise triangles, and a flipped cell turns the other way.
            for (const Cell& c : badCells) {
                auto v = [&](int i, int j) { return Vector2{static_cast<float>(grid.at(i, j).x), static_cast<float>(-grid.at(i, j).y)}; };
                Vector2 a = v(c.i, c.j), b = v(c.i + 1, c.j), cc = v(c.i + 1, c.j + 1), d = v(c.i, c.j + 1);
                DrawTriangle(a, b, cc, Fade(RED, 0.5f)); DrawTriangle(a, cc, b, Fade(RED, 0.5f));
                DrawTriangle(a, cc, d, Fade(RED, 0.5f)); DrawTriangle(a, d, cc, Fade(RED, 0.5f));
            }
            drawGrid(grid, BLUE, SKYBLUE, RED, px);
        }
        if (showAirfoil)                    // exact airfoil on top, thinner than the wall, so the wall shows around it
            DrawLineStrip(airfoil.data(), static_cast<int>(airfoil.size()), DARKGREEN);

        EndMode2D();

        // Buttons are drawn after EndMode2D, in screen pixels. They return their new state through the pointer.
        GuiToggle(Rectangle{right - S(230), S(15), S(104), S(24)}, "Computed", &showComputed);
        GuiToggle(Rectangle{right - S(124), S(15), S(104), S(24)}, "Profil NACA", &showAirfoil);

        int clicked = level;
        GuiToggleGroup(Rectangle{right - S(230), S(45), S(40), S(24)}, "8;16;32;64;128", &clicked);   // one item per level, side by side

        // NACA digits. A spinner returns 1 when its number box is clicked: that switches typing on or off.
        DrawText("NACA", static_cast<int>(right - S(275)), F(81), F(10), DARKGRAY);
        if (GuiSpinner(Rectangle{right - S(230), S(75), S(66), S(24)}, nullptr, &m, 0, 9, editM)) editM = !editM;
        if (GuiSpinner(Rectangle{right - S(158), S(75), S(66), S(24)}, nullptr, &p, 0, 9, editP)) editP = !editP;
        if (GuiSpinner(Rectangle{right - S(86), S(75), S(66), S(24)}, nullptr, &tt, 1, 40, editTT)) editTT = !editTT;

        bool valid = (m == 0) == (p == 0);  // camber M needs its position P, and P means nothing without M
        bool changed = clicked != level || m != digitM || p != digitP || tt != digitTT;
        if (changed && valid) {
            level = clicked;
            digitM = m; digitP = p; digitTT = tt;
            load();
        } else if (clicked != level) {          // nc can change even while the digits are invalid
            level = clicked;
            load();
        }

        DrawText(TextFormat("O-grid NACA %i%i%02i, nc = %i", digitM, digitP, digitTT, nc), F(20), F(20), F(20), DARKGRAY);
        DrawText("maillage", F(20), F(48), F(10), BLUE);
        DrawText("profil NACA exact", F(80), F(48), F(10), DARKGREEN);
        if (!valid)
            DrawText(TextFormat("NACA %i%i%02i: M and P must be both 0 or both > 0", m, p, tt), F(20), F(64), F(10), RED);

        // Validity panel: green when fine, red when not
        int bad = static_cast<int>(badCells.size());
        DrawText(TextFormat("cellules retournees : %i / %i", bad, nc * nc), F(20), F(84), F(10), bad == 0 ? DARKGREEN : RED);
        if (digitM == 0)
            DrawText(TextFormat("ecart paroi / NACA : %.1e", wallError), F(20), F(98), F(10), wallError < 1e-9 ? DARKGREEN : RED);
        else
            DrawText("ecart paroi / NACA : symetrique seulement", F(20), F(98), F(10), GRAY);
        DrawText(TextFormat("couture i=0 / i=nc : %.1e", seamGap), F(20), F(112), F(10), seamGap < 1e-9 ? DARKGREEN : RED);
        DrawText("wheel: zoom   right-drag: move   R: airfoil   F: far field", F(20), static_cast<int>(bottom - S(30)), F(10), GRAY);
        EndDrawing();
    }

    CloseWindow();
    return 0;
}
