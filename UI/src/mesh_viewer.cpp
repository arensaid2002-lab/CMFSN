// Mesh viewer (written with Claude): draws the O-grid from mesh::generateOGrid,
// and on top of it the reference grid of the same size (tests/data/x.9 ... x.5, NACA 0012 only).
//   buttons            Computed / Reference: show or hide each grid
//                      8 ... 128: number of cells nc (regenerates the grid, loads x.9 ... x.5)
//                      NACA M P TT: the airfoil digits, with - / + or click the number to type it
//   mouse wheel        zoom toward the mouse
//   right-click drag   move the view
//   R                  back to the airfoil view
//   F                  show the whole grid, out to the far field

#include <cmath>
#include <string>
#include "raylib.h"
#define RAYGUI_IMPLEMENTATION               // raygui is a single header: this line puts its code in this file
#include "raygui.h"
#include "OGrid.h"
#include "Plot3D.h"

int main()
{
    SetConfigFlags(FLAG_WINDOW_RESIZABLE);  // must come before InitWindow: lets you drag the window edges
    InitWindow(800, 450, "O-grid viewer");  // a second program, so a second window next to naca_viewer
    SetWindowMinSize(640, 360);             // small enough to shrink, big enough for the buttons
    SetTargetFPS(60);
    int screenW = GetScreenWidth(), screenH = GetScreenHeight();   // last known size, to measure a resize

    int digitM = 0, digitP = 0, digitTT = 12;   // NACA 0012, the airfoil of the reference grids
    int m = digitM, p = digitP, tt = digitTT;   // what the spinners show: can be invalid for a moment (2012 on the way to 2412)

    // Reference file x.(9 - level) has nc = 8 * 2^level cells: x.9 -> 8, x.8 -> 16, ... x.5 -> 128.
    const int levels[] = {8, 16, 32, 64, 128};
    int level = 3;                          // start at nc = 64
    int nc = 0;
    mesh::StructuredMesh2D grid;            // computed by generateOGrid
    mesh::StructuredMesh2D reference;       // read from the PLOT3D file
    std::string referenceFile;

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
        grid = mesh::generateOGrid(nacawing::fromDigits(digitM, digitP, digitTT), nc);
        referenceFile = "x." + std::to_string(9 - level);
        reference = mesh::readPlot3D(std::string(REF_DIR) + "/" + referenceFile);  // empty if the file is missing

        // the far field is ~150 chords away (more for thick or cambered airfoils): fit its radius in the window
        nacawing::Point2D far = grid.at(0, grid.nj - 1);
        farRadius = static_cast<float>(std::hypot(far.x - 0.5, far.y));
        farView.zoom = 0.45f * static_cast<float>(GetScreenHeight()) / farRadius;
    };
    load();

    bool showComputed = true;
    bool showReference = true;
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
        bool is0012 = (digitM == 0 && digitP == 0 && digitTT == 12);   // the only airfoil with reference grids

        BeginDrawing();
        ClearBackground(RAYWHITE);
        BeginMode2D(cam);

        if (showReference && is0012 && !reference.nodes.empty())
            drawGrid(reference, ORANGE, Fade(ORANGE, 0.5f), MAROON, 2.0f * px);   // under, and wider, so both stay visible
        if (showComputed)
            drawGrid(grid, BLUE, SKYBLUE, RED, px);

        EndMode2D();

        // Buttons are drawn after EndMode2D, in screen pixels. They return their new state through the pointer.
        GuiToggle(Rectangle{right - S(230), S(15), S(104), S(24)}, "Computed", &showComputed);
        GuiToggle(Rectangle{right - S(124), S(15), S(104), S(24)}, "Reference", &showReference);

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
        DrawText("computed", F(20), F(48), F(10), BLUE);
        if (!is0012)
            DrawText("reference grids exist for NACA 0012 only", F(80), F(48), F(10), GRAY);
        else if (reference.nodes.empty())
            DrawText(TextFormat("reference %s not found in %s", referenceFile.c_str(), REF_DIR), F(80), F(48), F(10), RED);
        else
            DrawText(TextFormat("reference %s", referenceFile.c_str()), F(80), F(48), F(10), ORANGE);
        if (!valid)
            DrawText(TextFormat("NACA %i%i%02i: M and P must be both 0 or both > 0", m, p, tt), F(20), F(64), F(10), RED);
        DrawText("wheel: zoom   right-drag: move   R: airfoil   F: far field", F(20), static_cast<int>(bottom - S(30)), F(10), GRAY);
        EndDrawing();
    }

    CloseWindow();
    return 0;
}
