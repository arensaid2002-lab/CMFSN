#include "raylib.h"
#define RAYGUI_IMPLEMENTATION               // raygui is a single header: this line puts its code in this file
#include "raygui.h"
#include "nacawing/NACA.hpp"
#include "nacawing/NACASurface.hpp"
#include "nacawing/3Dwing.hpp"
#include <cmath>
#include <vector>

// Wing viewer (written with Claude): draws the half wing from nacawing::wingSurface() in 3D.
// Axes are the ones of 3Dwing.hpp: x along the chord, y up, z along the span (root z = 0, tip z = b/2).
// raylib is y-up too, so a Point3D goes straight into a raylib Vector3 (x -> x, y -> y, z -> z).
//
// Controls: drag with the left mouse button to turn around the wing, mouse wheel to zoom.
// The panel on the left changes the airfoil (NACA digits) and the wing (span, chords, sweep, twist, dihedral, stations).


// Point3D (double, from CAD) -> Vector3 (float, what raylib draws). sign = -1 mirrors the point to the left wing (z -> -z).
Vector3 toVec(const nacawing::Point3D& q, float sign)
{
    return Vector3{static_cast<float>(q.x), static_cast<float>(q.y), sign * static_cast<float>(q.z)};
}

// Draws one half wing: wing[j][i] = station j, point i of the loop (see wingSurface).
// sign = +1 for the right half, -1 for its mirror image (the left half).
void drawHalfWing(const std::vector<std::vector<nacawing::Point3D>>& wing, float sign)
{
    const size_t N = wing[0].size();        // points per section, 2n - 1 (same for every section)
    const size_t le = (N - 1) / 2;          // the leading edge is the middle point of the loop (Selig order)

    // Spanwise lines: point i of section j joined to point i of section j + 1.
    // Every 4 neighbours (i, j), (i+1, j), (i+1, j+1), (i, j+1) make one quad of the surface mesh.
    for (size_t j = 0; j + 1 < wing.size(); ++j) {
        for (size_t i = 0; i < N; ++i) {
            bool edge = (i == 0 || i == le || i == N - 1);   // trailing edge (upper and lower) and leading edge
            DrawLine3D(toVec(wing[j][i], sign), toVec(wing[j + 1][i], sign), edge ? DARKGRAY : LIGHTGRAY);
        }
    }

    // Sections: each one is a loop TE -> upper -> LE -> lower -> TE, drawn point to point.
    for (size_t j = 0; j < wing.size(); ++j) {
        Color c = (j == 0) ? RED : BLUE;    // root section in red, the others in blue
        for (size_t i = 0; i + 1 < N; ++i)
            DrawLine3D(toVec(wing[j][i], sign), toVec(wing[j][i + 1], sign), c);
        DrawLine3D(toVec(wing[j][N - 1], sign), toVec(wing[j][0], sign), GREEN);   // close the (slightly open) trailing edge
    }
}


int main()
{
    SetConfigFlags(FLAG_WINDOW_RESIZABLE | FLAG_MSAA_4X_HINT);   // resizable window, smoother lines (both before InitWindow)
    InitWindow(1100, 620, "Wing viewer");
    SetWindowMinSize(800, 450);
    SetTargetFPS(60);

    // Airfoil: NACA digits M P TT (same spinners as naca_viewer), NACA 2412 to start.
    int digitM = 2, digitP = 4, digitTT = 12;
    bool editM = false, editP = false, editTT = false;   // true while typing in that number box
    int n = 41;                                          // points per surface -> 2n - 1 = 81 points per section
    nacawing::NacaParams airfoil = nacawing::fromDigits(digitM, digitP, digitTT);

    // Wing: starts as the doc's example wing (b = 8 m, c_r = 1 m, c_t = 0.5 m, 20 deg sweep, -4 deg twist).
    // raygui sliders work with float, so the inputs are floats here and become doubles in NacaWingParams.
    float span = 8.0f;          // [m]   b, full tip-to-tip span
    float rootChord = 1.0f;     // [m]   c_r
    float tipChord = 0.5f;      // [m]   c_t
    float sweep = 20.0f;        // [deg] leading-edge sweep
    float twist = -4.0f;        // [deg] tip twist (negative = washout, nose down at the tip)
    float dihedral = 0.0f;      // [deg]
    float stations = 11.0f;     // Ns, number of sections from root to tip (rounded to an int below)
    bool mirror = false;        // also draw the left half

    // Orbit camera: it sits on a sphere around the middle of the wing.
    // yaw turns around the vertical y axis, pitch looks from above (> 0) or below (< 0), distance is the zoom.
    float yaw = -0.8f, pitch = 0.5f, distance = 7.0f;
    Camera3D camera{};
    camera.up = Vector3{0.0f, 1.0f, 0.0f};
    camera.fovy = 45.0f;                         // field of view, degrees
    camera.projection = CAMERA_PERSPECTIVE;

    while (!WindowShouldClose())
    {
        // UI scale, as in naca_viewer: 1 for an 800 x 450 window. S() for rectangles, F() for font sizes.
        float sx = static_cast<float>(GetScreenWidth()) / 800.0f, sy = static_cast<float>(GetScreenHeight()) / 450.0f;
        float s = sx < sy ? sx : sy;
        auto S = [s](float v) { return v * s; };
        auto F = [s](int size) { return static_cast<int>(static_cast<float>(size) * s); };
        GuiSetStyle(DEFAULT, TEXT_SIZE, F(10));
        GuiSetStyle(SPINNER, SPIN_BUTTON_WIDTH, F(24));
        GuiSetStyle(SPINNER, SPIN_BUTTON_SPACING, F(2));
        Rectangle panel{0.0f, 0.0f, S(250), static_cast<float>(GetScreenHeight())};

        // 1. Geometry. Rebuilt every frame: it is only Ns x (2n - 1) points, so it is cheap,
        //    and it means every slider change shows up at once without extra bookkeeping.
        nacawing::NacaWingParams w{.b = span, .c_r = rootChord, .c_t = tipChord,
                                   .sweepLE = sweep, .twistTip = twist, .dihedral = dihedral};
        int Ns = static_cast<int>(std::lround(stations));
        std::vector<nacawing::Point2D> loop = nacawing::sectionLoop(airfoil.m, airfoil.p, airfoil.t, n);
        std::vector<std::vector<nacawing::Point3D>> wing = nacawing::wingSurface(w, loop, Ns);

        // 2. Camera. Drag outside the panel to turn, wheel to zoom.
        if (IsMouseButtonDown(MOUSE_BUTTON_LEFT) && !CheckCollisionPointRec(GetMousePosition(), panel)) {
            Vector2 d = GetMouseDelta();             // how far the mouse moved since the last frame, pixels
            yaw -= d.x * 0.01f;
            pitch += d.y * 0.01f;
            if (pitch > 1.5f) pitch = 1.5f;           // stop just before straight up/down (the view would flip)
            if (pitch < -1.5f) pitch = -1.5f;
        }
        distance -= GetMouseWheelMove() * 0.5f;
        if (distance < 1.0f) distance = 1.0f;

        // Look at the middle of the half wing (or of the whole wing when mirrored).
        float midZ = mirror ? 0.0f : 0.25f * span;
        float midX = 0.5f * static_cast<float>(wing.back()[0].x);   // halfway between the root LE (x = 0) and the tip TE
        camera.target = Vector3{midX, 0.0f, midZ};
        camera.position = Vector3{camera.target.x + distance * std::cos(pitch) * std::sin(yaw),
                                  camera.target.y + distance * std::sin(pitch),
                                  camera.target.z + distance * std::cos(pitch) * std::cos(yaw)};

        BeginDrawing();
        ClearBackground(RAYWHITE);

        // 3. The 3D scene.
        BeginMode3D(camera);
        DrawGrid(40, 0.5f);                          // floor grid in the x-z plane (y = 0), 0.5 m squares
        DrawLine3D(Vector3{0, 0, 0}, Vector3{1, 0, 0}, RED);     // x axis, 1 m long (chord direction)
        DrawLine3D(Vector3{0, 0, 0}, Vector3{0, 1, 0}, GREEN);   // y axis (up)
        DrawLine3D(Vector3{0, 0, 0}, Vector3{0, 0, 1}, BLUE);    // z axis (span)
        drawHalfWing(wing, 1.0f);
        if (mirror) drawHalfWing(wing, -1.0f);
        EndMode3D();

        // 4. The panel with the inputs.
        DrawRectangleRec(panel, Fade(LIGHTGRAY, 0.5f));
        DrawText("Wing viewer", F(12), F(10), F(20), DARKGRAY);

        DrawText(TextFormat("NACA %i%i%02i", digitM, digitP, digitTT), F(12), F(40), F(10), DARKGRAY);
        if (GuiSpinner(Rectangle{S(12), S(54), S(70), S(22)}, nullptr, &digitM, 0, 9, editM)) editM = !editM;
        if (GuiSpinner(Rectangle{S(88), S(54), S(70), S(22)}, nullptr, &digitP, 0, 9, editP)) editP = !editP;
        if (GuiSpinner(Rectangle{S(164), S(54), S(76), S(22)}, nullptr, &digitTT, 1, 40, editTT)) editTT = !editTT;
        bool valid = (digitM == 0) == (digitP == 0);   // camber M needs its position P, and P means nothing without M
        if (valid) airfoil = nacawing::fromDigits(digitM, digitP, digitTT);   // invalid digits: keep the last good airfoil
        else DrawText("M and P must be both 0 or both > 0", F(12), F(80), F(10), RED);

        // Sliders: label on the left, value on the right, the bar in between.
        float y = S(100);
        GuiSliderBar(Rectangle{S(70), y, S(110), S(14)}, "span b", TextFormat("%.1f m", span), &span, 1.0f, 20.0f);              y += S(24);
        GuiSliderBar(Rectangle{S(70), y, S(110), S(14)}, "root c_r", TextFormat("%.2f m", rootChord), &rootChord, 0.2f, 3.0f);  y += S(24);
        GuiSliderBar(Rectangle{S(70), y, S(110), S(14)}, "tip c_t", TextFormat("%.2f m", tipChord), &tipChord, 0.05f, 3.0f);    y += S(24);
        GuiSliderBar(Rectangle{S(70), y, S(110), S(14)}, "sweep LE", TextFormat("%.0f deg", sweep), &sweep, -20.0f, 60.0f);     y += S(24);
        GuiSliderBar(Rectangle{S(70), y, S(110), S(14)}, "tip twist", TextFormat("%.1f deg", twist), &twist, -10.0f, 10.0f);    y += S(24);
        GuiSliderBar(Rectangle{S(70), y, S(110), S(14)}, "dihedral", TextFormat("%.1f deg", dihedral), &dihedral, -10.0f, 15.0f); y += S(24);
        GuiSliderBar(Rectangle{S(70), y, S(110), S(14)}, "stations", TextFormat("%i", Ns), &stations, 2.0f, 41.0f);             y += S(28);
        GuiCheckBox(Rectangle{S(12), y, S(14), S(14)}, "mirror left half", &mirror);                                          y += S(30);

        // The wing numbers (reference values for the aero coefficients later).
        DrawText(TextFormat("taper ratio  %.3f", nacawing::taperRatio(w)), F(12), static_cast<int>(y), F(10), DARKGRAY);           y += S(16);
        DrawText(TextFormat("area S       %.3f m2", nacawing::area(w)), F(12), static_cast<int>(y), F(10), DARKGRAY);               y += S(16);
        DrawText(TextFormat("aspect ratio %.2f", nacawing::aspectRatio(w)), F(12), static_cast<int>(y), F(10), DARKGRAY);           y += S(16);
        DrawText(TextFormat("MAC          %.3f m", nacawing::meanAeroChord(w)), F(12), static_cast<int>(y), F(10), DARKGRAY);       y += S(16);
        DrawText(TextFormat("%i sections x %i points", Ns, static_cast<int>(loop.size())), F(12), static_cast<int>(y), F(10), DARKGRAY);

        DrawText("drag: turn   wheel: zoom", F(12), GetScreenHeight() - F(20), F(10), GRAY);
        EndDrawing();
    }

    CloseWindow();
    return 0;
}
