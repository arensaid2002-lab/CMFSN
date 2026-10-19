#include "raylib.h"
#define RAYGUI_IMPLEMENTATION               // raygui is a single header: this line puts its code in this file
#include "raygui.h"
#include "nacawing/NACA.hpp"
#include "nacawing/NACASurface.hpp"

int toPixelX(double x) { return static_cast<int>(50 + x * 700); }    // 50 px margin, chord = 700 px
int toPixelY(double y) { return static_cast<int>(225 - y * 700); }   // 225 = window middle, minus flips y


int main()
{
    
    InitWindow(800, 450, "NACA viewer");   // open a window: width, height in pixels, title
    SetTargetFPS(60);                       // run the loop at most 60 times per second

    double m = 0.02; // maximum camber
    double p = 0.4;  // location of maximum camber
    double t = 0.12; // maximum thickness
    int n = 101;     // number of points to generate
    
    std::vector<nacawing::Point2D> upper = nacawing::upperSurface(m, p, t, n); // NACA 
    std::vector<nacawing::Point2D> lower = nacawing::lowerSurface(m, p, t, n);

    // NACA buttons (added with Claude): the digits M P TT shown by the spinners, NACA 2412 to match m, p, t above.
    // They can be invalid for a moment (2012 on the way to 2412), so they are kept apart from m, p, t.
    int digitM = 2, digitP = 4, digitTT = 12;
    bool editM = false, editP = false, editTT = false;   // true while typing in that number box
    


    while (!WindowShouldClose())            // true when you click X or press Esc
    {
        BeginDrawing();                     // start the frame
        ClearBackground(RAYWHITE);          // erase last frame with a white-ish colour
        
        // join each point to the next one on the same surface (i + 1 < size so i+1 stays in range)
        for (size_t i = 0; i + 1 < upper.size(); ++i) {
            DrawLine(toPixelX(upper[i].x), toPixelY(upper[i].y), toPixelX(upper[i+1].x), toPixelY(upper[i+1].y), RED);
            DrawLine(toPixelX(lower[i].x), toPixelY(lower[i].y), toPixelX(lower[i+1].x), toPixelY(lower[i+1].y), BLUE);
            DrawCircle(toPixelX(upper[i].x), toPixelY(upper[i].y), 3, RED);
            DrawCircle(toPixelX(lower[i].x), toPixelY(lower[i].y), 3, BLUE);
        }
        DrawLine(toPixelX(lower[lower.size() - 1].x), toPixelY(lower[lower.size() - 1].y), toPixelX(upper[upper.size() - 1].x), toPixelY(upper[upper.size() - 1].y), GREEN); // close the trailing edge


        DrawText("NACA viewer", 20, 20, 20, DARKGRAY);   // text, x, y (pixels from top-left), size, colour

        // A spinner returns 1 when its number box is clicked: that switches typing on or off.
        DrawText(TextFormat("NACA %i%i%02i", digitM, digitP, digitTT), 500, 21, 10, DARKGRAY);
        if (GuiSpinner(Rectangle{570, 15, 66, 24}, nullptr, &digitM, 0, 9, editM)) editM = !editM;
        if (GuiSpinner(Rectangle{642, 15, 66, 24}, nullptr, &digitP, 0, 9, editP)) editP = !editP;
        if (GuiSpinner(Rectangle{714, 15, 66, 24}, nullptr, &digitTT, 1, 40, editTT)) editTT = !editTT;

        nacawing::NacaParams wanted = nacawing::fromDigits(digitM, digitP, digitTT);
        bool valid = (digitM == 0) == (digitP == 0);   // camber M needs its position P, and P means nothing without M
        if (valid && (wanted.m != m || wanted.p != p || wanted.t != t)) {
            m = wanted.m; p = wanted.p; t = wanted.t;
            upper = nacawing::upperSurface(m, p, t, n);
            lower = nacawing::lowerSurface(m, p, t, n);
        }
        if (!valid)
            DrawText("M and P must be both 0 or both > 0", 570, 45, 10, RED);
        EndDrawing();                       // show the frame on screen
    }

    CloseWindow();                          // clean up before the program ends
    return 0;
}