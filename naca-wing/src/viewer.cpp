#include "raylib.h"
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
    


    while (!WindowShouldClose())            // true when you click X or press Esc
    {
        BeginDrawing();                     // start the frame
        ClearBackground(RAYWHITE);          // erase last frame with a white-ish colour
        
        for (size_t i = 0; i < upper.size(); ++i) {
        
            DrawLine(toPixelX(upper[i].x), toPixelY(upper[i].y), toPixelX(lower[i].x), toPixelY(lower[i].y), BLACK);
        
        }

        DrawText("NACA viewer", 20, 20, 20, DARKGRAY);   // text, x, y (pixels from top-left), size, colour
        EndDrawing();                       // show the frame on screen
    }

    CloseWindow();                          // clean up before the program ends
    return 0;
}