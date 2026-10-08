#include "raylib.h"

int main()
{
    InitWindow(800, 450, "NACA viewer");   // open a window: width, height in pixels, title
    SetTargetFPS(60);                       // run the loop at most 60 times per second

    while (!WindowShouldClose())            // true when you click X or press Esc
    {
        BeginDrawing();                     // start the frame
        ClearBackground(RAYWHITE);          // erase last frame with a white-ish colour
        DrawText("NACA viewer", 20, 20, 20, DARKGRAY);   // text, x, y (pixels from top-left), size, colour
        EndDrawing();                       // show the frame on screen
    }

    CloseWindow();                          // clean up before the program ends
    return 0;
}