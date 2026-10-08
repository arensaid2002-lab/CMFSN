// Entry point of the NACA wing generator.
//
// For now it only proves the toolchain works. The next steps will add:
//   1. the NACA 4-digit airfoil equations (2D profile),
//   2. a basic primitive (a square) to test the geometry/output pipeline,
//   3. extrusion of the profile into a 3D wing,
//   4. export to a file a viewer can open.

#include <iostream>
#include "nacawing/NACASurface.hpp"
#include "nacawing/version.hpp"
#include "nacawing/NACA.hpp"
#include <fstream> // makes .csv file with x,y coordinates of upper and lower surfaces

int main()
{
    std::cout << "NACA wing generator v" << nacawing::version() << '\n';
    std::cout << "Toolchain OK: C++" << __cplusplus / 100 % 100 << " build is running.\n";


    // Test the NACA 4-digit airfoil equations
    std::cout << "NACA 4-digit airfoil thickness at x=0.5 for t=0.12 is: " << nacawing::y_t(0.12, 0.5) << '\n';
    std::cout << "NACA 4-digit airfoil y_c for m=0.02, p=0.4, x=0.4: " << nacawing::y_c(0.02, 0.4, 0.4) << '\n';
    std::cout << "NACA 4-digit airfoil y_c for m=0.0, p=0.0, x=0.5: " << nacawing::y_c(0.0, 0.0, 0.5) << '\n';


    // Creating the CSV file with the upper and lower surface coordinates of the NACA 2412 airfoil
    
    double m = 0.02; // maximum camber
    double p = 0.4;  // location of maximum camber
    double t = 0.12; // maximum thickness
    int n = 101;     // number of points to generate
    
    std::vector<nacawing::Point2D> upper = nacawing::upperSurface(m, p, t, n); // NACA 
    std::vector<nacawing::Point2D> lower = nacawing::lowerSurface(m, p, t, n);
    std::ofstream file("build/airfoil.csv");
    file << "x_upper,y_upper,x_lower,y_lower\n"; // header line of csv
        for (size_t i = 0; i < upper.size(); ++i) {
            
             file << upper[i].x << "," << upper[i].y << "," << lower[i].x << "," << lower[i].y << "\n"; // writes coordinates to csv file
            }

    return 0;
}
