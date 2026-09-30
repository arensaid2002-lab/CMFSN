// Entry point of the NACA wing generator.
//
// For now it only proves the toolchain works. The next steps will add:
//   1. the NACA 4-digit airfoil equations (2D profile),
//   2. a basic primitive (a square) to test the geometry/output pipeline,
//   3. extrusion of the profile into a 3D wing,
//   4. export to a file a viewer can open.

#include <iostream>

#include "nacawing/version.hpp"

int main()
{
    std::cout << "NACA wing generator v" << nacawing::version() << '\n';
    std::cout << "Toolchain OK: C++" << __cplusplus / 100 % 100 << " build is running.\n";
    return 0;
}
