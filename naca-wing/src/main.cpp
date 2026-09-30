// Entry point of the NACA wing generator.
//
// For now it only proves the toolchain works. The next steps will add:
//   1. the NACA 4-digit airfoil equations (2D profile),
//   2. a basic primitive (a square) to test the geometry/output pipeline,
//   3. extrusion of the profile into a 3D wing,
//   4. export to a file a viewer can open.

#include <iostream>
#include "nacawing/upperSurface.hpp"
#include "nacawing/version.hpp"
#include "nacawing/naca.hpp"

int main()
{
    std::cout << "NACA wing generator v" << nacawing::version() << '\n';
    std::cout << "Toolchain OK: C++" << __cplusplus / 100 % 100 << " build is running.\n";
    std::cout << "NACA 4-digit airfoil thickness at x=0.5 for t=0.12 is: " << nacawing::y_t(0.12, 0.5) << '\n';
    std::vector<nacawing::Point2D> list = nacawing::upperSurface(0.12, 11);
        for (size_t i = 0; i < list.size(); ++i) {
            std::cout << "Point " << i << ": (" << list[i].x << ", " << list[i].y << ")\n";
        }
    return 0;
}
