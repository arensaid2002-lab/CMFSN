// Entry point of the NACA wing generator.
//
// For now it only proves the toolchain works. The next steps will add:
//   1. the NACA 4-digit airfoil equations (2D profile),
//   2. a basic primitive (a square) to test the geometry/output pipeline,
//   3. extrusion of the profile into a 3D wing,
//   4. export to a file a viewer can open.


#include "nacawing/naca.hpp"

#include <cmath>

namespace nacawing {

    double y_t(double t, double x)
        {
            return 5 * t * (0.2969 * std::sqrt(x) - 0.1260 * x - 0.3516 * std::pow(x, 2) + 0.2843 * std::pow(x, 3) - 0.1015 * std::pow(x, 4));
        }

    } // namespace nacawing
