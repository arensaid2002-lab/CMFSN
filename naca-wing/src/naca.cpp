// Entry point of the NACA wing generator.
//
// For now it only proves the toolchain works. The next steps will add:
//   1. the NACA 4-digit airfoil equations (2D profile),
//   2. a basic primitive (a square) to test the geometry/output pipeline,
//   3. extrusion of the profile into a 3D wing,
//   4. export to a file a viewer can open.


#include "nacawing/NACA.hpp"

#include <cmath>

namespace nacawing {

    double y_t(double t, double x)
        {
            return 5 * t * (0.2969 * std::sqrt(x) - 0.1260 * x - 0.3516 * x * x + 0.2843 * x * x * x - 0.1015 * x * x *x * x);
        }

    double y_c(double m, double p, double x)
        {
            if (m == 0)
            {
                return 0.0;
            }

            if (x < p)
            {
                return m / (p * p) * (2 * p * x - x * x);
            }
            else
            {
                return m / ((1 - p) * (1 - p)) * (1 - 2 * p + 2 * p * x - x * x);
            }
        }


        double dyc_dx(double m, double p, double x)
        {

            if (x < p)
            {
                return 2*m / (p * p) * (p - x);
            }
            else
            {
                return 2*m / ((1 - p) * (1 - p)) * (p - x);
            }
        }

    } // namespace nacawing
