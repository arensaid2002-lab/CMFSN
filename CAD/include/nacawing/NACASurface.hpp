#pragma once
#include <vector>
#include "nacawing/NACA.hpp"   // NacaParams {m, p, t}

namespace nacawing {

// builds the results from NACA equations into a vector of points for the upper surface of the airfoil

    struct Point2D {
        double x;
        double y;
    };
    double cosineX(int i, int n, double xEnd = 1.0);
    // naca = the airfoil shape (m, p, t, fractions of the chord, see fromDigits).
    // n = number of points per surface, xEnd = last x (1 = TE). n and xEnd are the sampling, not the shape,
    // so they stay separate arguments.
    std::vector<Point2D> upperSurface(const NacaParams& naca, int n, double xEnd = 1.0);
    std::vector<Point2D> lowerSurface(const NacaParams& naca, int n, double xEnd = 1.0);

} // namespace nacawing
