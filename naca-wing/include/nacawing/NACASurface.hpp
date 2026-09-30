#pragma once
#include <vector>

namespace nacawing {

// builds the results from NACA equations into a vector of points for the upper surface of the airfoil

    struct Point2D {
        double x;
        double y;
    };

    std::vector<Point2D> upperSurface(double t, int n);
    std::vector<Point2D> lowerSurface(double t, int n);

} // namespace nacawing
