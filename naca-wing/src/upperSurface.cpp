
#include "nacawing/naca.hpp"
#include "nacawing/upperSurface.hpp"
#include <vector>

// Creates a vector of points for the upper surface of the NACA airfoil based on the thickness distribution.

namespace nacawing {


    std::vector<Point2D> upperSurface(double t, int n)
        {

            std::vector<Point2D> wing;
            for (int i = 0; i < n; i++) {   
                double x = static_cast<double>(i)/(n-1); // x goes from 0 to 1
                double y = nacawing::y_t(t, x);
                wing.push_back(Point2D{x, y});
            }
            
            return wing;
        }

    } // namespace nacawing
