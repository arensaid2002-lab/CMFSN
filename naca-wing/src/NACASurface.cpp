
#include "nacawing/naca.hpp"
#include "nacawing/NACASurface.hpp"
#include <vector>

// Creates a vector of points for the upper surface of the NACA airfoil based on the thickness distribution.

namespace nacawing {


    std::vector<Point2D> upperSurface(double t, int n)
        {

            std::vector<Point2D> upper;
            for (int i = 0; i < n; i++) {   
                double x_upper = static_cast<double>(i)/(n-1); // x goes from 0 to 1
                double y = nacawing::y_t(t, x_upper);
                double y_upper= 0.0 + y; // camber line y_c, 0 for now
                upper.push_back(Point2D{x_upper, y_upper});
            }
            
            return upper;
        }


    std::vector<Point2D> lowerSurface(double t, int n)
        {

            std::vector<Point2D> lower;
            for (int i = 0; i < n; i++) {   
                double x_lower = static_cast<double>(i)/(n-1); // x goes from 0 to 1
                double y = nacawing::y_t(t, x_lower);
                double y_lower= 0.0 - y; // camber line y_c, 0 for now
                lower.push_back(Point2D{x_lower, y_lower});
            }
            
            return lower;
        }

    } // namespace nacawing

