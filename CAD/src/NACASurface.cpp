
#include "nacawing/NACA.hpp"
#include "nacawing/NACASurface.hpp"
#include <cmath>
#include <vector>

// Creates a vector of points for the upper surface of the NACA airfoil based on the thickness distribution.

namespace nacawing {


    std::vector<Point2D> upperSurface(double m, double p, double t, int n)
        {

            std::vector<Point2D> upper;
            for (int i = 0; i < n; i++) {   
                double x = static_cast<double>(i)/(n-1); // x goes from 0 to 1
                double y_t = nacawing::y_t(t, x);
                double y_c= nacawing::y_c(m, p, x); // camber line y_c
                
                double theta = std::atan(nacawing::dyc_dx(m, p, x));
                double x_upper = x - y_t * std::sin(theta); // upper surface x coordinate
                double y_upper= y_c + y_t * std::cos(theta); // upper surface = camber line + half-thickness
                upper.push_back(Point2D{x_upper, y_upper});
            }
            
            return upper;
        }


    std::vector<Point2D> lowerSurface(double m, double p, double t, int n)
        {

            std::vector<Point2D> lower;
            for (int i = 0; i < n; i++) {   
                double x = static_cast<double>(i)/(n-1); // x goes from 0 to 1
                double y_t = nacawing::y_t(t, x);
                double y_c= nacawing::y_c(m, p, x); // camber line y_c
                
                double theta = std::atan(nacawing::dyc_dx(m, p, x));
                double x_lower = x + y_t * std::sin(theta);  // lower surface x coordinate
                double y_lower= y_c - y_t * std::cos(theta); // lower surface = camber line - half-thickness
                lower.push_back(Point2D{x_lower, y_lower});
            }
            
            return lower;
        }

    } // namespace nacawing

