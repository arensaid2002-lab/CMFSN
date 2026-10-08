#pragma once

namespace nacawing {

// Returns the result from NACA
    double y_t(double t, double x);
    double y_c(double m, double p, double x);
    double dyc_dx(double m, double p, double x);
    
} // namespace nacawing
