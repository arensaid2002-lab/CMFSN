#pragma once

namespace nacawing {

// Returns the result from NACA
    double y_t(double t, double x);
    double y_c(double m, double p, double x);
    double dyc_dx(double m, double p, double x);
    struct NacaParams {
        double m; // maximum camber
        double p; // location of maximum camber
        double t; // maximum thickness
    };
    NacaParams fromDigits(int M, int P, int TT);   
} // namespace nacawing
