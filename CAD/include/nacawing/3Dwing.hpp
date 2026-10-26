#pragma once
#include "nacawing/NACASurface.hpp"


namespace nacawing {

// z axis representes the extrusion, one wing goes from (0 ≤ z ≤ b/2), mirror for the second
// (z = 0) to the tip (z = b/2)

    struct NacaWingParams {
        double b;                   // [m] full tip-to-tip span; we build z from 0 to b/2
        double c_r;                 // [m] chord at the root, z = 0
        double c_t;                 // [m] chord at the tip, z = b/2

        double sweepLE;             // [deg]  \Lambda_{LE}
        double twistTip;            // [deg]  \varepsilon_t,
        double dihedral;            // [deg]  \Gamma$
    };

    double taperRatio (const NacaWingParams& w);

    double area(const NacaWingParams& w);

    double aspectRatio(const NacaWingParams& w);

    double meanAeroChord(const NacaWingParams& w);

    double chordAt(const NacaWingParams& w, double z);

    struct Point3D {
        double x;
        double y;
        double z;
    };

    std::vector<Point2D> sectionLoop(double m, double p, double t, int n);

    std::vector<Point3D> sectionAt(const NacaWingParams& w, const std::vector<Point2D>& loop, double z);

    std::vector<std::vector<Point3D>>  wingSurface(const NacaWingParams& w, const std::vector<Point2D>& loop, int Ns);

    
} // namespace nacawing
