#pragma once
#include "nacawing/NACASurface.hpp"


namespace nacawing {

// z axis representes the extrusion, one wing goes from (0 ≤ z ≤ b/2), mirror for the second
// (z = 0) to the tip (z = b/2)

    struct WingParams {
        double b;                   // [m] full tip-to-tip span; we build z from 0 to b/2
        double c_r;                 // [m] chord at the root, z = 0
        double c_t;                 // [m] chord at the tip, z = b/2

        double sweepLE;             // [deg]  \Lambda_{LE}
        double twistTip;            // [deg]  \varepsilon_t,
        double dihedral;            // [deg]  phi (not Gamma: in the VLM, Gamma is the circulation)
    };
    //For post traitement
    double taperRatio (const WingParams& w);

    double area(const WingParams& w);

    double aspectRatio(const WingParams& w);

    double meanAeroChord(const WingParams& w);
    
    //Making the wing
    double chordAt(const WingParams& w, double z);

    struct Point3D {
        double x;
        double y;
        double z;
    };

    std::vector<Point2D> sectionLoop(double m, double p, double t, int n);

    std::vector<Point3D> sectionAt(const WingParams& w, const std::vector<Point2D>& loop, double z);

    // How the Ns stations are spread along the span (see wingSurface).
    //   Uniform: z_j = (b/2) j/(Ns-1), equal gaps.
    //   SineTip: z_j = (b/2) sin(pi j / (2(Ns-1))), gaps shrink toward the tip, where the VLM
    //            circulation drops to 0 and changes fastest. Use it for the VLM and the beam.
    enum class SpanSpacing { Uniform, SineTip };

    // Surface[i][j] = point i of the loop at station j (root j = 0, tip j = Ns-1).
    // j (span) varies fastest, so flattening row by row gives k = i*Ns + j.
    std::vector<std::vector<Point3D>>  wingSurface(const WingParams& w, const std::vector<Point2D>& loop, int Ns,
                                                   SpanSpacing spacing = SpanSpacing::Uniform);

    
} // namespace nacawing
