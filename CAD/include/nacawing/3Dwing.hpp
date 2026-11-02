#pragma once
#include "nacawing/NACASurface.hpp"


namespace nacawing {

// Axes (right-handed): x along the chord (LE -> TE), y along the span (extrusion), z up.
// One wing goes from the root (y = 0) to the tip (y = b/2); mirror for the second.
// The 2D profile (x, y_profile) is placed in the x-z plane: y_profile becomes z.

    struct WingParams {
        double b;                   // [m] full tip-to-tip span; we build y from 0 to b/2
        double c_r;                 // [m] chord at the root, y = 0
        double c_t;                 // [m] chord at the tip, y = b/2

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
    double chordAt(const WingParams& w, double y);

    struct Point3D {
        double x;
        double y;
        double z;
    };

    std::vector<Point2D> sectionLoop(double m, double p, double t, int n);

    std::vector<Point3D> sectionAt(const WingParams& w, const std::vector<Point2D>& loop, double y);

    // How the Ns stations are spread along the span (see wingSurface).
    //   Uniform: y_j = (b/2) j/(Ns-1), equal gaps.
    //   SineTip: y_j = (b/2) sin(pi j / (2(Ns-1))), gaps shrink toward the tip, where the VLM
    //            circulation drops to 0 and changes fastest. Use it for the VLM and the beam.
    enum class SpanSpacing { Uniform, SineTip };

    // Surface[i][j] = point i of the loop at station j (root j = 0, tip j = Ns-1).
    // j (span) varies fastest, so flattening row by row gives k = i*Ns + j.
    std::vector<std::vector<Point3D>>  wingSurface(const WingParams& w, const std::vector<Point2D>& loop, int Ns,
                                                   SpanSpacing spacing = SpanSpacing::Uniform);

    
} // namespace nacawing
