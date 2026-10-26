#pragma once



namespace nacawing {

// z axis representes the extrusion, one wing goes from (0 ≤ z ≤ b/2), mirror for the second
// (z = 0) to the tip (z = b/2)

    struct NacaWingParams {
        double b;                   // [m] full tip-to-tip span; we build z from 0 to b/2
        double c_r;                 // [m] chord at the root, z = 0
        double c_t;                 // [m] chord at the tip, z = b/2

        double sweepLE;             // [deg] 
        double twistTip;            // [deg] 
        double dihedral;            // [deg] 
    };

    double taperRatio (const NacaWingParams& w);

    double area(const NacaWingParams& w);

    double aspectRatio(const NacaWingParams& w);

    double meanAeroChord(const NacaWingParams& w);

    double chordAt(const NacaWingParams& w, double z);

    
} // namespace nacawing
