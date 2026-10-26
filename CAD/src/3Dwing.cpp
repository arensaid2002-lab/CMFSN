#include "nacawing/3Dwing.hpp"
#include <cmath>

namespace nacawing {

    double taperRatio(const NacaWingParams& w)      //lambda ; taper ratio: 1 = rectangular, 0 = pointed tip
        {
            return w.c_t / w.c_r;
        }     

    double area(const NacaWingParams& w)            //S is the full wing area (both halves)
        {
            double lambda = taperRatio(w);
            return 0.5*w.b*w.c_r*(1+lambda);
        }
        
    double aspectRatio(const NacaWingParams& w)     //AR the aspect ratio
        {
            double S = area(w);
            return w.b*w.b/S;
        }
        
    double meanAeroChord(const NacaWingParams& w)  //c̄ the mean aerodynamic chord (corde aérodynamique moyenne), the reference length for the pitching moment.
        {
            double lambda = taperRatio(w);
            return 2.0/3.0*w.c_r*(1+lambda+lambda*lambda)/(1+lambda);
        }
    
    double chordAt(const NacaWingParams& w, double z)
        {
            double lambda = taperRatio(w);
            double eta = z / (w.b / 2.0);
            return w.c_r * (1-(1-lambda)*eta);
        }

    } // namespace nacawing
