#include "nacawing/3Dwing.hpp"
#include <cmath>
#include <vector>
#include <algorithm>
#include <numbers>

namespace nacawing {
    //For post traitement
    double taperRatio(const NacaWingParams& w)      //lambda ; taper ratio: 1 = rectangular, 0 = pointed tip
        {
            return w.c_t / w.c_r;
        }     

    double area(const NacaWingParams& w)            //S is the full wing area (both halves)
        {
            double lambda = taperRatio(w);
            return 0.5*w.b*w.c_r*(1.0+lambda);
        }
        
    double aspectRatio(const NacaWingParams& w)     //AR the aspect ratio
        {
            double S = area(w);
            return w.b*w.b/S;
        }
        
    double meanAeroChord(const NacaWingParams& w)  //c̄ the mean aerodynamic chord (corde aérodynamique moyenne), the reference length for the pitching moment.
        {
            double lambda = taperRatio(w);
            return 2.0/3.0*w.c_r*(1.0+lambda+lambda*lambda)/(1.0+lambda);
        }
    
    //Making the wing
    double chordAt(const NacaWingParams& w, double z)
        {
            double lambda = taperRatio(w);
            double eta = z / (w.b / 2.0);
            return w.c_r * (1.0-(1.0-lambda)*eta);
        }

    std::vector<Point2D> sectionLoop(double m, double p, double t, int n)
        {
            // Point 0 is the TE upeer passes to LE then lower to its TE
            std::vector<Point2D> upper = upperSurface(m, p, t, n);
            std::vector<Point2D> lower = lowerSurface(m, p, t, n);
            std::reverse(upper.begin(), upper.end());

            upper.insert(upper.end(), lower.begin() + 1, lower.end());
            return upper;
        }

    std::vector<Point3D> sectionAt(const NacaWingParams& w, const std::vector<Point2D>& loop, double z)
        {

            auto AngleToRad = [&](double angle) {
                return angle * std::numbers::pi / 180.0;};

            double eta = z / (w.b / 2.0);
            double c = chordAt(w, z);
            double epsilon = AngleToRad(w.twistTip * eta);
            double x_LE = z * std::tan(AngleToRad(w.sweepLE));
            double x_p = c * 0.25;

            std::vector<Point3D> section;
            section.reserve(loop.size());          // allocate 2n−1 slots up front
            for (const auto& q : loop) {                // q.x, q.y = unit loop point
                // 4. scale
                double x  = q.x * c;
                double y  = q.y * c;
                // 5. twist
                double xr = x_p + (x-x_p)*std::cos(epsilon) + y*std::sin(epsilon);
                double yr = -(x-x_p)*std::sin(epsilon) + y*std::cos(epsilon);
                // 6. translate
                section.push_back(Point3D{x_LE  + xr, yr + z*std::tan(AngleToRad(w.dihedral)), z });      
                }

                return section;
        }

     std::vector<std::vector<Point3D>>  wingSurface(const NacaWingParams& w, const std::vector<Point2D>& loop, int Ns)
         {
            std::vector<std::vector<Point3D>> Surface;
            Surface.reserve(Ns);
                        
            for (int i = 0; i < Ns ; i++) {
                double z = 0.5 * w.b * i / (Ns - 1);
                Surface.push_back(sectionAt(w, loop, z));
            }
             return Surface;
         }
        
    } // namespace nacawing
