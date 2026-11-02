#include "nacawing/3Dwing.hpp"
#include <cmath>
#include <vector>
#include <algorithm>
#include <numbers>

namespace nacawing {
    //For post traitement
    double area(const WingParams& w)            //S is the full wing area (both halves)
        {
            double lambda = taperRatio(w);
            return 0.5*w.b*w.c_r*(1.0+lambda);
        }
        
    double aspectRatio(const WingParams& w)     //AR the aspect ratio
        {
            double S = area(w);
            return w.b*w.b/S;
        }
        
    double meanAeroChord(const WingParams& w)  //c̄ the mean aerodynamic chord (corde aérodynamique moyenne), the reference length for the pitching moment.
        {
            double lambda = taperRatio(w);
            return 2.0/3.0*w.c_r*(1.0+lambda+lambda*lambda)/(1.0+lambda);
        }
    
    //Making the wing
    double taperRatio(const WingParams& w)      //lambda ; taper ratio: 1 = rectangular, 0 = pointed tip
        {
            return w.c_t / w.c_r;
        }   

    double chordAt(const WingParams& w, double z)
        {
            double lambda = taperRatio(w);
            double eta = z / (w.b / 2.0);
            return w.c_r * (1.0-(1.0-lambda)*eta);
        }

    std::vector<Point2D> sectionLoop(double m, double p, double t, int n)
        {
            // Point 0 is the TE upper passes to LE then lower to its TE
            std::vector<Point2D> upper = upperSurface(m, p, t, n);
            std::vector<Point2D> lower = lowerSurface(m, p, t, n);
            std::reverse(upper.begin(), upper.end());

            upper.insert(upper.end(), lower.begin() + 1, lower.end());
            return upper;
        }

    std::vector<Point3D> sectionAt(const WingParams& w, const std::vector<Point2D>& loop, double z)
        {

            auto AngleToRad = [](double angle) {
                return angle * std::numbers::pi / 180.0;};

            double eta = z / (w.b / 2.0);
            double c = chordAt(w, z);
            double epsilon = AngleToRad(w.twistTip * eta);
            double x_LE = z * std::tan(AngleToRad(w.sweepLE));
            double x_p = c * 0.25;

            std::vector<Point3D> section;
            section.reserve(loop.size());          // allocate 2n−1 slots up front
            double di = std::tan(AngleToRad(w.dihedral));
            for (const auto& q : loop) {                // q.x, q.y = unit loop point
                // 4. scale
                double x  = q.x * c;
                double y  = q.y * c;
                // 5. twist
                double xr = x_p + (x-x_p)*std::cos(epsilon) + y*std::sin(epsilon);
                double yr = -(x-x_p)*std::sin(epsilon) + y*std::cos(epsilon);
                // 6. translate
                section.push_back(Point3D{x_LE  + xr, yr + z*di, z });      
                }

                return section;
        }

     std::vector<std::vector<Point3D>>  wingSurface(const WingParams& w, const std::vector<Point2D>& loop, int Ns,
                                                    SpanSpacing spacing)
         {
            // Surface[i][j]: one row per loop point i, one column per station j (see 3Dwing.hpp).
            std::vector<std::vector<Point3D>> Surface(loop.size(), std::vector<Point3D>(Ns));

            for (int j = 0; j < Ns ; j++) {
                double eta = static_cast<double>(j) / (Ns - 1);         // 0 at the root, 1 at the tip
                if (spacing == SpanSpacing::SineTip)
                    eta = std::sin(0.5 * std::numbers::pi * eta);       // sin(pi j / (2(Ns-1)))
                double z = 0.5 * w.b * eta;

                std::vector<Point3D> section = sectionAt(w, loop, z);   // the whole loop at this z
                for (size_t i = 0; i < loop.size(); i++)
                    Surface[i][j] = section[i];                         // goes in column j
            }
            return Surface;
         }
        
    } // namespace nacawing
