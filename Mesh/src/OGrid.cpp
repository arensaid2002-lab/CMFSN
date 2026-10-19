#include "OGrid.h"
#include <complex>
#include <cmath>
#include <numbers>   //for pi
#include <algorithm> //for min

//Construction plan

    //INPUTS, naca.m, naca.p, naca.t, and nc (number of cells per direction, so ni = nj = nc + 1)
   

    // 1. constants: x_TE, delta (abs!), P, r_LE, z1, z2, zm, zeta1, zeta2

    // 2. the two maps, as lambdas so they see the constants:
    //    auto toCircle   = [&](cplx z)    { ... };
    //    auto fromCircle = [&](cplx zeta) { ... };

    // 3. zetaLE = toCircle(0), zc = (zetaLE + zeta1) / 2

    // 4. R1: sample the airfoil densely, sum |toCircle(p[k+1]) - toCircle(p[k])|, divide by 2*pi

    // 5. wall: for each i, theta_i = -2*pi*i/nc, find x by bisection, keep r_i0

    // 6. for each j, for each i: R_j, r_ij (eq. 14),
        // zeta = zc + r_ij * std::polar(1.0, theta_i), z = fromCircle(zeta)
        // r_ij avec l'éq. 14, en utilisant R[0], R[j], R[nc]
        // zeta = zeta_c + r_ij * std::polar(1.0, theta_i)
        // z = fromCircle(zeta)
        // mesh.at(i, j) = {z.real(), z.imag()};


namespace mesh {
    //  SOURCE : In Pursuit of Grid Convergence for Two-Dimensional Euler Solutions

    StructuredMesh2D generateOGrid(const nacawing::NacaParams& naca, int nc)
    {
        StructuredMesh2D mesh;
        double x_TE = 1.0089304115;                                     // trailing edge x coordinate constante for NACA 4 digits airfoils
        double delta = 2.0 * std::atan(std::abs(nacawing::dy_t(naca.t, x_TE)));     //  trailing-edge included angle

        double P = std::numbers::pi / (2*std::numbers::pi - delta);     //(éq. 9).   power of Karman–Trefftz mapping

        double rho_LE = 1.1019 * naca.t * naca.t;                      //éq. 7;  radius of curvatur
        std::complex<double> z_1 = std::complex<double>(x_TE,0.0);         //(éq. 10). physical plane z  x  iy
        std::complex<double> z_2 = std::complex<double>(0.5 * rho_LE, 0.0);

        //Given z_1 and z_2, we can compute zeta_1 and zeta_2 using the Karman–Trefftz mapping. The mapping is defined as follows:
        //zeta = (zeta_1 - zeta_2 * ((z - z_1) / (z - z_2))^P) / (1 - ((z - z_1) / (z - z_2))^P)
        //ζ1, ζ2 = zm ± P(z1 − z2)/2, avec zm = (z1 + z2)/2.                    //Equation to compute zeta_1 and zeta_2, but I don't understand it. So I just use the values from Claude's code.
        //std::complex<double> zeta_1 = std::complex<double>(0.77043505,0.0);        //Claude gives a equation to compute zeta_1 and zeta_2, but I don't understand it. So I just use the values from Claude's code.
        //std::complex<double> zeta_2 = std::complex<double>(0.24642903,0.0);              //
        //std::complex<double> zeta_LE = std::complex<double>(0.20139626,0.0);             //
        std::complex<double> zm = (z_1 + z_2) * 0.5;
        std::complex<double> zeta_1 = zm + P * (z_1 - z_2) * 0.5;
        std::complex<double> zeta_2 = zm - P * (z_1 - z_2) * 0.5;
        
    

        auto toCircle = [&](std::complex<double> z) {
            std::complex<double> W = std::pow((z - z_1) / (z - z_2), P);
            return (zeta_1 - zeta_2 * W) / (1.0 - W);
        };

        auto fromCircle = [&](std::complex<double> zeta) {
            std::complex<double> w = std::pow(((zeta - zeta_1) / (zeta - zeta_2)), 1.0 / P);
            return (z_1 - z_2 * w) / (1.0 - w);
        };

        
        std::complex<double> zeta_TE = zeta_1;
        std::complex<double> zeta_LE = toCircle(std::complex<double>(0.0, 0.0)); // Leading edge in zeta-plane
        std::complex<double> zeta_c = 0.5 * (zeta_LE + zeta_TE);

        //4.Calculating R1 from perimeter 
        //TOP
        int M = 2000;
        std::vector<nacawing::Point2D> upper = nacawing::upperSurface(naca.m, naca.p, naca.t, M , x_TE);
        double L = 0.0;
        std::complex<double> prev = zeta_LE;   // le LE, sans passer par toCircle

        for (int k = 1; k < M; ++k) {
            // calcule le point z_k du profil
            std::complex<double> z_k = std::complex<double>(upper[k].x, upper[k].y);
            std::complex<double> cur = toCircle(z_k);
            L += std::abs(cur - prev);
            prev = cur;
            }
        //BOTTOM
        std::vector<nacawing::Point2D> lower = nacawing::lowerSurface(naca.m, naca.p, naca.t, M , x_TE);
        prev = zeta_LE;   // le LE, sans passer par toCircle
        for (int k = 1; k < M; ++k) {
            // calcule le point z_k du profil
            std::complex<double> z_k = std::complex<double>(lower[k].x, lower[k].y);
            std::complex<double> cur = toCircle(z_k);
            L += std::abs(cur - prev);
            prev = cur;
            }

        double R1 = L / (2.0 * std::numbers::pi); // radius of the circle in zeta-plane



        //5. Calculating the wall points

        auto surfacePoint = [&](double x_m, bool isUpper) {
            double s = isUpper ? 1.0 : -1.0;
            double y_t = nacawing::y_t(naca.t, x_m);
            double y_c = nacawing::y_c(naca.m, naca.p, x_m); // camber line y_c
            double theta = std::atan(nacawing::dyc_dx(naca.m, naca.p, x_m));
            // Calculate the surface coordinates
            double x_surf = x_m - s * y_t * std::sin(theta); // surface x coordinate
            double y_surf = y_c + s * y_t * std::cos(theta); // surface = camber line + half-thickne
            std::complex<double> z = std::complex<double>(x_surf, y_surf);
            return z;
        };

        auto angleAt = [&](double x, bool isUpper) {
            return std::abs(std::arg(toCircle(surfacePoint(x, isUpper)) - zeta_c));
        };


            
        // bissection method to find x such that angleAt(x, isUpper) = phi
        std::vector<std::complex<double>> wall(nc + 1);
        wall[0] = wall[nc] = z_1;

        for (int i = 1; i < nc; ++i) {          
            bool isUpper = (i > nc / 2);
            double phi = 2.0 * std::numbers::pi * std::min(i, nc - i) / nc;
            // bisection on x in [0, x_TE] with angleAt(x, isUpper)
            double phi_real = 0.0;
            double x1 = 0.0;
            double x2 = x_TE;
            double xm = 0.0;
            while (std::abs(phi - phi_real) > 1e-15) {

                xm = (x2+x1) * 0.5;
                phi_real = angleAt(xm, isUpper);
                if (phi_real > phi) {
                    x1 = xm;
                } else {
                    x2 = xm;
                }

            } 
            wall[i] = surfacePoint(xm, isUpper);

            }

    // 6. for each j, for each i: R_j, r_ij (eq. 14),
    //    zeta = zc + r_ij * std::polar(1.0, theta_i), z = fromCircle(zeta)
        std::vector<double> R;
        for (int j = 0; j <= nc; ++j) {
            // eq. 13 - perfect circle of radius R; j starts at 0 so (j +1 -1)*2*pi
            R.push_back(R1 * std::exp((static_cast<double>(j)) * 2.0 * std::numbers::pi / nc ));
        }

    // eq. 14 - radial coordinate
        mesh.ni = mesh.nj = nc + 1;
        mesh.nodes.resize((nc + 1) * (nc + 1));   // sinon at() écrit hors du vecteur

        for (int i = 0; i <= nc; ++i) {
            double theta_i = -2.0 * std::numbers::pi * i / nc;
            double r_i0 = std::abs(toCircle(wall[i]) - zeta_c);   // rayon du quasi-cercle sur ce rayon
            
            for (int j = 0; j <= nc; ++j) {
                double r_ij = (r_i0 * (R[nc] - R[0]) + R[nc] * (R[j] - R[0])) / (R[nc] - R[0]);  // radial coordinate ; eq.14 - r_ij avec l'éq. 14, en utilisant R[0], R[j], R[nc]
                std::complex<double> zeta = zeta_c + r_ij * std::polar(1.0, theta_i);
                std::complex<double> z = fromCircle(zeta);
                mesh.at(i, j) = {z.real(), z.imag()};
            }
        }

        return mesh;    
        }
} // namespace mesh