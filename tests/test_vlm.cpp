#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <numbers>
#include <fstream>
#include <filesystem>
#include "VLMmesh.h"
#include "nacawing/NACA.hpp"

// Tests for the VLM grid (Mesh/VLMmesh.h).
// Axes as in 3Dwing.hpp: x along the chord, y along the span (root y = 0, tip y = b/2), z up.
// Grid order (option B, chosen by the team): k = j*(ni+1) + i, the chord index i varies fastest.
// Expected values are computed by hand from the formulas (NACA Report 460 for y_c).

// Compares two doubles with a tolerance (never use == on doubles).
// Prints the result and returns 1 if the check failed, 0 if it passed.
int check(const std::string& name, double got, double expected, double tol = 1e-9)
{
    if (std::abs(got - expected) > tol) {
        std::cerr << name << " FAILED: got " << got << ", expected " << expected << "\n";
        return 1;
    }
    std::cout << name << " ok (" << got << ")\n";
    return 0;
}

int main()
{
    int failures = 0;

    // ---------------------------------------------------------------------
    // 1. nacaCambrure: camber line y_c(x) at n points, cosine spacing from the LE to the TE.
    // ---------------------------------------------------------------------
    std::vector<nacawing::Point2D> c2412 = mesh::nacaCambrure(0.02, 0.4, 5);   // NACA 2412: m = 0.02, p = 0.4
    failures += check("nacaCambrure number of points", static_cast<double>(c2412.size()), 5.0);
    failures += check("nacaCambrure LE x", c2412.front().x, 0.0);
    failures += check("nacaCambrure LE y", c2412.front().y, 0.0);
    failures += check("nacaCambrure TE x", c2412.back().x, 1.0);
    failures += check("nacaCambrure TE y", c2412.back().y, 0.0);
    // Middle point: cosine spacing puts it at x = 0.5. Since x > p:
    // y_c = m/(1-p)^2 * ((1-2p) + 2p x - x^2) = 0.02/0.36 * (0.2 + 0.4 - 0.25) = 0.019444...
    failures += check("nacaCambrure 2412 y_c(0.5)", c2412[2].y, 0.02 / 0.36 * 0.35);

    // NACA 0012 has no camber (m = p = 0): every y must be 0, with no division by zero (NaN).
    std::vector<nacawing::Point2D> c0012 = mesh::nacaCambrure(0.0, 0.0, 5);
    double maxY = 0.0;
    for (const auto& q : c0012) maxY = std::max(maxY, std::abs(q.y));
    failures += check("nacaCambrure 0012 all y = 0", maxY, 0.0);

    // ---------------------------------------------------------------------
    // 2. cambrureDuProfil: camber taken back from a full airfoil loop.
    //    For a NACA loop the midpoint of (upper i, lower i) is EXACTLY the camber point:
    //    (x_u + x_l)/2 = x_c and (y_u + y_l)/2 = y_c, because the +/- y_t terms cancel.
    // ---------------------------------------------------------------------
    const int n = 9;                                                        // points per surface
    std::vector<nacawing::Point2D> loop = nacawing::sectionLoop(0.02, 0.4, 0.12, n);   // 2n - 1 = 17 points
    std::vector<nacawing::Point2D> cam = mesh::cambrureDuProfil(loop);
    failures += check("cambrureDuProfil number of points", static_cast<double>(cam.size()), static_cast<double>(n));
    failures += check("cambrureDuProfil starts at the LE (x)", cam.front().x, 0.0);
    failures += check("cambrureDuProfil ends at the TE (x)", cam.back().x, 1.0);
    double maxErr = 0.0;                                                    // largest |y - y_c(x)| over all points
    for (const auto& q : cam) maxErr = std::max(maxErr, std::abs(q.y - nacawing::y_c(0.02, 0.4, q.x)));
    failures += check("cambrureDuProfil matches y_c at every point", maxErr, 0.0, 1e-12);

    // ---------------------------------------------------------------------
    // 3. vlmGrid: rectangular wing b = 8 m (b/2 = 4), c = 1 m, no sweep/twist/dihedral,
    //    ni = 4 chord panels (5 camber points), nj = 6 span panels (7 stations), SineTip.
    // ---------------------------------------------------------------------
    nacawing::WingParams rect{.b = 8.0, .c_r = 1.0, .c_t = 1.0,
                              .sweepLE = 0.0, .twistTip = 0.0, .dihedral = 0.0};
    mesh::VlmGrid g = mesh::vlmGrid(rect, c2412, 6);
    failures += check("vlmGrid ni", g.ni, 4.0);
    failures += check("vlmGrid nj", g.nj, 6.0);
    failures += check("vlmGrid P.size = (ni+1)(nj+1)", static_cast<double>(g.P.size()), 35.0);

    // Option B order: k = j*(ni+1) + i.
    failures += check("vlmGrid P[1] is i=1 at the root (y)", g.P[1].y, 0.0);
    failures += check("vlmGrid P[1] is i=1 at the root (x)", g.P[1].x, c2412[1].x);   // c = 1, so x = s
    failures += check("vlmGrid P[ni+1] is the LE of station 1 (x)", g.P[5].x, 0.0);
    failures += check("vlmGrid last point = tip TE (x)", g.P.back().x, 1.0);
    failures += check("vlmGrid last point = tip TE (y)", g.P.back().y, 4.0);
    failures += check("vlmGrid last point = tip TE (z)", g.P.back().z, 0.0);
    // Camber is kept: point (i = 2, j = 0) has z = y_c(0.5) * c.
    failures += check("vlmGrid keeps the camber", g.P[2].z, 0.02 / 0.36 * 0.35);

    // SineTip: y_j = (b/2) sin(pi j / (2 nj)). Station 1: 4 sin(pi/12) = 1.03528 m.
    failures += check("vlmGrid SineTip station 1", g.P[5].y, 4.0 * std::sin(std::numbers::pi / 12.0));
    double firstGap = g.P[1 * 5].y - g.P[0].y;                              // root gap, station 1 - station 0
    double lastGap  = g.P[6 * 5].y - g.P[5 * 5].y;                          // tip gap, station 6 - station 5
    failures += check("vlmGrid SineTip closer at the tip", lastGap < firstGap ? 1.0 : 0.0, 1.0);

    // ---------------------------------------------------------------------
    // 4. writeTecplot: write the file, then read it back.
    //    Tecplot counts POINTS: I = ni+1 (the direction that varies fastest), J = nj+1.
    // ---------------------------------------------------------------------
    std::string file = (std::filesystem::temp_directory_path() / "test_vlm_grid.dat").string();
    mesh::writeTecplot(g, file);
    std::ifstream in(file);
    std::vector<std::string> lines;
    for (std::string line; std::getline(in, line); ) lines.push_back(line);
    failures += check("writeTecplot number of lines (2 header + 35 points)", static_cast<double>(lines.size()), 37.0);
    bool zoneOk = lines.size() > 1 && lines[1] == "ZONE I=5, J=7, DATAPACKING=POINT";
    if (!zoneOk && lines.size() > 1) std::cerr << "  line 2 is: " << lines[1] << "\n";
    failures += check("writeTecplot ZONE line", zoneOk ? 1.0 : 0.0, 1.0);

    return failures == 0 ? 0 : 1;
}
