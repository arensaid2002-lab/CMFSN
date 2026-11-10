#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <numbers>
#include <fstream>
#include <sstream>
#include <filesystem>
#include "BeamMesh.h"

// Tests for the beam mesh (Mesh/BeamMesh.h).
// Axes as in 3Dwing.hpp: x along the chord, y along the span (root y = 0, tip y = b/2), z up.
// Idea of BeamMesh: one camber point at x = f (fraction of the chord) is extruded with
// wingSurface, so each node gets the local chord, sweep, twist and dihedral for free.
// Expected values are computed by hand from the formulas (NACA Report 460 for y_c,
// 3Dwing.hpp for the planform).

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
    const double deg = std::numbers::pi / 180.0;

    // NACA 2412 camber (m = 0.02, p = 0.4), beam at f = 0.38 of the chord.
    const double m = 0.02, p = 0.4, f = 0.38;
    // f < p, so y_c = m/p^2 * (2 p f - f^2) = 0.125 * (0.304 - 0.1444) = 0.01995
    const double yc = m / (p * p) * (2.0 * p * f - f * f);

    // ---------------------------------------------------------------------
    // 1. nacaCambrurePoint: ONE point (f, y_c(f)) on the unit-chord camber line.
    // ---------------------------------------------------------------------
    std::vector<nacawing::Point2D> pt = mesh::nacaCambrurePoint(m, p, f);
    failures += check("nacaCambrurePoint number of points", static_cast<double>(pt.size()), 1.0);
    failures += check("nacaCambrurePoint x", pt[0].x, f);
    failures += check("nacaCambrurePoint y = y_c(0.38)", pt[0].y, yc);
    failures += check("nacaCambrurePoint 0012 y = 0", mesh::nacaCambrurePoint(0.0, 0.0, f)[0].y, 0.0);

    // ---------------------------------------------------------------------
    // 2. BeamMesh, no twist, no dihedral, default SineTip spacing.
    //    c_r = 1, c_t = 0.5, b = 8 (so b/2 = 4), sweep 20 deg, nj = 6 elements -> 7 nodes.
    // ---------------------------------------------------------------------
    nacawing::WingParams w{};
    w.b = 8.0;  w.c_r = 1.0;  w.c_t = 0.5;
    w.sweepLE = 20.0;  w.twistTip = 0.0;  w.dihedral = 0.0;
    const int nj = 6;

    std::vector<mesh::NastranData> g = mesh::BeamMesh(w, pt, nj);
    failures += check("BeamMesh number of nodes = nj + 1", static_cast<double>(g.size()), nj + 1.0);
    failures += check("BeamMesh first id = 1", g.front().id, 1.0);
    failures += check("BeamMesh last id = nj + 1", g.back().id, nj + 1.0);

    // Root node (y = 0): x = f * c_r, z = y_c * c_r.
    failures += check("root x = f c_r", g[0].x, f * 1.0);
    failures += check("root y = 0", g[0].y, 0.0);
    failures += check("root z = y_c c_r", g[0].z, yc * 1.0);

    // Tip node (y = b/2 = 4): c = c_t = 0.5, x = x_LE + f c = 4 tan(20 deg) + 0.19.
    failures += check("tip x = 4 tan20 + f c_t", g[nj].x, 4.0 * std::tan(20.0 * deg) + f * 0.5);
    failures += check("tip y = b/2", g[nj].y, 4.0);
    failures += check("tip z = y_c c_t", g[nj].z, yc * 0.5);

    // Middle node j = 3, SineTip: y = 4 sin(pi/2 * 3/6) = 4 sin(pi/4) = 2.828427...
    // chord c = c_r + (c_t - c_r) * y/(b/2) (linear taper).
    const double y3 = 4.0 * std::sin(std::numbers::pi / 4.0);
    const double c3 = 1.0 + (0.5 - 1.0) * y3 / 4.0;
    failures += check("node 3 y (SineTip)", g[3].y, y3);
    failures += check("node 3 x = y tan20 + f c", g[3].x, y3 * std::tan(20.0 * deg) + f * c3);

    // Every node must be farther out than the previous one (no duplicates, no reversal).
    bool increasing = true;
    for (std::size_t j = 1; j < g.size(); j++)
        if (g[j].y <= g[j - 1].y) increasing = false;
    failures += check("y strictly increasing root -> tip", increasing ? 1.0 : 0.0, 1.0);

    // ---------------------------------------------------------------------
    // 3. Uniform spacing: y_j = (b/2) j / nj.
    // ---------------------------------------------------------------------
    std::vector<mesh::NastranData> gu = mesh::BeamMesh(w, pt, nj, nacawing::SpanSpacing::Uniform);
    failures += check("Uniform node 2 y = 4 * 2/6", gu[2].y, 4.0 * 2.0 / 6.0);

    // ---------------------------------------------------------------------
    // 4. Twist and dihedral at the tip: tip twist -4 deg about c/4, dihedral 5 deg.
    //    Local point before twist: x = f c, z = y_c c, pivot x_p = 0.25 c (c = c_t = 0.5).
    //    xr = x_p + (x - x_p) cos(eps) + z sin(eps)
    //    zr = -(x - x_p) sin(eps) + z cos(eps)
    //    then x = x_LE + xr, z = zr + y tan(dihedral).
    // ---------------------------------------------------------------------
    nacawing::WingParams wt = w;
    wt.twistTip = -4.0;  wt.dihedral = 5.0;
    std::vector<mesh::NastranData> gt = mesh::BeamMesh(wt, pt, nj);
    {
        const double c = 0.5, eps = -4.0 * deg, xp = 0.25 * c;
        const double x = f * c, z = yc * c;
        const double xr = xp + (x - xp) * std::cos(eps) + z * std::sin(eps);
        const double zr = -(x - xp) * std::sin(eps) + z * std::cos(eps);
        failures += check("twisted tip x", gt[nj].x, 4.0 * std::tan(20.0 * deg) + xr);
        failures += check("twisted tip z", gt[nj].z, zr + 4.0 * std::tan(5.0 * deg));
    }
    // The root has no twist and no dihedral lift (eta = 0, y = 0): same as before.
    failures += check("twisted root x unchanged", gt[0].x, g[0].x);
    failures += check("twisted root z unchanged", gt[0].z, g[0].z);

    // ---------------------------------------------------------------------
    // 5. writeNastran: structure-team format, ONE line "id0,x0,y0,z0,id1,x1,..."
    //    -> 4 values per node, no leading or trailing comma.
    // ---------------------------------------------------------------------
    std::string file = (std::filesystem::temp_directory_path() / "test_beam_mesh.txt").string();
    mesh::writeNastran(g, file);

    std::ifstream in(file);
    std::string line;
    std::getline(in, line);
    std::string extra;
    bool onlyOneLine = !std::getline(in, extra) || extra.empty();
    failures += check("writeNastran one line", onlyOneLine ? 1.0 : 0.0, 1.0);
    failures += check("writeNastran no leading comma", (!line.empty() && line.front() != ',') ? 1.0 : 0.0, 1.0);
    failures += check("writeNastran no trailing comma", (!line.empty() && line.back() != ',') ? 1.0 : 0.0, 1.0);

    std::vector<double> values;
    std::stringstream ss(line);
    std::string token;
    while (std::getline(ss, token, ',')) values.push_back(std::stod(token));
    failures += check("writeNastran 4 values per node", static_cast<double>(values.size()), 4.0 * (nj + 1));
    if (values.size() == static_cast<std::size_t>(4 * (nj + 1))) {
        failures += check("writeNastran first id", values[0], 1.0);
        failures += check("writeNastran first x", values[1], f, 1e-5);    // 6 significant digits by default
        failures += check("writeNastran last id", values[4 * nj], nj + 1.0);
        failures += check("writeNastran last y", values[4 * nj + 2], 4.0, 1e-5);
    }
    std::filesystem::remove(file);

    std::cout << (failures == 0 ? "ALL TESTS PASSED\n" : "SOME TESTS FAILED\n");
    return failures == 0 ? 0 : 1;
}
