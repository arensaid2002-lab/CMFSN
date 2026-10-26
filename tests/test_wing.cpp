#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <numbers>
#include "nacawing/3Dwing.hpp"

// Tests for the trapezoidal wing helpers (3Dwing.hpp).
// Reference wing: c_r = 1 m, c_t = 0.5 m, b = 8 m (so b/2 = 4 m).
// Expected values are computed by hand from the formulas in
// D. Scholz, Aircraft Design ch. 7 "Wing Design", eqs 7.6 to 7.12:
// https://www.fzt.haw-hamburg.de/pers/Scholz/HOOU/AircraftDesign_7_WingDesign.pdf

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

    // Angles are not used by these functions yet, so they are set to 0.
    nacawing::NacaWingParams w{.b = 8.0, .c_r = 1.0, .c_t = 0.5,
                               .sweepLE = 0.0, .twistTip = 0.0, .dihedral = 0.0};

    // lambda = c_t / c_r = 0.5 / 1 = 0.5
    failures += check("taperRatio", nacawing::taperRatio(w), 0.5);

    // c(eta) = c_r * (1 - (1 - lambda) * eta), with eta = z / (b/2)
    failures += check("chordAt root (z = 0)", nacawing::chordAt(w, 0.0), 1.0);   // eta = 0   -> c_r
    failures += check("chordAt mid (z = 2)",  nacawing::chordAt(w, 2.0), 0.75);  // eta = 0.5 -> 0.75
    failures += check("chordAt tip (z = 4)",  nacawing::chordAt(w, 4.0), 0.5);   // eta = 1   -> lambda * c_r

    // S = (b/2) * c_r * (1 + lambda) = 4 * 1 * 1.5 = 6 m^2
    failures += check("area", nacawing::area(w), 6.0);

    // AR = b^2 / S = 64 / 6 = 10.666...
    failures += check("aspectRatio", nacawing::aspectRatio(w), 64.0 / 6.0);

    // c_bar = (2/3) * c_r * (1 + lambda + lambda^2) / (1 + lambda) = (2/3) * 1.75 / 1.5 = 0.777...
    failures += check("meanAeroChord", nacawing::meanAeroChord(w), (2.0 / 3.0) * 1.75 / 1.5);

    // Rectangular wing (lambda = 1): the chord must be the same everywhere.
    nacawing::NacaWingParams rect{.b = 6.0, .c_r = 1.0, .c_t = 1.0,
                                  .sweepLE = 0.0, .twistTip = 0.0, .dihedral = 0.0};
    failures += check("rectangular chordAt tip", nacawing::chordAt(rect, 3.0), 1.0);
    failures += check("rectangular meanAeroChord", nacawing::meanAeroChord(rect), 1.0); // c_bar = c for a rectangle

    // ---------------------------------------------------------------------
    // sectionLoop: NACA 2412, n = 5 -> 2n - 1 = 9 points in Selig order
    // (upper TE -> upper surface -> LE -> lower surface -> lower TE).
    // ---------------------------------------------------------------------
    const int n = 5;
    std::vector<nacawing::Point2D> loop = nacawing::sectionLoop(0.02, 0.4, 0.12, n);
    failures += check("sectionLoop size = 2n-1", static_cast<double>(loop.size()), 2.0 * n - 1.0);
    failures += check("sectionLoop[0].x  upper TE", loop[0].x, 1.0001, 1e-4);
    failures += check("sectionLoop[0].y  upper TE above chord", loop[0].y, 0.0013, 1e-4);   // open TE (-0.1015)
    failures += check("sectionLoop[n-1].x  LE", loop[n - 1].x, 0.0);                    // LE kept once, in the middle
    failures += check("sectionLoop[n-1].y  LE", loop[n - 1].y, 0.0);
    failures += check("sectionLoop[last].y  lower TE below chord", loop.back().y, -0.0013, 1e-4);

    // ---------------------------------------------------------------------
    // sectionAt: the doc's example wing, b = 8, sweep 20 deg, tip twist -4 deg.
    // ---------------------------------------------------------------------
    nacawing::NacaWingParams doc{.b = 8.0, .c_r = 1.0, .c_t = 0.5,
                                 .sweepLE = 20.0, .twistTip = -4.0, .dihedral = 0.0};

    // Worked example (doc section 5), computed by hand from a single TE point (1, 0):
    // z = 2 -> eta = 0.5, c = 0.75, eps = -2 deg, x_LE = 2 tan(20 deg) = 0.7279,
    // x' = 0.1875 + 0.5625 cos(-2 deg) = 0.7497, y' = -0.5625 sin(-2 deg) = 0.0196
    std::vector<nacawing::Point2D> te{{1.0, 0.0}};
    nacawing::Point3D p = nacawing::sectionAt(doc, te, 2.0)[0];
    failures += check("sectionAt worked example X", p.x, 1.4776, 1e-4);
    failures += check("sectionAt worked example Y", p.y, 0.0196, 1e-4);   // TE goes UP: washout = nose down
    failures += check("sectionAt worked example Z", p.z, 2.0);

    // At the root nothing moves (c = c_r = 1, eps = 0, x_LE = 0): output = loop with z = 0.
    std::vector<nacawing::Point3D> root = nacawing::sectionAt(doc, loop, 0.0);
    failures += check("sectionAt root keeps size", static_cast<double>(root.size()), static_cast<double>(loop.size()));
    failures += check("sectionAt root x unchanged", root[2].x, loop[2].x);
    failures += check("sectionAt root y unchanged", root[2].y, loop[2].y);

    // Twist pivots about c/4, so the quarter-chord point (0.25, 0) only moves by the sweep.
    // At the tip: c = 0.5 -> c/4 = 0.125, x_LE = 4 tan(20 deg) = 1.4559
    std::vector<nacawing::Point2D> qc{{0.25, 0.0}};
    nacawing::Point3D q = nacawing::sectionAt(doc, qc, 4.0)[0];
    failures += check("sectionAt c/4 does not move with twist (x)", q.x, 4.0 * std::tan(20.0 * std::numbers::pi / 180.0) + 0.125);
    failures += check("sectionAt c/4 does not move with twist (y)", q.y, 0.0);

    // Dihedral lifts the section by z tan(Gamma): LE at the tip with Gamma = 10 deg.
    // Twist is set to 0 here, otherwise the LE would also move by the twist about c/4.
    nacawing::NacaWingParams dih = doc;
    dih.dihedral = 10.0;
    dih.twistTip = 0.0;
    std::vector<nacawing::Point2D> le{{0.0, 0.0}};
    failures += check("sectionAt dihedral lift", nacawing::sectionAt(dih, le, 4.0)[0].y,
                      4.0 * std::tan(10.0 * std::numbers::pi / 180.0));

    // ---------------------------------------------------------------------
    // wingSurface: Ns sections from the root (z = 0) to the tip (z = b/2).
    // ---------------------------------------------------------------------
    const int Ns = 5;
    std::vector<std::vector<nacawing::Point3D>> wing = nacawing::wingSurface(doc, loop, Ns);
    failures += check("wingSurface number of sections", static_cast<double>(wing.size()), Ns);
    failures += check("wingSurface points per section", static_cast<double>(wing[3].size()), static_cast<double>(loop.size()));
    failures += check("wingSurface first section at root", wing.front()[0].z, 0.0);
    failures += check("wingSurface last section at tip", wing.back()[0].z, doc.b / 2.0);
    failures += check("wingSurface stations evenly spaced", wing[1][0].z, 1.0);     // (b/2) / (Ns - 1) = 1
    failures += check("wingSurface[2] == sectionAt(z = 2)", wing[2][0].x, nacawing::sectionAt(doc, loop, 2.0)[0].x);

    // Same wing, different span: the tip must follow b/2 (catches "z = i" bugs).
    nacawing::NacaWingParams wide = doc;
    wide.b = 16.0;
    failures += check("wingSurface tip follows the span", nacawing::wingSurface(wide, loop, Ns).back()[0].z, 8.0);


    return failures == 0 ? 0 : 1;
}
