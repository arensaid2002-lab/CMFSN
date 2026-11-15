#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <numbers>
#include <algorithm>
#include "nacawing/3Dwing.hpp"

// Tests for the wing geometry (3Dwing.hpp).
// Axes: x along the chord (LE -> TE), y along the span (root y = 0, tip y = b/2), z up.
// Reference wing: c_r = 1 m, c_t = 0.5 m, b = 8 m (so b/2 = 4 m), eta = y / (b/2).
// Planform formulas (Linear chord law):
//   D. Scholz, Aircraft Design ch. 7 "Wing Design", eqs 7.6 to 7.12:
//   https://www.fzt.haw-hamburg.de/pers/Scholz/HOOU/AircraftDesign_7_WingDesign.pdf
// Planform law used by sectionAt:
//   Trapeze: c = c_r [1 - (1 - lambda) eta]
//   Ellipse: c = c_r sqrt(1 - eta^2)                      (c_t is not used)
//   x_LE    = y tan(sweep) + xRef (c_r - c)               (the xRef line is straight, swept by "sweep")
//   twist   eps = twistTip eta^p, about the point x_p = twistPivot c   (twistPivot and xRef are separate)

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

double deg(double a) { return a * std::numbers::pi / 180.0; }   // degrees -> radians

// Local chord measured on the 3D section: x(TE) - x(LE) of the chord line (0,0) -> (1,0).
// Only valid with twistTip = 0 (a twisted chord line is shorter in x).
// (There is no chordAt any more: the chord law is written inside sectionAt.)
double chord(const nacawing::WingParams& w, double y)
{
    std::vector<nacawing::Point2D> line{{0.0, 0.0}, {1.0, 0.0}};
    std::vector<nacawing::Point3D> s = nacawing::sectionAt(w, line, y);
    return s[1].x - s[0].x;
}

// x of the point at fraction f of the local chord, at span station y.
double xAt(const nacawing::WingParams& w, double f, double y)
{
    std::vector<nacawing::Point2D> pt{{f, 0.0}};
    return nacawing::sectionAt(w, pt, y)[0].x;
}

int main()
{
    int failures = 0;

    // ---------------------------------------------------------------------
    // Planform helpers, Trapeze (default chord law). Angles are not used here.
    // ---------------------------------------------------------------------
    nacawing::WingParams w{.b = 8.0, .c_r = 1.0, .c_t = 0.5,
                           .sweep = 0.0, .twistTip = 0.0, .dihedral = 0.0};

    // lambda = c_t / c_r = 0.5 / 1 = 0.5
    failures += check("taperRatio", nacawing::taperRatio(w), 0.5);

    // c(eta) = c_r * (1 - (1 - lambda) * eta)
    failures += check("Trapeze chord root (y = 0)", chord(w, 0.0), 1.0);    // eta = 0   -> c_r
    failures += check("Trapeze chord mid (y = 2)",  chord(w, 2.0), 0.75);   // eta = 0.5 -> 0.75
    failures += check("Trapeze chord tip (y = 4)",  chord(w, 4.0), 0.5);    // eta = 1   -> c_t

    // S = (b/2) * c_r * (1 + lambda) = 4 * 1 * 1.5 = 6 m^2
    failures += check("area", nacawing::area(w), 6.0);
    // AR = b^2 / S = 64 / 6 = 10.666...
    failures += check("aspectRatio", nacawing::aspectRatio(w), 64.0 / 6.0);
    // c_bar = (2/3) * c_r * (1 + lambda + lambda^2) / (1 + lambda) = (2/3) * 1.75 / 1.5 = 0.777...
    failures += check("meanAeroChord", nacawing::meanAeroChord(w), (2.0 / 3.0) * 1.75 / 1.5);

    // Rectangular wing (lambda = 1): the chord must be the same everywhere.
    nacawing::WingParams rect{.b = 6.0, .c_r = 1.0, .c_t = 1.0,
                              .sweep = 0.0, .twistTip = 0.0, .dihedral = 0.0};
    failures += check("rectangular chord tip", chord(rect, 3.0), 1.0);
    failures += check("rectangular meanAeroChord", nacawing::meanAeroChord(rect), 1.0);   // c_bar = c for a rectangle

    // ---------------------------------------------------------------------
    // Ellipse chord law: c = c_r sqrt(1 - eta^2), c_t ignored.
    // ---------------------------------------------------------------------
    nacawing::WingParams ell = w;
    ell.type = nacawing::ChordLaw::Ellipse;
    failures += check("Ellipse chord root (y = 0)", chord(ell, 0.0), 1.0);
    failures += check("Ellipse chord y = b/4", chord(ell, 2.0), std::sqrt(0.75));   // eta = 0.5 -> sqrt(1 - 1/4)
    double cTip = chord(ell, 4.0);                                                 // eta = 1 -> sqrt(0) = 0
    failures += check("Ellipse chord tip is a number (not NaN)", std::isnan(cTip) ? 1.0 : 0.0, 0.0);
    failures += check("Ellipse chord tip = 0 exactly", cTip, 0.0, 0.0);

    // c_t must not change anything for the ellipse: same section with c_t = 0.5 and c_t = 2.
    nacawing::WingParams ell2 = ell;
    ell2.c_t = 2.0;
    failures += check("Ellipse ignores c_t (chord y = 3)", chord(ell2, 3.0), chord(ell, 3.0));
    failures += check("Ellipse ignores c_t (c/4 x, y = 3)", xAt(ell2, 0.25, 3.0), xAt(ell, 0.25, 3.0));

    // TODO (other team): area/aspectRatio/meanAeroChord still use the trapeze formulas.
    // For the ellipse they should give S = pi b c_r / 4, AR = 4 b / (pi c_r), MAC = 8 c_r / (3 pi).

    // ---------------------------------------------------------------------
    // xRef: the line x = xRef c is straight, swept by "sweep". Twist = 0 here
    // (twist would move the points that are not on the pivot).
    // Check: x(xRef c, y) = xRef c_r + y tan(sweep) at several stations.
    // ---------------------------------------------------------------------
    nacawing::WingParams sw{.b = 8.0, .c_r = 1.0, .c_t = 0.5,
                            .sweep = 20.0, .twistTip = 0.0, .dihedral = 0.0};
    const double t20 = std::tan(deg(20.0));
    const double ys[] = {1.0, 2.5, 4.0};

    sw.xRef = 0.0;      // the LE is straight (the old sweepLE behaviour)
    for (double y : ys)
        failures += check("xRef = 0: LE straight, y = " + std::to_string(y), xAt(sw, 0.0, y), y * t20);

    sw.xRef = 0.25;     // the quarter-chord line is straight
    for (double y : ys)
        failures += check("xRef = 0.25: c/4 straight, y = " + std::to_string(y), xAt(sw, 0.25, y), 0.25 * 1.0 + y * t20);

    sw.xRef = 1.0;      // the TE is straight
    for (double y : ys)
        failures += check("xRef = 1: TE straight, y = " + std::to_string(y), xAt(sw, 1.0, y), 1.0 * 1.0 + y * t20);

    // Same for the ellipse: with xRef = 0.25 the c/4 line of an elliptic wing is straight too.
    nacawing::WingParams swEll = sw;
    swEll.type = nacawing::ChordLaw::Ellipse;
    swEll.xRef = 0.25;
    for (double y : ys)
        failures += check("Ellipse xRef = 0.25: c/4 straight, y = " + std::to_string(y), xAt(swEll, 0.25, y), 0.25 + y * t20);

    // ---------------------------------------------------------------------
    // sectionLoop: NACA 2412, n = 5 -> 2n - 1 = 9 points in Selig order
    // (upper TE -> upper surface -> LE -> lower surface -> lower TE).
    // ---------------------------------------------------------------------
    const int n = 5;
    std::vector<nacawing::Point2D> loop = nacawing::sectionLoop(0.02, 0.4, 0.12, n);
    failures += check("sectionLoop size = 2n-1", static_cast<double>(loop.size()), 2.0 * n - 1.0);
    failures += check("sectionLoop[0].x  upper TE", loop[0].x, 1.0001, 1e-4);
    failures += check("sectionLoop[0].y  upper TE above chord", loop[0].y, 0.0013, 1e-4);   // open TE (-0.1015)
    failures += check("sectionLoop[n-1].x  LE", loop[n - 1].x, 0.0);                       // LE kept once, in the middle
    failures += check("sectionLoop[n-1].y  LE", loop[n - 1].y, 0.0);
    failures += check("sectionLoop[last].y  lower TE below chord", loop.back().y, -0.0013, 1e-4);

    // ---------------------------------------------------------------------
    // Regression: the doc's example wing (old code: LE sweep 20 deg, twist -4 deg about c/4).
    // The same wing in the new parameters: xRef = 0.25 and the quarter-chord sweep
    //   tan(L_c/4) = tan(L_LE) - (1/A) (1 - lambda) / (1 + lambda)    (Scholz eq. 7.9, n = 1/4, m = 0)
    // with A = b^2 / S = 64/6 = 10.667 and lambda = 0.5 -> L_c/4 = 18.40 deg.
    // Old worked example (doc section 5), one TE point (1, 0) at y = 2:
    //   eta = 0.5, c = 0.75, eps = -2 deg, x_LE = 2 tan(20 deg) = 0.7279,
    //   x' = 0.1875 + 0.5625 cos(-2 deg) = 0.7497 -> X = 1.4776, z' = -0.5625 sin(-2 deg) = 0.0196
    // ---------------------------------------------------------------------
    const double A = 64.0 / 6.0, lambda = 0.5;
    const double sweepQC = std::atan(t20 - (1.0 / A) * (1.0 - lambda) / (1.0 + lambda)) * 180.0 / std::numbers::pi;
    failures += check("quarter-chord sweep of the doc wing", sweepQC, 18.40, 1e-2);

    nacawing::WingParams doc{.b = 8.0, .c_r = 1.0, .c_t = 0.5,
                             .sweep = sweepQC, .twistTip = -4.0, .dihedral = 0.0, .xRef = 0.25};
    std::vector<nacawing::Point2D> te{{1.0, 0.0}};
    nacawing::Point3D p = nacawing::sectionAt(doc, te, 2.0)[0];
    failures += check("regression old wing TE X (y = 2)", p.x, 1.4776, 1e-4);
    failures += check("regression old wing TE Z (up)", p.z, 0.0196, 1e-4);   // TE goes UP: washout = nose down
    failures += check("regression old wing TE Y (span)", p.y, 2.0);

    // At the root nothing moves (c = c_r = 1, eps = 0, x_LE = 0): output = loop with y = 0 (profile y -> z).
    std::vector<nacawing::Point3D> root = nacawing::sectionAt(doc, loop, 0.0);
    failures += check("sectionAt root keeps size", static_cast<double>(root.size()), static_cast<double>(loop.size()));
    failures += check("sectionAt root x unchanged", root[2].x, loop[2].x);
    failures += check("sectionAt root profile y -> z", root[2].z, loop[2].y);

    // Dihedral lifts the section by y tan(phi): LE at the tip with phi = 10 deg (no twist).
    nacawing::WingParams dih = doc;
    dih.dihedral = 10.0;
    dih.twistTip = 0.0;
    std::vector<nacawing::Point2D> le{{0.0, 0.0}};
    failures += check("sectionAt dihedral lift", nacawing::sectionAt(dih, le, 4.0)[0].z, 4.0 * std::tan(deg(10.0)));

    // ---------------------------------------------------------------------
    // Ellipse tip: c = 0, so every loop point collapses onto one point
    //   (xRef c_r + (b/2) tan(sweep), b/2, (b/2) tan(dihedral)), even with twist.
    // ---------------------------------------------------------------------
    nacawing::WingParams tipW{.b = 8.0, .c_r = 1.0, .c_t = 0.5,
                              .sweep = 10.0, .twistTip = -4.0, .dihedral = 5.0, .xRef = 0.25,
                              .type = nacawing::ChordLaw::Ellipse};
    std::vector<nacawing::Point3D> tip = nacawing::sectionAt(tipW, loop, 4.0);
    double worst = 0.0;     // largest distance from the expected point, over the whole loop
    for (const auto& q : tip) {
        double d = std::abs(q.x - (0.25 + 4.0 * std::tan(deg(10.0))))
                 + std::abs(q.y - 4.0)
                 + std::abs(q.z - 4.0 * std::tan(deg(5.0)));
        if (std::isnan(d)) d = 1.0e9;   // a NaN must fail, not be silently ignored
        worst = std::max(worst, d);
    }
    failures += check("Ellipse tip: all points collapse to one point", worst, 0.0, 1e-12);

    // ---------------------------------------------------------------------
    // Twist exponent: eps = eps_t eta^p. p = 2, eps_t = -8 deg, eta = 0.5 -> eps = -8/4 = -2 deg.
    // Measured from the chord line (0,0) -> (1,0): dx = c cos(eps), dz = -c sin(eps).
    // ---------------------------------------------------------------------
    nacawing::WingParams tw{.b = 8.0, .c_r = 1.0, .c_t = 0.5,
                            .sweep = 0.0, .twistTip = -8.0, .dihedral = 0.0, .twistExponent = 2.0};
    std::vector<nacawing::Point2D> line{{0.0, 0.0}, {1.0, 0.0}};
    std::vector<nacawing::Point3D> s = nacawing::sectionAt(tw, line, 2.0);
    double eps = std::atan2(-(s[1].z - s[0].z), s[1].x - s[0].x) * 180.0 / std::numbers::pi;
    failures += check("twist exponent p = 2 at eta = 0.5: eps = eps_t / 4", eps, -8.0 / 4.0);

    // ---------------------------------------------------------------------
    // Pivot: the point x = twistPivot c is the twist axis, so twist must not move it.
    // twistPivot and xRef are separate fields, so each case uses a DIFFERENT xRef
    // (this checks that the twist really turns about twistPivot, not about xRef).
    // ---------------------------------------------------------------------
    const double pivots[][2] = {{0.25, 0.25},    // {twistPivot, xRef}: both at the default
                                {0.4,  0.25},    // pivot moved, xRef kept
                                {0.25, 0.7}};    // xRef moved, pivot kept
    for (const auto& pr : pivots) {
        nacawing::WingParams flat{.b = 8.0, .c_r = 1.0, .c_t = 0.5,
                                  .sweep = 15.0, .twistTip = 0.0, .dihedral = 0.0,
                                  .twistPivot = pr[0], .xRef = pr[1]};
        nacawing::WingParams twisted = flat;
        twisted.twistTip = -6.0;
        std::vector<nacawing::Point2D> piv{{pr[0], 0.0}};           // the pivot point on the unit chord
        nacawing::Point3D a = nacawing::sectionAt(flat, piv, 3.0)[0];
        nacawing::Point3D b = nacawing::sectionAt(twisted, piv, 3.0)[0];
        std::string tag = "pivot " + std::to_string(pr[0]) + " c (xRef " + std::to_string(pr[1]) + ") unchanged by twist";
        failures += check(tag + " (x)", b.x, a.x);
        failures += check(tag + " (z)", b.z, a.z);
    }

    // ---------------------------------------------------------------------
    // wingSurface: wing[i][j] = point i of the loop at station j, Ns stations from the root (y = 0) to the tip (y = b/2).
    // ---------------------------------------------------------------------
    const int Ns = 5;
    std::vector<std::vector<nacawing::Point3D>> wing = nacawing::wingSurface(doc, loop, Ns);
    failures += check("wingSurface one row per loop point", static_cast<double>(wing.size()), static_cast<double>(loop.size()));
    failures += check("wingSurface Ns stations per row", static_cast<double>(wing[3].size()), Ns);
    failures += check("wingSurface first station at root", wing[0].front().y, 0.0);
    failures += check("wingSurface last station at tip", wing[0].back().y, doc.b / 2.0);
    failures += check("wingSurface stations evenly spaced", wing[0][1].y, 1.0);     // (b/2) / (Ns - 1) = 1
    failures += check("wingSurface column 2 == sectionAt(y = 2)", wing[0][2].x, nacawing::sectionAt(doc, loop, 2.0)[0].x);

    // Same wing, different span: the tip must follow b/2 (catches "y = i" bugs).
    nacawing::WingParams wide = doc;
    wide.b = 16.0;
    failures += check("wingSurface tip follows the span", nacawing::wingSurface(wide, loop, Ns)[0].back().y, 8.0);

    // SineTip (now a field of WingParams): y_j = (b/2) sin(pi j / (2(Ns-1))). b/2 = 4, Ns = 5: y_1 = 4 sin(pi/8).
    nacawing::WingParams sineW = doc;
    sineW.spacing = nacawing::SpanSpacing::SineTip;
    std::vector<std::vector<nacawing::Point3D>> sine = nacawing::wingSurface(sineW, loop, Ns);
    failures += check("SineTip root and tip", sine[0].front().y + sine[0].back().y, 0.0 + doc.b / 2.0);
    failures += check("SineTip station 1", sine[0][1].y, 4.0 * std::sin(std::numbers::pi / 8.0));
    // Last gap (tip) smaller than the first gap (root): 1 = true.
    failures += check("SineTip closer at the tip",
                      (sine[0][4].y - sine[0][3].y) < (sine[0][1].y - sine[0][0].y) ? 1.0 : 0.0, 1.0);

    std::cout << failures << " check(s) failed\n";
    return failures == 0 ? 0 : 1;
}
