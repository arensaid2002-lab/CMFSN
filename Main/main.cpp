// cmfsn: main program. Inputs -> geometry -> meshes -> files in Results/.
//
// Axes (same as 3Dwing.hpp, right-handed): x along the chord (LE -> TE), y along the span
// (root y = 0, tip y = b/2), z up. Lengths in metres, angles in degrees.
//
// Output tree (CMFSN_RESULTS_DIR is set by Main/CMakeLists.txt to <repo>/Results):
//   Results/Geometry        airfoil coordinates (.csv)
//   Results/Mesh/Euler      2D O-grid around the airfoil (PLOT3D .x)
//   Results/Mesh/VLM        vortex-lattice panels on the camber surface (Tecplot .dat)
//   Results/Mesh/Structure  beam nodes on the VLM span stations (one line id,x,y,z,..., structure-team format)
//   Results/Solution        solver results (later)

#include "nacawing/NACA.hpp"
#include "nacawing/NACASurface.hpp"
#include "nacawing/3Dwing.hpp"
#include "OGrid.h"
#include "Plot3D.h"
#include "VLMmesh.h"
#include "BeamMesh.h"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace fs = std::filesystem;

int main()
{
    // =====================================================================
    // 1. INPUTS
    // =====================================================================

    // 1.1 Airfoil: NACA M P TT -> m = M/100 (max camber), p = P/10 (its position), t = TT/100 (thickness),
    //     all as fractions of the chord.

    const int M = 0, P = 0, TT = 12;                       // NACA 
    const int nSurface = 1001;                              // points per surface for the geometry .csv

    
    const nacawing::NacaParams airfoil = nacawing::fromDigits(M, P, TT);
    const std::string airfoilName = "NACA" + std::to_string(M) + std::to_string(P) + (TT < 10 ? "0" : "") + std::to_string(TT);

    // 1.2 Wing (see 3Dwing.hpp). Every field is named so the order cannot be mixed up.
    const nacawing::WingParams wing{
        .b             = 8.0,    // [m]   full span, tip to tip (we build one half, y = 0 .. b/2)
        .c_r           = 1.0,    // [m]   root chord
        .c_t           = 0.5,    // [m]   tip chord
        .sweepLE       = 20.0,   // [deg] leading-edge sweep
        .twistTip      = -4.0,   // [deg] tip twist (negative = washout, nose down)
        .dihedral      = 0.0,    // [deg] dihedral
        .twistPivot    = 0.25,   // [-]   twist turns around x = 0.25 c (quarter chord)
        .twistExponent = 1.0,    // [-]   eps(eta) = eps_t * eta^p, 1 = linear
    };

    // 1.3 Mesh sizes
    //  VLM: ni panels along the chord, nj panels along the half span.
    const int vlmNi = 8;
    const int vlmNj = 12;
    const nacawing::SpanSpacing vlmSpacing = nacawing::SpanSpacing::SineTip;   // panels shrink toward the tip

    //  Euler: nc cells in each direction (ni = nj = nc + 1 nodes). nc must be even. PRENDRE PUISSANCE DE 2
    //  generateOGrid only supports symmetric profiles for now, so the Euler mesh uses NACA 0012.
    const nacawing::NacaParams eulerAirfoil = nacawing::fromDigits(0, 0, 12);
    const int eulerNc = 4096;

    //  Structure: beam nodes on the SAME span stations as the VLM (vlmNj + 1 nodes, vlmNj elements,
    //  vlmSpacing), so every node sits on a VLM chord line (load transfer, Jasa et al. 2018, sec. 3.2).
    //  beamXc = chordwise position of the node, fraction of the local chord on the camber line
    //  (0 = LE, 1 = TE). 0.38 is a proposal (OAS wingbox formula); the structure team decides.
    const double beamXc = 0.38;

    // 1.4 Output folders (created if they are missing; nothing happens if they exist).
    const fs::path results    = CMFSN_RESULTS_DIR;
    const fs::path geomDir    = results / "Geometry";
    const fs::path eulerDir   = results / "Mesh" / "Euler";
    const fs::path vlmDir     = results / "Mesh" / "VLM";
    const fs::path structDir  = results / "Mesh" / "Structure";
    const fs::path solDir     = results / "Solution";
    for (const fs::path& dir : {geomDir, eulerDir, vlmDir, structDir, solDir})
        fs::create_directories(dir);

    std::vector<std::string> summary;   // one line per file written, printed at the end

    // =====================================================================
    // 2. GEOMETRY
    // =====================================================================
    // Upper and lower surfaces at unit chord, same points as CAD/src/main.cpp.
    {
        std::vector<nacawing::Point2D> upper = nacawing::upperSurface(airfoil.m, airfoil.p, airfoil.t, nSurface);
        std::vector<nacawing::Point2D> lower = nacawing::lowerSurface(airfoil.m, airfoil.p, airfoil.t, nSurface);
        const fs::path file = geomDir / ("airfoil_" + airfoilName + ".csv");
        std::ofstream out(file);
        out << "x_upper,y_upper,x_lower,y_lower\n";
        for (std::size_t i = 0; i < upper.size(); ++i)
            out << upper[i].x << "," << upper[i].y << "," << lower[i].x << "," << lower[i].y << "\n";
        summary.push_back(file.generic_string() + "  (" + std::to_string(upper.size()) + " points per surface)");
    }

    // =====================================================================
    // 3. MESH
    // =====================================================================

    // 3.1 Euler: 2D O-grid around the airfoil (Vassberg-Jameson conformal mapping).
    {
        mesh::StructuredMesh2D grid = mesh::generateOGrid(eulerAirfoil, eulerNc);
        std::cout << "Euler O-grid: ni = " << grid.ni << ", nj = " << grid.nj << " nodes\n";
        const fs::path file = eulerDir / ("ogrid_NACA0012_nc" + std::to_string(eulerNc) + ".x");
        // PLOT3D 2D, one block: "ni nj", then all ni*nj x, then all ni*nj y (readPlot3D reads it back).
        mesh::writePlot3D(grid, file.string());
        summary.push_back(file.generic_string() + "  (" + std::to_string(grid.ni) + " x " + std::to_string(grid.nj)
                          + " nodes, " + std::to_string(grid.nodes.size()) + " points)");
    }

    // 3.2 VLM: panels on the camber surface of the wing.
    //     nacaCambrure gives ni + 1 points on the camber line (cosine spacing in x),
    //     vlmGrid places that line at nj + 1 span stations with chord, sweep, twist and dihedral.
    {
        std::vector<nacawing::Point2D> camber = mesh::nacaCambrure(airfoil.m, airfoil.p, vlmNi + 1);
        mesh::VlmGrid grid = mesh::vlmGrid(wing, camber, vlmNj, vlmSpacing);
        const fs::path file = vlmDir / ("vlm_" + airfoilName + "_ni" + std::to_string(vlmNi)
                                        + "_nj" + std::to_string(vlmNj) + ".dat");
        mesh::writeTecplot(grid, file.string());
        summary.push_back(file.generic_string() + "  (" + std::to_string(grid.ni) + " x " + std::to_string(grid.nj)
                          + " panels, " + std::to_string(grid.P.size()) + " points)");
    }

    // 3.3 Structure: beam nodes. nacaCambrurePoint gives ONE camber point (beamXc, y_c), and BeamMesh
    //     extrudes it with wingSurface on the VLM stations (same vlmNj and vlmSpacing), so it gets
    //     the local chord, sweep, twist and dihedral. Ids start at 1; element e joins nodes e and e + 1.
    {
        std::vector<nacawing::Point2D> beamPoint = mesh::nacaCambrurePoint(airfoil.m, airfoil.p, beamXc);
        std::vector<mesh::NastranData> nodes = mesh::BeamMesh(wing, beamPoint, vlmNj, vlmSpacing);
        const fs::path file = structDir / ("beam_" + airfoilName + "_nj" + std::to_string(vlmNj) + ".txt");
        mesh::writeNastran(nodes, file.string());
        summary.push_back(file.generic_string() + "  (" + std::to_string(nodes.size()) + " nodes, "
                          + std::to_string(nodes.size() - 1) + " elements, x/c = " + std::to_string(beamXc) + ")");
    }

    // =====================================================================
    // 4. SUMMARY
    // =====================================================================
    std::cout << "\nFiles written:\n";
    for (const std::string& line : summary)
        std::cout << "  " << line << "\n";

    return 0;
}
