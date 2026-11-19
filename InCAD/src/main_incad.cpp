// cmfsn_incad : sections de la Phase 2 -> aile CAD (OpenCASCADE) -> maillage Gmsh.
//
// ATTENTION : les entrees de la section 1 doivent rester IDENTIQUES a celles de Main/main.cpp
// (profil, aile, nombre de stations du VLM), sinon l'aile CAD n'est pas celle du maillage
// out-of-CAD.
//
// Fichiers ecrits (CMFSN_RESULTS_DIR = <depot>/Results) :
//   Results/Geometry/wing_NACA0012.step         la geometrie CAD (ouvrable dans tout logiciel CAD)
//   Results/Mesh/InCAD/wing_NACA0012_incad.msh  le maillage de surface (format Gmsh)
//   Results/Mesh/InCAD/incad_report.txt         les chiffres de verification

#include "incad/InCAD.h"
#include "nacawing/3Dwing.hpp"
#include "nacawing/NACA.hpp"
#include "nacawing/NACASurface.hpp"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

namespace fs = std::filesystem;

int main()
{
    // =====================================================================
    // 1. ENTREES (copie de Main/main.cpp)
    // =====================================================================
    const int M = 0, P = 0, TT = 12;                                    // NACA 0012
    const nacawing::NacaParams airfoil = nacawing::fromDigits(M, P, TT);
    const std::string airfoilName = "NACA" + std::to_string(M) + std::to_string(P)
                                    + (TT < 10 ? "0" : "") + std::to_string(TT);

    const nacawing::WingParams wing{
        .b             = 8.0,    // [m]   envergure totale (on construit la demi-aile y = 0 .. b/2)
        .c_r           = 1.0,    // [m]   corde a l'emplanture
        .c_t           = 0.5,    // [m]   corde au saumon
        .sweepLE       = 20.0,   // [deg] fleche du bord d'attaque
        .twistTip      = -4.0,   // [deg] vrillage au saumon
        .dihedral      = 0.0,    // [deg] diedre
        .twistPivot    = 0.25,   // [-]   le vrillage tourne autour de c/4
        .twistExponent = 1.0,    // [-]   vrillage lineaire
    };
    const int vlmNj = 12;                                               // panneaux en envergure
    const nacawing::SpanSpacing spacing = nacawing::SpanSpacing::SineTip;

    // Specifique a InCAD
    const int nProfile = 101;      // points par surface (extrados / intrados) de chaque section
    incad::MeshOptions mesh;       // tailles de maille : voir InCAD.h
    mesh.sizeMax = 0.10;
    mesh.sizeMin = 0.002;
    mesh.curvaturePoints = 60;
    mesh.dim = 2;

    const fs::path results  = CMFSN_RESULTS_DIR;
    const fs::path geomDir  = results / "Geometry";
    const fs::path meshDir  = results / "Mesh" / "InCAD";
    fs::create_directories(geomDir);
    fs::create_directories(meshDir);

    // =====================================================================
    // 2. IMPORT des sections depuis la Phase 2
    // =====================================================================
    const std::vector<nacawing::Point2D> loop = nacawing::sectionLoop(airfoil.m, airfoil.p, airfoil.t, nProfile);
    const std::vector<incad::Section> sections = incad::sectionsFromPhase2(wing, loop, vlmNj + 1, spacing);

    // =====================================================================
    // 3. NOYAU CAD : B-splines + loft, puis maillage
    // =====================================================================
    incad::initGmsh(false);
    const incad::WingCad cad = incad::buildWing(sections);
    incad::defineGroups(cad);

    const double dev     = incad::maxDeviation(sections, cad);
    const double volRef  = incad::expectedVolume(wing, loop);

    incad::writeFile((geomDir / ("wing_" + airfoilName + ".step")).string());

    incad::meshWing(mesh);
    const incad::MeshStats st = incad::meshStats(cad);
    const fs::path mshFile = meshDir / ("wing_" + airfoilName + "_incad.msh");
    incad::writeFile(mshFile.string());
    incad::finalizeGmsh();

    // =====================================================================
    // 4. RAPPORT
    // =====================================================================
    std::ostringstream rep;
    rep << std::setprecision(8);
    rep << "InCAD report - " << airfoilName << "\n";
    rep << "sections importees (Phase 2)      : " << sections.size() << " x " << loop.size() << " points\n";
    rep << "surfaces de peau                  : " << cad.skin.size() << "\n";
    rep << "ecart max points importes/surface : " << dev << " m\n";
    rep << "volume (a partir du maillage)      : " << st.volume << " m3\n";
    rep << "volume attendu (formule)          : " << volRef << " m3\n";
    rep << "ecart relatif de volume           : " << (st.volume - volRef) / volRef << "\n";
    rep << "aire de peau (maillage)           : " << st.skinArea << " m2\n";
    rep << "aire de la face d'emplanture      : " << st.rootArea << " m2\n";
    rep << "maillage : " << st.nodes << " noeuds, " << st.triangles << " triangles\n";
    std::cout << rep.str();
    std::ofstream(meshDir / "incad_report.txt") << rep.str();

    std::cout << "\nFichiers ecrits :\n  " << (geomDir / ("wing_" + airfoilName + ".step")).generic_string()
              << "\n  " << mshFile.generic_string() << "\n";
    return 0;
}
