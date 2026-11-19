#include "incad/InCAD.h"
#include "BeamMesh.h"
#include "nacawing/3Dwing.hpp"
#include "nacawing/NACA.hpp"
#include "nacawing/NACASurface.hpp"

#include <gmsh.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

// Tests du module InCAD. Meme style que les tests d'Eddie : check() renvoie 1 si ca echoue.
// Ces tests utilisent la meme aile que Main/main.cpp (NACA 0012, b = 8, c_r = 1, c_t = 0.5,
// fleche 20 deg, vrillage -4 deg autour de c/4, 12 panneaux SineTip).

int check(const std::string& name, bool ok, double value)
{
    if (!ok) {
        std::cerr << name << " FAILED (" << value << ")\n";
        return 1;
    }
    std::cout << name << " ok (" << value << ")\n";
    return 0;
}

int main()
{
    int failures = 0;

    const nacawing::WingParams wing{
        .b = 8.0, .c_r = 1.0, .c_t = 0.5, .sweepLE = 20.0, .twistTip = -4.0, .dihedral = 0.0,
        .twistPivot = 0.25, .twistExponent = 1.0};
    const int nj = 12;
    const nacawing::NacaParams naca = nacawing::fromDigits(0, 0, 12);
    const std::vector<nacawing::Point2D> loop = nacawing::sectionLoop(naca.m, naca.p, naca.t, 101);

    // 1. Import : les stations de InCAD sont celles de la poutre d'Eddie (meme formule).
    const std::vector<incad::Section> sections =
        incad::sectionsFromPhase2(wing, loop, nj + 1, nacawing::SpanSpacing::SineTip);
    failures += check("nombre de sections = nj + 1", sections.size() == static_cast<std::size_t>(nj + 1),
                      static_cast<double>(sections.size()));
    failures += check("premiere section a y = 0", std::abs(sections.front().y) < 1e-12, sections.front().y);
    failures += check("derniere section a y = b/2", std::abs(sections.back().y - 4.0) < 1e-12, sections.back().y);

    const std::vector<nacawing::Point2D> beamPoint = mesh::nacaCambrurePoint(naca.m, naca.p, 0.38);
    const std::vector<mesh::NastranData> beam = mesh::BeamMesh(wing, beamPoint, nj, nacawing::SpanSpacing::SineTip);
    double worstY = 0.0;
    for (int j = 0; j <= nj; j++)
        worstY = std::max(worstY, std::abs(beam[static_cast<std::size_t>(j)].y - sections[static_cast<std::size_t>(j)].y));
    failures += check("stations InCAD = stations de la poutre", worstY < 1e-12, worstY);

    // 2. Noyau CAD
    incad::initGmsh(false);
    const incad::WingCad cad = incad::buildWing(sections);
    failures += check("un solide cree", cad.volume > 0, cad.volume);
    failures += check("face d'emplanture trouvee", cad.rootCap > 0, cad.rootCap);
    failures += check("face de saumon trouvee", cad.tipCap > 0, cad.tipCap);
    failures += check("surfaces de peau (le loft lisse + la face du bord de fuite)", cad.skin.size() == 2,
                      static_cast<double>(cad.skin.size()));

    // 3. La surface passe par les points importes (tolerance 0.1 mm)
    const double dev = incad::maxDeviation(sections, cad);
    failures += check("ecart max points importes / surface < 1e-4 m", dev < 1e-4, dev);

    // 4. Maillage de surface, puis volume et aire de l'emplanture compares aux FORMULES
    //    (volume attendu : voir expectedVolume ; aire d'emplanture = a0 = aire du profil).
    incad::defineGroups(cad);
    incad::MeshOptions opt;
    incad::meshWing(opt);
    const incad::MeshStats st = incad::meshStats(cad);
    failures += check("maillage : triangles > 0", st.triangles > 0, static_cast<double>(st.triangles));

    const double vRef = incad::expectedVolume(wing, loop);
    const double rel = std::abs(st.volume - vRef) / vRef;
    failures += check("volume (maillage) = formule (ecart relatif < 5e-3)", rel < 5e-3, rel);

    double twiceA0 = 0.0;
    for (std::size_t i = 0; i < loop.size(); i++) {
        const nacawing::Point2D& p = loop[i];
        const nacawing::Point2D& q = loop[(i + 1) % loop.size()];
        twiceA0 += p.x * q.y - q.x * p.y;
    }
    const double a0 = std::abs(twiceA0) * 0.5;
    const double relRoot = std::abs(st.rootArea - a0 * wing.c_r * wing.c_r) / (a0 * wing.c_r * wing.c_r);
    failures += check("aire de l'emplanture = a0 * c_r^2 (ecart relatif < 5e-3)", relRoot < 5e-3, relRoot);

    // 5. Les noeuds de la poutre (hors saumon, qui est sur la face) sont DANS l'aile CAD
    int inside = 0;
    for (int j = 0; j < nj; j++) {
        const mesh::NastranData& n = beam[static_cast<std::size_t>(j)];
        if (gmsh::model::isInside(3, cad.volume, {n.x, n.y, n.z}) == 1)
            inside++;
    }
    failures += check("noeuds de la poutre dans l'aile CAD", inside == nj, inside);

    incad::finalizeGmsh();

    if (failures == 0)
        std::cout << "ALL OK\n";
    return failures;
}
