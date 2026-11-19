#include "incad/InCAD.h"

#include <gmsh.h>

#include <cmath>
#include <stdexcept>
#include <utility>

namespace incad {

// ---------------------------------------------------------------------------
// 1. IMPORT depuis la Phase 2
// ---------------------------------------------------------------------------
std::vector<Section> sectionsFromPhase2(const nacawing::WingParams& w,
                                        const std::vector<nacawing::Point2D>& loop,
                                        int nStations,
                                        nacawing::SpanSpacing spacing)
{
    if (nStations < 2)
        throw std::invalid_argument("sectionsFromPhase2 : il faut au moins 2 stations");

    std::vector<Section> sections;
    for (int j = 0; j < nStations; j++) {
        // Meme formule de station que nacawing::wingSurface (3Dwing.cpp), pour tomber
        // exactement sur les stations y_j du VLM et de la poutre.
        double eta = static_cast<double>(j) / (nStations - 1);          // 0 = emplanture, 1 = saumon
        if (spacing == nacawing::SpanSpacing::SineTip)
            eta = std::sin(0.5 * 3.14159265358979323846 * eta);          // sin(pi j / (2 (Ns - 1)))
        const double y = 0.5 * w.b * eta;

        Section s;
        s.y = y;
        s.points = nacawing::sectionAt(w, loop, y);   // corde, fleche, vrillage, diedre inclus
        sections.push_back(s);
    }
    return sections;
}

// ---------------------------------------------------------------------------
// 2. Gmsh
// ---------------------------------------------------------------------------
void initGmsh(bool verbose)
{
    gmsh::initialize();
    gmsh::option::setNumber("General.Terminal", verbose ? 1 : 0);
    gmsh::model::add("wing");
}

void finalizeGmsh()
{
    gmsh::finalize();
}

// ---------------------------------------------------------------------------
// 3. Sections -> B-splines -> loft
// ---------------------------------------------------------------------------
WingCad buildWing(const std::vector<Section>& sections)
{
    if (sections.size() < 2)
        throw std::invalid_argument("buildWing : il faut au moins 2 sections");

    std::vector<int> wires;
    for (std::size_t s = 0; s < sections.size(); s++) {
        const std::vector<nacawing::Point3D>& pts = sections[s].points;

        // Un point OpenCASCADE par point de la section.
        std::vector<int> pointTags;
        for (std::size_t k = 0; k < pts.size(); k++)
            pointTags.push_back(gmsh::model::occ::addPoint(pts[k].x, pts[k].y, pts[k].z));

        // Bord de fuite : distance entre le premier et le dernier point (= epaisseur du BF).
        const nacawing::Point3D& a = pts.front();
        const nacawing::Point3D& b = pts.back();
        const double gap = std::sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y)
                                     + (a.z - b.z) * (a.z - b.z));
        if (gap < 1e-9)
            throw std::runtime_error("buildWing : bord de fuite ferme, non gere (utiliser un NACA standard)");

        // Une seule B-spline interpolante passe par TOUS les points : BF extrados -> BA -> BF intrados.
        // Au bord d'attaque elle est lisse (un seul morceau), il n'y a pas d'angle artificiel.
        const int spline = gmsh::model::occ::addSpline(pointTags);

        // Segment droit qui ferme le contour au bord de fuite (de l'intrados vers l'extrados).
        const int closing = gmsh::model::occ::addLine(pointTags.back(), pointTags.front());

        wires.push_back(gmsh::model::occ::addWire({spline, closing}));
    }

    // Loft : une surface lisse qui passe par tous les contours, fermee par deux faces planes
    // (makeSolid = true), ce qui donne un solide.
    gmsh::vectorpair out;
    gmsh::model::occ::addThruSections(wires, out, -1, true, false);
    gmsh::model::occ::synchronize();

    WingCad cad;
    for (std::size_t i = 0; i < out.size(); i++)
        if (out[i].first == 3)
            cad.volume = out[i].second;
    if (cad.volume < 0)
        throw std::runtime_error("buildWing : le loft n'a pas donne de solide");

    // Classement des faces du solide d'apres leur boite englobante :
    //   y_max - y_min ~ 0  et y ~ 0      -> face d'emplanture
    //   y_max - y_min ~ 0  et y ~ b/2    -> face de saumon
    //   sinon                            -> peau
    const double yTip = sections.back().y;
    const double tol = 1e-6;
    gmsh::vectorpair faces;
    gmsh::model::getBoundary({{3, cad.volume}}, faces, false, false, false);
    for (std::size_t i = 0; i < faces.size(); i++) {
        const int tag = faces[i].second;
        double xmin, ymin, zmin, xmax, ymax, zmax;
        gmsh::model::getBoundingBox(2, tag, xmin, ymin, zmin, xmax, ymax, zmax);
        if (ymax - ymin < tol && std::abs(ymax) < tol)
            cad.rootCap = tag;
        else if (ymax - ymin < tol && std::abs(ymin - yTip) < tol)
            cad.tipCap = tag;
        else
            cad.skin.push_back(tag);
    }
    return cad;
}

void defineGroups(const WingCad& cad)
{
    gmsh::model::addPhysicalGroup(2, cad.skin, 1, "Skin");
    if (cad.rootCap > 0)
        gmsh::model::addPhysicalGroup(2, {cad.rootCap}, 2, "Root");
    if (cad.tipCap > 0)
        gmsh::model::addPhysicalGroup(2, {cad.tipCap}, 3, "Tip");
    gmsh::model::addPhysicalGroup(3, {cad.volume}, 4, "Wing");
}

void meshWing(const MeshOptions& options)
{
    gmsh::option::setNumber("Mesh.MeshSizeMax", options.sizeMax);
    gmsh::option::setNumber("Mesh.MeshSizeMin", options.sizeMin);
    // Raffinement automatique la ou la surface se courbe (bord d'attaque, saumon).
    gmsh::option::setNumber("Mesh.MeshSizeFromCurvature", options.curvaturePoints);
    gmsh::model::mesh::generate(options.dim);
}

// ---------------------------------------------------------------------------
// 4. Verifications
// ---------------------------------------------------------------------------
double maxDeviation(const std::vector<Section>& sections, const WingCad& cad)
{
    std::vector<int> faces = cad.skin;
    faces.push_back(cad.rootCap);
    faces.push_back(cad.tipCap);

    double worst = 0.0;
    for (std::size_t s = 0; s < sections.size(); s++) {
        const std::vector<nacawing::Point3D>& pts = sections[s].points;
        std::vector<double> coord;
        for (std::size_t k = 0; k < pts.size(); k++) {
            coord.push_back(pts[k].x);
            coord.push_back(pts[k].y);
            coord.push_back(pts[k].z);
        }

        // Pour chaque point : distance a la face la plus proche (projection orthogonale).
        std::vector<double> best(pts.size(), 1e30);
        for (std::size_t f = 0; f < faces.size(); f++) {
            if (faces[f] < 0)
                continue;
            std::vector<double> closest, param;
            gmsh::model::getClosestPoint(2, faces[f], coord, closest, param);
            for (std::size_t k = 0; k < pts.size(); k++) {
                const double dx = coord[3 * k] - closest[3 * k];
                const double dy = coord[3 * k + 1] - closest[3 * k + 1];
                const double dz = coord[3 * k + 2] - closest[3 * k + 2];
                const double d = std::sqrt(dx * dx + dy * dy + dz * dz);
                if (d < best[k])
                    best[k] = d;
            }
        }
        for (std::size_t k = 0; k < pts.size(); k++)
            if (best[k] > worst)
                worst = best[k];
    }
    return worst;
}

double expectedVolume(const nacawing::WingParams& w, const std::vector<nacawing::Point2D>& loop)
{
    // a0 : aire du polygone (formule du lacet). La fermeture dernier -> premier point est
    // justement le segment de bord de fuite.
    double twiceArea = 0.0;
    const std::size_t n = loop.size();
    for (std::size_t i = 0; i < n; i++) {
        const nacawing::Point2D& p = loop[i];
        const nacawing::Point2D& q = loop[(i + 1) % n];
        twiceArea += p.x * q.y - q.x * p.y;
    }
    const double a0 = std::abs(twiceArea) * 0.5;

    const double lambda = w.c_t / w.c_r;                                   // effilement
    return a0 * w.c_r * w.c_r * (0.5 * w.b) * (1.0 + lambda + lambda * lambda) / 3.0;
}

// Somme, sur les triangles d'une face, de
//   - la contribution au volume : p0 . (p1 x p2) / 6   (theoreme de la divergence)
//   - l'aire : |(p1 - p0) x (p2 - p0)| / 2
// xyz[3 * tag + c] donne la coordonnee c du noeud de numero tag.
static void integrateFace(int faceTag, const std::vector<double>& xyz, double& volume, double& area)
{
    std::vector<int> types;
    std::vector<std::vector<std::size_t>> elemTags, nodeTags;
    gmsh::model::mesh::getElements(types, elemTags, nodeTags, 2, faceTag);

    for (std::size_t t = 0; t < types.size(); t++) {
        if (types[t] != 2)                       // 2 = triangle a 3 noeuds
            continue;
        for (std::size_t e = 0; e < elemTags[t].size(); e++) {
            const double* p0 = &xyz[3 * nodeTags[t][3 * e]];
            const double* p1 = &xyz[3 * nodeTags[t][3 * e + 1]];
            const double* p2 = &xyz[3 * nodeTags[t][3 * e + 2]];

            // p1 x p2
            const double cx = p1[1] * p2[2] - p1[2] * p2[1];
            const double cy = p1[2] * p2[0] - p1[0] * p2[2];
            const double cz = p1[0] * p2[1] - p1[1] * p2[0];
            volume += (p0[0] * cx + p0[1] * cy + p0[2] * cz) / 6.0;

            // (p1 - p0) x (p2 - p0)
            const double ux = p1[0] - p0[0], uy = p1[1] - p0[1], uz = p1[2] - p0[2];
            const double vx = p2[0] - p0[0], vy = p2[1] - p0[1], vz = p2[2] - p0[2];
            const double nx = uy * vz - uz * vy;
            const double ny = uz * vx - ux * vz;
            const double nz = ux * vy - uy * vx;
            area += 0.5 * std::sqrt(nx * nx + ny * ny + nz * nz);
        }
    }
}

MeshStats meshStats(const WingCad& cad)
{
    MeshStats st;

    // Tous les noeuds, ranges par numero pour les retrouver vite.
    std::vector<std::size_t> nodeTags;
    std::vector<double> coord, param;
    gmsh::model::mesh::getNodes(nodeTags, coord, param);
    st.nodes = nodeTags.size();
    std::size_t maxTag = 0;
    for (std::size_t i = 0; i < nodeTags.size(); i++)
        if (nodeTags[i] > maxTag)
            maxTag = nodeTags[i];
    std::vector<double> xyz(3 * (maxTag + 1), 0.0);
    for (std::size_t i = 0; i < nodeTags.size(); i++)
        for (int c = 0; c < 3; c++)
            xyz[3 * nodeTags[i] + static_cast<std::size_t>(c)] = coord[3 * i + static_cast<std::size_t>(c)];

    // Comptage des elements
    std::vector<int> types;
    std::vector<std::vector<std::size_t>> elemTags, nodeTagsPerType;
    gmsh::model::mesh::getElements(types, elemTags, nodeTagsPerType, 2);
    for (std::size_t i = 0; i < types.size(); i++)
        if (types[i] == 2)
            st.triangles += elemTags[i].size();

    types.clear();
    elemTags.clear();
    nodeTagsPerType.clear();
    gmsh::model::mesh::getElements(types, elemTags, nodeTagsPerType, 3);
    for (std::size_t i = 0; i < types.size(); i++)
        if (types[i] == 4)                       // 4 = tetraedre a 4 noeuds
            st.tetrahedra += elemTags[i].size();

    // Volume = somme sur toutes les faces du solide ; aire de peau = faces de peau seulement.
    // (Gmsh oriente les triangles vers l'exterieur du solide ; le resultat est verifie
    //  contre la formule dans tests/test_incad.cpp.)
    double volume = 0.0, area = 0.0;
    for (std::size_t i = 0; i < cad.skin.size(); i++) {
        double v = 0.0, a = 0.0;
        integrateFace(cad.skin[i], xyz, v, a);
        volume += v;
        area += a;
    }
    st.skinArea = area;
    if (cad.rootCap > 0) {
        double v = 0.0, a = 0.0;
        integrateFace(cad.rootCap, xyz, v, a);
        volume += v;
        st.rootArea = a;
    }
    if (cad.tipCap > 0) {
        double v = 0.0, a = 0.0;
        integrateFace(cad.tipCap, xyz, v, a);
        volume += v;
    }
    st.volume = std::abs(volume);
    return st;
}

void writeFile(const std::string& filename)
{
    gmsh::write(filename);
}

}  // namespace incad
