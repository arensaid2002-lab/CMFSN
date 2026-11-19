#pragma once
// InCAD : noyau CAD (OpenCASCADE, via l'API C++ de Gmsh).
//
// Idee en 4 etapes :
//   1. Phase 2 (out-of-CAD) : nacawing::sectionLoop + sectionAt donnent, a chaque station y,
//      les points du profil NACA deja a la bonne corde, avec fleche, vrillage et diedre.
//   2. IMPORT : chaque section devient une B-spline interpolante dans OpenCASCADE.
//   3. LOFT : les B-splines sont reliees en un solide (l'aile, demi-envergure y = 0 .. b/2).
//   4. Gmsh maille ce solide : c'est Gmsh qui place les noeuds, pas nous.
//
// Axes (comme 3Dwing.hpp) : x le long de la corde (BA -> BF), y le long de l'envergure
// (emplanture y = 0, saumon y = b/2), z vers le haut. Longueurs en metres.

#include "nacawing/3Dwing.hpp"
#include "nacawing/NACASurface.hpp"

#include <string>
#include <vector>

namespace incad {

// Une section = le contour du profil place a une station y de l'envergure.
// Ordre des points (Selig, comme sectionLoop) : BF extrados -> BA -> BF intrados.
struct Section {
    double y = 0.0;                          // station d'envergure [m]
    std::vector<nacawing::Point3D> points;   // points du contour, deja en 3D
};

// IMPORT depuis la Phase 2 : calcule les sections aux nStations stations de l'envergure.
// Avec les memes nStations et spacing que le VLM (vlmNj + 1, SineTip), l'aile CAD passe
// exactement par les memes stations y_j que le VLM et la poutre.
std::vector<Section> sectionsFromPhase2(const nacawing::WingParams& w,
                                        const std::vector<nacawing::Point2D>& loop,
                                        int nStations,
                                        nacawing::SpanSpacing spacing);

// Ce que buildWing() cree dans Gmsh (numeros d'entites OpenCASCADE).
struct WingCad {
    int volume = -1;               // le solide de l'aile
    std::vector<int> skin;         // surfaces de peau : extrados, intrados, bord de fuite
    int rootCap = -1;              // face plane a y = 0 (plan de symetrie / encastrement)
    int tipCap = -1;               // face plane a y = b/2 (saumon)
};

// Options du maillage de surface.
struct MeshOptions {
    double sizeMax = 0.10;         // taille maximale des mailles [m]
    double sizeMin = 0.002;        // taille minimale [m]
    int curvaturePoints = 60;      // nb de mailles pour un tour de courbure complet : raffine le BA
    int dim = 2;                   // 2 = maillage de surface (triangles), 3 = volume (tetraedres)
};

// Demarre / arrete Gmsh (verbose = true pour voir les messages de Gmsh).
void initGmsh(bool verbose);
void finalizeGmsh();

// Sections -> B-splines -> loft. Il faut avoir appele initGmsh() avant.
// Le bord de fuite du NACA 4 chiffres est epais (~0,25 % de corde) : on le ferme par un
// segment droit. Un bord de fuite a epaisseur nulle n'est pas gere (exception).
WingCad buildWing(const std::vector<Section>& sections);

// Place les groupes physiques (Skin, Root, Tip, Wing) pour que le .msh soit exploitable.
void defineGroups(const WingCad& cad);

// Maille l'aile.
void meshWing(const MeshOptions& options);

// ---------------------------------------------------------------------------
// Verifications (utilisees par le programme et par le test)
// ---------------------------------------------------------------------------

// Plus grande distance entre un point importe (Phase 2) et la surface CAD.
// Si l'import est fidele, c'est ~0 : la surface passe par tous les points.
double maxDeviation(const std::vector<Section>& sections, const WingCad& cad);

// Volume attendu par la formule (sans Gmsh) :
//   V = a0 * c_r^2 * (b/2) * (1 + lambda + lambda^2) / 3
// avec a0 = aire du profil a corde 1 (formule du lacet / shoelace sur la boucle)
// et lambda = c_t / c_r. Explication : l'aire d'une section vaut a0 * c(y)^2, la fleche,
// le diedre et le vrillage ne changent pas cette aire, et c(y) est lineaire en y.
double expectedVolume(const nacawing::WingParams& w, const std::vector<nacawing::Point2D>& loop);

// Chiffres du maillage. Le volume et l'aire sont calcules a partir des triangles de surface
// (theoreme de la divergence : V = somme de p0 . (p1 x p2) / 6). Ils dependent donc un peu
// de la finesse du maillage (~0,1 % avec les options par defaut).
//
// Pourquoi pas gmsh::model::occ::getMass ? Sur cette aile, il donne un volume faux de 19 %
// (et une face d'emplanture trop grande de 11 %), alors que la geometrie est correcte :
// l'integration d'OpenCASCADE est imprecise sur ces faces B-spline a tres nombreux noeuds.
struct MeshStats {
    std::size_t nodes = 0;
    std::size_t triangles = 0;
    std::size_t tetrahedra = 0;
    double volume = 0.0;      // volume enferme par les surfaces maillees [m^3]
    double skinArea = 0.0;    // aire de la peau (sans les deux faces planes) [m^2]
    double rootArea = 0.0;    // aire de la face d'emplanture [m^2]
};
MeshStats meshStats(const WingCad& cad);

// Ecrit le maillage (.msh) ou la geometrie (.step) ; le format est choisi par l'extension.
void writeFile(const std::string& filename);

}  // namespace incad
