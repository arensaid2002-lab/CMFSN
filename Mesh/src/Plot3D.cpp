#include "Plot3D.h"
#include <fstream>
#include <iostream>

namespace mesh {

// readPlot3D : lit un maillage structuré 2D au format Plot3D ASCII (un seul bloc).
//
// IN
//   path : chemin du fichier (ex. ".../NACA_0012_mesh_PLOT3D_format/x.9").
//          Passé en const& : pas de copie, et la fonction ne peut pas le modifier.
//   Contenu attendu du fichier :
//          ni nj
//          ni*nj valeurs de X
//          ni*nj valeurs de Y
//          (>> saute espaces et retours à la ligne : le nombre de valeurs par ligne n'importe pas)
//
// OUT
//   Retourne un StructuredMesh2D (par valeur) :
//     ni, nj : dimensions lues dans l'en-tête
//     nodes  : ni*nj points {x, y}, le point (i, j) est à l'indice k = i + j*ni (voir at(i, j))
//   En cas d'erreur (ouverture, en-tête invalide, fichier tronqué) :
//     un message est écrit dans std::cerr
//     un maillage VIDE est retourné (ni = 0, nj = 0, nodes vide)
//     -> l'appelant teste l'échec avec mesh.nodes.empty()
StructuredMesh2D readPlot3D(const std::string& path)
{
    StructuredMesh2D mesh;

        std::ifstream file(path);       // 1. ouvrir le fichier (std::ifstream)
        if (!file) {
            std::cerr << "Impossible d'ouvrir " << path << "\n";        // 1b. si l'ouverture rate : message dans std::cerr, puis return mesh;
            // si le fichier n'a pas pu être ouvert, on retourne un maillage vide
            return mesh;
        }

        file >> mesh.ni >> mesh.nj;                             // 2. lire ni et nj avec >>
        if (!file || mesh.ni <= 0 || mesh.nj <= 0) {            // 2b. en-tête illisible ou dimensions invalides
            std::cerr << "En-tete invalide dans " << path << " (ni = " << mesh.ni << ", nj = " << mesh.nj << ")\n";
            return StructuredMesh2D{};                          // maillage vide
        }

        mesh.nodes.resize(mesh.ni * mesh.nj);                   // 3. resize de nodes à ni * nj

        for (int k = 0; k < mesh.ni * mesh.nj; ++k) {
            file >> mesh.nodes[k].x;                    // 4a. boucle k : lire nodes[k].x
        }

        for (int k = 0; k < mesh.ni * mesh.nj; ++k) {
            file >> mesh.nodes[k].y;                    // 4b. boucle k : lire nodes[k].y
        }

        if (!file) {                                            // 5. une lecture a raté : fichier tronqué ou valeur non numérique
            std::cerr << "Lecture incomplete de " << path << " (attendu " << 2 * mesh.ni * mesh.nj << " valeurs)\n";
            return StructuredMesh2D{};                          // maillage vide
        }

        return mesh;
}

} // namespace mesh
