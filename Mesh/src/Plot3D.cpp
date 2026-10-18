#include "Plot3D.h"
#include <fstream>
#include <iostream>

namespace mesh {

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
        mesh.nodes.resize(mesh.ni * mesh.nj);                   // 3. resize de nodes à ni * nj

        for (int k = 0; k < mesh.ni * mesh.nj; ++k) {
            file >> mesh.nodes[k].x;                    // 4a. boucle k : lire nodes[k].x
        }

        for (int k = 0; k < mesh.ni * mesh.nj; ++k) {
            file >> mesh.nodes[k].y;                    // 4b. boucle k : lire nodes[k].y
        }

        return mesh;            //
}

} // namespace mesh