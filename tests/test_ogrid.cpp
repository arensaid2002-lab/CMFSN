#include <algorithm>
#include <cmath>
#include <iostream>
#include "OGrid.h"
#include "Plot3D.h"
#include "nacawing/NACA.hpp"

// TEST GÉNÉRÉ AVEC CLAUDE : génère la grille en O NACA0012 pour NC = 8 et la compare au maillage
// de référence x.9. x.9 diffère de notre grille d'environ 1 % à la paroi (décalage d'angle) et
// 0.25 % au loin, donc on compare avec une tolérance relative de 2 %.

double distance(const nacawing::Point2D& a, const nacawing::Point2D& b)
{
    return std::hypot(a.x - b.x, a.y - b.y);
}

int main(int argc, char* argv[])
{
    // Le chemin de x.9 est passé par CMake (voir add_test)
    if (argc < 2) {
        std::cerr << "Usage : test_ogrid <chemin vers x.9>\n";
        return 1;
    }

    const int nc = 8;
    int failures = 0;
    mesh::StructuredMesh2D m = mesh::generateOGrid(nacawing::fromDigits(0, 0, 12), nc);
    mesh::StructuredMesh2D ref = mesh::readPlot3D(argv[1]);

    // 1. Taille : (nc + 1) x (nc + 1), comme x.9
    if (m.ni != nc + 1 || m.nj != nc + 1 || static_cast<int>(m.nodes.size()) != (nc + 1) * (nc + 1)) {
        std::cerr << "Taille : " << m.ni << " x " << m.nj << " (attendu 9 x 9)\n";
        return 1;   // on arrête ici, sinon at() lirait hors du vecteur
    }

    // 2. Chaque nœud proche de x.9 (erreur relative < 2 %, au moins 0.1 de rayon pour la paroi)
    for (int j = 0; j < m.nj; ++j) {
        for (int i = 0; i < m.ni; ++i) {
            double rayon = std::hypot(ref.at(i, j).x, ref.at(i, j).y);
            double erreur = distance(m.at(i, j), ref.at(i, j)) / std::max(rayon, 0.1);
            if (erreur > 0.02) {
                std::cerr << "Nœud (" << i << "," << j << ") : " << m.at(i, j).x << ", " << m.at(i, j).y
                          << " vs x.9 " << ref.at(i, j).x << ", " << ref.at(i, j).y << "\n";
                failures++;
            }
        }
    }

    // 3. Paroi (j = 0) sur le profil : |y| = y_t(x)
    for (int i = 0; i < m.ni; ++i) {
        double x = m.at(i, 0).x;
        if (std::abs(std::abs(m.at(i, 0).y) - nacawing::y_t(0.12, x)) > 1e-6) {
            std::cerr << "Nœud de paroi " << i << " hors du profil\n";
            failures++;
        }
    }

    // 4. Couture : at(0,j) et at(nc,j) sont le même point
    for (int j = 0; j < m.nj; ++j) {
        if (distance(m.at(0, j), m.at(nc, j)) > 1e-9) {
            std::cerr << "Couture incorrecte à j=" << j << "\n";
            failures++;
        }
    }

    // 5. Symétrie haut/bas : (x, y) à i  <->  (x, -y) à nc - i
    for (int j = 0; j < m.nj; ++j) {
        for (int i = 0; i < m.ni; ++i) {
            nacawing::Point2D miroir{m.at(nc - i, j).x, -m.at(nc - i, j).y};
            if (distance(m.at(i, j), miroir) > 1e-9) {
                std::cerr << "Symétrie incorrecte en (" << i << "," << j << ")\n";
                failures++;
            }
        }
    }

    // 6. Aucune cellule retournée : l'aire signée a le même signe partout
    int positives = 0, negatives = 0;
    for (int j = 0; j < m.nj - 1; ++j) {
        for (int i = 0; i < m.ni - 1; ++i) {
            const auto& a = m.at(i, j);
            const auto& b = m.at(i + 1, j);
            const auto& c = m.at(i + 1, j + 1);
            const auto& d = m.at(i, j + 1);
            double aire = 0.5 * ((c.x - a.x) * (d.y - b.y) - (d.x - b.x) * (c.y - a.y));
            if (aire > 0) positives++; else negatives++;
        }
    }
    if (positives != 0 && negatives != 0) {
        std::cerr << "Cellules retournées : " << positives << " positives, " << negatives << " négatives\n";
        failures++;
    }

    std::cout << (failures == 0 ? "OGrid : OK\n" : "OGrid : ECHEC\n");
    return failures == 0 ? 0 : 1;
}
