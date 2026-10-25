#include <algorithm>
#include <cmath>
#include <iostream>
#include "OGrid.h"
#include "nacawing/NACA.hpp"

// TEST GÉNÉRÉ AVEC CLAUDE : vérifie que la paroi (anneau j = 0) de la grille en O est toujours
// le profil NACA 00xx, pour plusieurs épaisseurs et tailles de grille. Pas de fichier de référence :
// on compare directement à l'équation d'épaisseur nacawing::y_t de CAD.
// Profils symétriques seulement (m = 0) : la méthode n'est pas encore valable avec cambrure.

int main()
{
    const double x_TE = 1.0089304115;   // bord de fuite : x où y_t(x) = 0 (NACA 4 chiffres, bord ouvert à x = 1)
    const double tol = 1e-9;            // tolérance sur les coordonnées (corde = 1)
    const int thicknesses[] = {6, 12, 15, 24};      // NACA 0006, 0012, 0015, 0024
    const int sizes[] = {8, 16, 32, 64, 128};       // nc, comme dans le viewer
    int failures = 0;

    for (int tt : thicknesses) {
        for (int nc : sizes) {
            const double t = tt / 100.0;    // épaisseur relative (12 -> 0.12)
            mesh::StructuredMesh2D m = mesh::generateOGrid(nacawing::fromDigits(0, 0, tt), nc);
            const int le = nc / 2;          // indice du bord d'attaque
            auto fail = [&](int i, const char* what) {
                std::cerr << "NACA 00" << (tt < 10 ? "0" : "") << tt << ", nc=" << nc << ", nœud " << i
                          << " (" << m.at(i, 0).x << ", " << m.at(i, 0).y << ") : " << what << "\n";
                failures++;
            };

            for (int i = 0; i <= nc; ++i) {
                const nacawing::Point2D& p = m.at(i, 0);   // paroi : j = 0

                // a. Sur le profil : |y| = y_t(x).
                //    y_t contient sqrt(x), donc NaN si x < 0. L'arrondi peut placer le bord d'attaque
                //    à x = -1e-17 : on vérifie d'abord que x est dans [0, x_TE] à 1e-12 près
                //    (sinon c'est une vraie erreur), puis on ramène x dans [0, x_TE] avant d'appeler y_t.
                if (p.x < -1e-12 || p.x > x_TE + 1e-12) {
                    fail(i, "x hors de [0, x_TE]");
                    continue;
                }
                double x = std::clamp(p.x, 0.0, x_TE);
                //    Le bord d'attaque (i = nc/2) est exclu ici : sqrt(x) est vertical en x = 0, donc un
                //    arrondi de 1e-17 sur x donne 1e-9 sur y_t. Ce nœud est vérifié par le test b.
                if (i != le && std::abs(std::abs(p.y) - nacawing::y_t(t, x)) > tol) {
                    fail(i, "hors du profil, |y| != y_t(x)");
                }

                // d. Côté : intrados (i < nc/2) sous l'axe, extrados (i > nc/2) au-dessus
                //    (même convention que isUpper = i > nc/2 dans OGrid.cpp)
                if (i < le && p.y > 0.0) fail(i, "intrados avec y > 0");
                if (i > le && p.y < 0.0) fail(i, "extrados avec y < 0");
            }

            // b. Bord d'attaque : le nœud nc/2 est en (0, 0)
            if (std::hypot(m.at(le, 0).x, m.at(le, 0).y) > tol) fail(le, "bord d'attaque différent de (0, 0)");

            // c. Bord de fuite : les nœuds 0 et nc sont en (x_TE, 0)
            if (std::hypot(m.at(0, 0).x - x_TE, m.at(0, 0).y) > tol) fail(0, "bord de fuite différent de (x_TE, 0)");
            if (std::hypot(m.at(nc, 0).x - x_TE, m.at(nc, 0).y) > tol) fail(nc, "bord de fuite différent de (x_TE, 0)");

            // e. x diminue du bord de fuite au bord d'attaque (i = 0 -> nc/2), puis augmente (nc/2 -> nc)
            for (int i = 0; i < le; ++i) {
                if (m.at(i + 1, 0).x >= m.at(i, 0).x) fail(i + 1, "x ne diminue pas vers le bord d'attaque");
            }
            for (int i = le; i < nc; ++i) {
                if (m.at(i + 1, 0).x <= m.at(i, 0).x) fail(i + 1, "x n'augmente pas vers le bord de fuite");
            }
        }
    }

    std::cout << (failures == 0 ? "OGrid wall : OK\n" : "OGrid wall : ECHEC\n");
    return failures == 0 ? 0 : 1;
}
