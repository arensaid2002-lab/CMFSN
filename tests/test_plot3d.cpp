#include <cmath>
#include <iostream>
#include "Plot3D.h"   

// Les valeurs du fichier ont 6 décimales, donc on compare avec une tolérance.
// TEST GÉNÉRÉ AVEC CLAUDE, COMPARE SI LES VALEURS LUES SONT PROCHES DE CELLES ATTENDUES (1e-6)
bool proche(double a, double b, double tol = 1e-6) { return std::abs(a - b) < tol; }

int main(int argc, char* argv[])
{
    // Le chemin du fichier est passé par CMake (voir add_test)
    if (argc < 2) {
        std::cerr << "Usage : test_plot3d <chemin vers x.9>\n";
        return 1;
    }

    int failures = 0;
    mesh::StructuredMesh2D m = mesh::readPlot3D(argv[1]);

    // 1. En-tête et taille
    if (m.ni != 9 || m.nj != 9 || m.nodes.size() != 81) {
        std::cerr << "Taille lue : " << m.ni << " x " << m.nj << " (attendu 9 x 9)\n";
        return 1;   // on arrête ici, sinon at() lirait hors du vecteur
    }

    // 2. Premier nœud = bord de fuite
    if (!proche(m.at(0, 0).x, 1.008930) || !proche(m.at(0, 0).y, 0.0)) {
        std::cerr << "at(0,0) incorrect\n";
        failures++;
    }

    // 3. at(1,0) ≠ at(0,1) : vérifie qu'on n'a pas inversé i et j
    if (!proche(m.at(1, 0).x, 0.832879) || !proche(m.at(1, 0).y, -0.022517)) {
        std::cerr << "at(1,0) incorrect (i et j inversés ?)\n";
        failures++;
    }

    // 4. Dernier nœud : vérifie qu'on a bien lu tous les X puis tous les Y
    if (!proche(m.at(8, 8).x, 149.042511) || !proche(m.at(8, 8).y, 0.0)) {
        std::cerr << "at(8,8) incorrect (X et Y mélangés ?)\n";
        failures++;
    }

    // 5. Couture : at(0,j) et at(ni-1,j) sont le même point pour chaque j
    for (int j = 0; j < m.nj; ++j) {
        if (!proche(m.at(0, j).x, m.at(m.ni - 1, j).x) || !proche(m.at(0, j).y, m.at(m.ni - 1, j).y)) {
            std::cerr << "Couture incorrecte à j=" << j << "\n";
            failures++;
        }

    }

    std::cout << (failures == 0 ? "Plot3D : OK\n" : "Plot3D : ECHEC\n");
    return failures == 0 ? 0 : 1;
}