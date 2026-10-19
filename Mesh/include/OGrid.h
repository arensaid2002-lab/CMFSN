#pragma once
#include "StructuredMesh2D.h"
#include "nacawing/NACA.hpp"


namespace mesh {

    // Grille en O conforme (Vassberg-Jameson) autour d'un NACA 4 chiffres.
    // nc = nombre de cellules par direction, donc ni = nj = nc + 1.
    StructuredMesh2D generateOGrid(const nacawing::NacaParams& naca, int nc);

} // namespace mesh