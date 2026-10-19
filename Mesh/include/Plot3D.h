#pragma once
#include <string>
#include "StructuredMesh2D.h"

// Fait la lecture d'un fichier Plot3D et retourne un maillage structuré 2D.
namespace mesh {

    StructuredMesh2D readPlot3D(const std::string& path);
        
        }