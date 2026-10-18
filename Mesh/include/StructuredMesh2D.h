#pragma once
#include "nacawing/NACASurface.hpp"
#include <vector>

namespace mesh {

        struct StructuredMesh2D {
            int ni = 0; 
            int nj = 0;
            std::vector<nacawing::Point2D> nodes;
            nacawing::Point2D& at(int i, int j) { return nodes[i + j * ni]; }
        };

        

} // namespace mesh

