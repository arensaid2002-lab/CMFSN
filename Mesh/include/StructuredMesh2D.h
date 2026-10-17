#pragma once
#include "nacawing/NACASurface.hpp"
#include <vector>

namespace mesh {

        struct StructuredMesh2D {
            int ni; 
            int nj;
            std::vector<nacawing::Point2D> nodes;
            nacawing::Point2D& at(int i, int j) { return nodes[i + j * ni]; }
        };

        

} // namespace mesh

