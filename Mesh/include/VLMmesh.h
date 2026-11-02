#pragma once
#include "nacawing/NACASurface.hpp"
#include "nacawing/3Dwing.hpp"
#include <vector>
#include <string>


namespace mesh {

        std::vector<nacawing::Point2D> nacaCambrure(double m, double p, int n); //t not used reusing y_c from NACA.cpp

        std::vector<nacawing::Point2D> cambrureDuProfil(const std::vector<nacawing::Point2D>& profil); //expects: a closed loop at unit chord in Selig order (TE → upper → LE → lower → TE), like sectionLoop.

        struct VlmGrid { int ni, nj; std::vector<nacawing::Point3D> P; };  // k = j(n_i+1) + i;
        VlmGrid vlmGrid(const nacawing::WingParams& w, const std::vector<nacawing::Point2D>& cambrure, int nj, 
                                                nacawing::SpanSpacing spacing = nacawing::SpanSpacing::SineTip);  //ni is then just cambrure.size() - 1

        void writeTecplot(const VlmGrid& grid, const std::string& filename); //location to save the plot .dat,


} // namespace mesh

