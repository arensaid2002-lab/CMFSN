#pragma once
#include "nacawing/NACASurface.hpp"
#include "nacawing/3Dwing.hpp"
#include <vector>
#include <string>


namespace mesh {

        std::vector<nacawing::Point2D> nacaCambrurePoint(const nacawing::NacaParams& naca, double x); //x = fraction of chord; uses m and p only (t not used), y_c from NACA.cpp

        struct NastranData {
            int id;
            double x;
            double y;
            double z;
        };

        // The span spacing (Uniform / SineTip) comes from w.spacing.
        std::vector<NastranData> BeamMesh(const nacawing::WingParams& w, const std::vector<nacawing::Point2D>& pointPoutre, int nj);  //ni is then just cambrure.size() - 1

        void writeNastran(const std::vector<NastranData>& grid, const std::string& filename); //location to save the plot .bdf,


} // namespace mesh

