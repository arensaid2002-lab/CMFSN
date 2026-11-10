#pragma once
#include "nacawing/NACASurface.hpp"
#include "nacawing/3Dwing.hpp"
#include <vector>
#include <string>


namespace mesh {

        std::vector<nacawing::Point2D> nacaCambrurePoint(double m, double p, double x); //t not used reusing y_c from NACA.cpp

        struct NastranData {
            int id;
            double x;
            double y;
            double z;
        };

        std::vector<NastranData> BeamMesh(const nacawing::WingParams& w, const std::vector<nacawing::Point2D>& pointPoutre, int nj, 
                                                nacawing::SpanSpacing spacing = nacawing::SpanSpacing::SineTip);  //ni is then just cambrure.size() - 1

        void writeNastran(const std::vector<NastranData>& grid, const std::string& filename); //location to save the plot .bdf,


} // namespace mesh

