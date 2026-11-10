#include "BeamMesh.h"
#include "nacawing/NACA.hpp"
#include "nacawing/NACASurface.hpp"
#include "nacawing/3Dwing.hpp"
#include <iostream>
#include <fstream> 
#include <iomanip>

namespace mesh {

    std::vector<nacawing::Point2D> nacaCambrurePoint(double m, double p, double x){     //nx is the chosen point on the cambrure line
        std::vector<nacawing::Point2D> cambrure;
            double y_c= nacawing::y_c(m, p, x);
            cambrure.push_back(nacawing::Point2D{x, y_c});
        return cambrure; 
        }


    std::vector<NastranData> BeamMesh(const nacawing::WingParams& w, const std::vector<nacawing::Point2D>& pointPoutre, int nj, nacawing::SpanSpacing spacing){
        std::vector<std::vector<nacawing::Point3D>>    surface = nacawing::wingSurface(w, pointPoutre, nj + 1, spacing);
        std::vector<NastranData> grid{};

        for (int j = 0; j <= nj; j++){
                int id = j + 1;
                grid.push_back({id, surface[0][j].x ,surface[0][j].y ,surface[0][j].z });
                }
       
        return grid;

    }   

    void writeNastran(const std::vector<NastranData>& grid, const std::string& filename){ //location to save the plot .dat,
        std::ofstream file(filename);
        file << std::setprecision(10);

        file << grid[0].id << "," << grid[0].x << "," << grid[0].y << "," << grid[0].z;
        for (std::size_t i = 1; i < grid.size(); i++){
            file << "," << grid[i].id << "," << grid[i].x << "," << grid[i].y << "," << grid[i].z; // line 3 to size
            }
        file << "\n";
    }

//    void writeNastran(std::vector<NastranData>& grid, const std::string& filename){ //location to save the plot .dat,
//        std::ofstream file(filename);
//        file << "GRID";
//
//        for (int i = 0; i < grid.size(); i++){
//            file << "," << grid[i].id << "," << grid[i].x << "," << grid[i].y << "," << grid[i].z; // line 3 to size
//            }
//        file << "\n";
//     }


} //namespace