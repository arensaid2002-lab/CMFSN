#include "VLMmesh.h"
#include "nacawing/NACA.hpp"
#include "nacawing/NACASurface.hpp"
#include "nacawing/3Dwing.hpp"
#include <iostream>
#include <fstream> 
#include <algorithm> //std::reverse.

namespace mesh {

    std::vector<nacawing::Point2D> nacaCambrure(double m, double p, int n){     //n = number of points; ni = n-1
        std::vector<nacawing::Point2D> cambrure;
        for (int i = 0; i < n; i++) {   
            double x = nacawing::cosineX(i, n, 1);
            double y_c= nacawing::y_c(m, p, x);
            cambrure.push_back(nacawing::Point2D{x, y_c});
        }
        return cambrure; 
    }

//expects: a closed loop at unit chord in Selig order (TE → upper → LE → lower → TE), like sectionLoop.
    std::vector<nacawing::Point2D> cambrureDuProfil(const std::vector<nacawing::Point2D>& profil){
        int n = static_cast<int>(profil.size());
        std::vector<nacawing::Point2D> cambrure;
        for (int i = 0; i <= (n-1)/2; i++) {
            double x_c = (profil[i].x + profil[n - 1 - i].x) * 0.5;
            double y_c = (profil[i].y + profil[n - 1 - i].y) * 0.5;
            cambrure.push_back(nacawing::Point2D{x_c, y_c});
        }
        std::reverse(cambrure.begin(), cambrure.end());
        return cambrure; 
    }

    VlmGrid vlmGrid(const nacawing::WingParams& w, const std::vector<nacawing::Point2D>& cambrure, int nj, nacawing::SpanSpacing spacing){
        std::vector<std::vector<nacawing::Point3D>>  surface = nacawing::wingSurface(w, cambrure, nj + 1, spacing);
        VlmGrid grid{};
        grid.ni = static_cast<int>(cambrure.size()) - 1;   // 3 points → 2 panels
        grid.nj = nj;                                      // 3 panels (4 stations)
        for (int j = 0; j <= grid.nj; j++){
            for (int i = 0; i <= grid.ni; i++){
                grid.P.push_back(surface[i][j]);
                }
            }       
        return grid;
    }   

     void writeTecplot(const VlmGrid& grid, const std::string& filename){ //location to save the plot .dat,
        std::ofstream file(filename);
        file << "VARIABLES = \"X\" \"Y\" \"Z\"\n"; // line 1
        file << "ZONE I="<< grid.ni + 1 << ", J=" << grid.nj + 1 << ", DATAPACKING=POINT\n"; // line 2
        
        
        for (std::size_t i = 0; i < grid.P.size(); i++){
            file << grid.P[i].x << " " << grid.P[i].y << " " << grid.P[i].z << "\n"; // line 3 to size
            }
    }
}// namespace mesh