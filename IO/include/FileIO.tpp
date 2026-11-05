#pragma once

#include <Kokkos_Core.hpp>
#include <iostream>
#include <fstream>
#include <limits>

#include <MeshHandler.h>

using namespace std;
using namespace Mesh;

namespace IO
{
    template <typename T, MeshFileFormat Format>
    void Import_Mesh(string path, MeshData_Structured<T>& MeshData)
    {
        try
        {
            ifstream Meshfile(path);
            if (Meshfile.is_open())
            {
                MeshData_Structured<T> imported;
                if (!(Meshfile >> imported.Ni >> imported.Nj) || imported.Ni < 2 || imported.Nj < 2)
                    throw runtime_error("Invalid mesh dimensions");
                if (imported.Ni > numeric_limits<int>::max() / imported.Nj)
                    throw runtime_error("Mesh dimensions exceed the supported index range");
                imported.N_nodes = imported.Ni * imported.Nj;
                imported.N_Cells = (imported.Ni-1) * (imported.Nj-1);
 
                Kokkos::View<T*, Kokkos::LayoutLeft > X("X", imported.N_nodes);
                Kokkos::View<T*, Kokkos::LayoutLeft > Y("Y", imported.N_nodes);

                auto X_host = Kokkos::create_mirror_view(X);
                auto Y_host = Kokkos::create_mirror_view(Y);

                for(int i = 0; i < imported.N_nodes; i++)
                {
                    if(!(Meshfile >> X_host(i)))
                    {
                        throw runtime_error("X Coordinate reading failed");
                    }
                }

                for(int i = 0; i < imported.N_nodes; i++)
                {
                    if(!(Meshfile >> Y_host(i)))
                    {
                        throw runtime_error("Y Coordinate reading failed");
                    }
                }
                Kokkos::deep_copy(X, X_host);
                Kokkos::deep_copy(Y, Y_host);

                imported.X = X;
                imported.Y = Y; 
                MeshData = imported;
            }
            else
            {
                throw runtime_error("Error: Mesh file could not be opened");
            }
            Meshfile.close();
        }
        catch (exception& Error)
        {
            cerr << Error.what() << endl;
            throw;
        }
    }
    
}