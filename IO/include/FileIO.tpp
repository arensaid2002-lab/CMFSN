#pragma once

#include <Kokkos_Core.hpp>
#include <iostream>
#include <fstream>

#include <MeshHandler.h>

using namespace std;
using namespace Mesh;

namespace IO
{
    template <typename T, typename MeshFileFormat>
    void Import_Mesh(string path, MeshData_Structured<T>& MeshData)
    {
        try
        {
            fstream Meshfile(path);
            if (Meshfile.is_open())
            {
                Meshfile >> MeshData.Ni >> MeshData.Nj;
                MeshData.N_nodes = MeshData.Ni * MeshData.Nj;
                MeshData.N_Cells = (MeshData.Ni-1) * (MeshData.Nj-1);
 
                Kokkos::View<T*, Kokkos::LayoutLeft > X("X", MeshData.N_nodes);
                Kokkos::View<T*, Kokkos::LayoutLeft > Y("X", MeshData.N_nodes);

                auto X_host = Kokkos::create_mirror_view(X);
                auto Y_host = Kokkos::create_mirror_view(Y);

                for(int i = 0; i < MeshData.N_nodes; i++)
                {
                    if(!(Meshfile >> X_host(i)))
                    {
                        throw runtime_error("X Coordinate reading failed");
                    }
                }

                for(int i = 0; i < MeshData.N_nodes; i++)
                {
                    if(!(Meshfile >> Y_host(i)))
                    {
                        throw runtime_error("Y Coordinate reading failed");
                    }
                }
                Kokkos::deep_copy(X, X_host);
                Kokkos::deep_copy(Y, Y_host);

                MeshData.X = X;
                MeshData.Y = Y; 
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
        }
    }
    
}