#pragma once

#include <Kokkos_Core.hpp>
#include <iostream>
#include <fstream>

#include <MeshHandler.h>

using namespace std;
using namespace Mesh;

namespace IO
{
    template <typename T>
    void Import_Mesh(string path, MeshData_Structured<T>& Mesh_Euler, typename MeshHandler<T>::MeshFileFormat Format)
    {
        try
        {
            fstream Meshfile(path);
            if (Meshfile.is_open())
            {
                Meshfile >> Mesh_Euler.Ni >> Mesh_Euler.Nj;
                Mesh_Euler.N_nodes = Mesh_Euler.Ni * Mesh_Euler.Nj;
                Mesh_Euler.N_Cells = (Mesh_Euler.Ni-1) * (Mesh_Euler.Nj-1);
 
                Kokkos::View<T*, Kokkos::LayoutLeft > X("X", Mesh_Euler.N_nodes);
                Kokkos::View<T*, Kokkos::LayoutLeft > Y("X", Mesh_Euler.N_nodes);

                auto X_host = Kokkos::create_mirror_view(X);
                auto Y_host = Kokkos::create_mirror_view(Y);

                for(int i = 0; i < Mesh_Euler.N_nodes; i++)
                {
                    if(!(Meshfile >> X_host(i)))
                    {
                        throw runtime_error("X Coordinate reading failed");
                    }
                }

                for(int i = 0; i < Mesh_Euler.N_nodes; i++)
                {
                    if(!(Meshfile >> Y_host(i)))
                    {
                        throw runtime_error("Y Coordinate reading failed");
                    }
                }
                Kokkos::deep_copy(X, X_host);
                Kokkos::deep_copy(Y, Y_host);

                Mesh_Euler.X = X;
                Mesh_Euler.Y = Y; 
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