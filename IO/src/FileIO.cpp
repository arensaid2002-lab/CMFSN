#include <Kokkos_Core.hpp>
#include <iostream>
#include <fstream>

#include <MeshHandler.h>

using namespace std;
using namespace Mesh;

namespace IO
{
    void Import_Mesh(string path, MeshData_Euler& Mesh_Euler, MeshHandler::MeshFileFormat Format)
    {
        try
        {
            fstream Meshfile(path);
            if (Meshfile.is_open())
            {
                Meshfile >> Mesh_Euler.Ni >> Mesh_Euler.Nj;
                Mesh_Euler.N_nodes = Mesh_Euler.Ni * Mesh_Euler.Nj;
                
 
                Mesh_Euler.X = Kokkos::View<float* >("X", Mesh_Euler.N_nodes);
                Mesh_Euler.Y = Kokkos::View<float* >("X", Mesh_Euler.N_nodes);

                auto X_host = Kokkos::create_mirror_view(Mesh_Euler.X);
                auto Y_host = Kokkos::create_mirror_view(Mesh_Euler.Y);

                for(int i = 0; i < Mesh_Euler.N_nodes; i++)
                {
                    if(!(Meshfile >> Mesh_Euler.X(i)))
                    {
                        throw runtime_error("X Coordinate reading failed");
                    }
                }

                for(int i = 0; i < Mesh_Euler.N_nodes; i++)
                {
                    if(!(Meshfile >> Mesh_Euler.Y(i)))
                    {
                        throw runtime_error("Y Coordinate reading failed");
                    }
                }
                Kokkos::deep_copy(Mesh_Euler.X, X_host);
                Kokkos::deep_copy(Mesh_Euler.Y, Y_host);
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