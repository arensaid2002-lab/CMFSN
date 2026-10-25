#pragma once

#include <iostream>
#include <Kokkos_Core.hpp>

using namespace std;
namespace Mesh
{
    struct MeshData_Euler
    {
        int Ni{}, Nj{}, N_nodes{}, N_edge{4}, Mesh_N_Cells{};

        Kokkos::View<float* > X, Y;
    };
    struct MeshData_VLM
    {
        
    };
}