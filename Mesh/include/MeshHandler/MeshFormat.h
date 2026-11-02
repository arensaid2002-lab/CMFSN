#pragma once

#include <iostream>
#include <Kokkos_Core.hpp>

using namespace std;
namespace Mesh
{
    template <typename T = float>
    struct MeshData_Euler
    {
        int Ni{}, Nj{}, N_nodes{}, N_edges{ 4 }, N_Cells{};

        Kokkos::View<const T*, Kokkos::LayoutLeft> X, Y;
    };
    struct MeshData_VLM
    {
        
    };
}