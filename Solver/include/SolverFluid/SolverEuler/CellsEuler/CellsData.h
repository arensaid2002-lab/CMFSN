#pragma once

#include "MeshFormat.h"

using namespace std;
using namespace Mesh;

template <typename T = float>
using View1D = Kokkos::View<T*, Kokkos::LayoutLeft>;
template <typename T = float>
using View2D = Kokkos::View<T**, Kokkos::LayoutLeft>;
template <typename T = float>
using View3D = Kokkos::View<T***, Kokkos::LayoutLeft>;

namespace Solver
{
    template <typename T = float>
    struct CellsData
    {
        View1D<T> Cell_Areas;
        View1D<T> velocity_u, velocity_v, pressure;
        View1D<T> rho_old, Momentum_x_old, Momentum_y_old, E_old;
        View1D<T> rho_new, Momentum_x_new, Momentum_y_new, E_new;
        View2D<T> Flux_mass, Flux_Momentum_x, Flux_Momentum_y, Flux_E, D;
        View3D<T> Normals;

        int N_cells_total, Nitotal, Njtotal, N_gCell_layers;
    };
}
