#pragma once

#include "MeshFormat.h"

using namespace std;
using namespace Mesh;

using View1D = Kokkos::View<float*, Kokkos::LayoutLeft>;
using View2D = Kokkos::View<float**, Kokkos::LayoutLeft>;
using View3D = Kokkos::View<float***, Kokkos::LayoutLeft>;

namespace Solver
{
    struct CellsData
    {
        View1D Cell_Areas;
        View1D velocity_x, velocity_y, pressure;
        View1D rho_old, Momentum_x_old, Momentum_y_old, E_old;
        View1D rho_new, Momentum_x_new, Momentum_y_new, E_new;
        View2D Flux_mass, Flux_Momentum_x, Flux_Momentum_y, Flux_E, D;
        View3D Normals;

        int N_cells_total{}, Nitotal{}, Njtotal{}, N_gCell_layers{};
    };
}
