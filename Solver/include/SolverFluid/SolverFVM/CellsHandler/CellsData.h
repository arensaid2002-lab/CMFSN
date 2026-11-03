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
		CellsData(int N_ghostlayers_, const MeshData_Structured<T>& meshdata) :
			N_i_mesh(meshdata.Ni), N_j_mesh(meshdata.Nj), N_edges(meshdata.N_edges), N_cells_mesh(meshdata.N_Cells)
        {
            N_ghostlayers   = N_ghostlayers_;
            Nitotal         = MeshData.Ni + 2 * N_ghostlayers;
            Njtotal         = MeshData.Nj + 2 * N_ghostlayers;
            N_cells_total   = (Nitotal - 1) * (Njtotal - 1);

            Normals     = View3D<T>("Normals", N_cells_total, MeshData.N_edges, 2);
            Cell_Areas  = View1D<T>("Cell_Areas", N_cells_total);

            velocity_u      = View1D<T>("velocity_u", N_cells_total);
            velocity_v      = View1D<T>("velocity_v", N_cells_total);
            rho_old         = View1D<T>("rho_old", N_cells_total);
            Momentum_u_old  = View1D<T>("Momentum_u_old", N_cells_total);
            Momentum_v_old  = View1D<T>("Momentum_v_old", N_cells_total);
            E_old           = View1D<T>("E_old", N_cells_total);

            rho_new         = View1D<T>("rho_new", N_cells_total);
            Momentum_u_new  = View1D<T>("Momentum_u_new", N_cells_total);
            Momentum_v_new  = View1D<T>("Momentum_v_new", N_cells_total);
            E_new           = View1D<T>("E_new", N_cells_total);

            Flux_mass       = View2D<T>("Flux_mass", N_cells_total, MeshData.N_edges);
            Flux_Momentum_u = View2D<T>("Flux_Momentum_u", N_cells_total, MeshData.N_edges);
            Flux_Momentum_v = View2D<T>("Flux_Momentum_v", N_cells_total, MeshData.N_edges);
            Flux_E          = View2D<T>("Flux_E", N_cells_total, MeshData.N_edges);
            D               = View2D<T>("D", N_cells_total, MeshData.N_edges);
        }
        View1D<T> Cell_Areas;
        View1D<T> velocity_u, velocity_v, pressure;
        View1D<T> rho_old, Momentum_u_old, Momentum_v_old, E_old;
        View1D<T> rho_new, Momentum_u_new, Momentum_v_new, E_new;
        View2D<T> Flux_mass, Flux_Momentum_u, Flux_Momentum_v, Flux_E, D;
        View3D<T> Normals;

        int N_i_mesh, N_j_mesh, N_edges, N_cells_mesh;
        int N_cells_total, Nitotal, Njtotal, N_ghostlayers;
    };

    template <typename T = float>
    struct initialization_data
    {
        initialization_data(T Mach_inf, T gamma, T alpha)
        {
            gamma_ = gamma;
            Mach_inf_ = Mach_inf;
            alpha_ = alpha; // Radians
            velocity_u_inf = Mach_inf * sqrt(gamma_) * cos(alpha_);
            velocity_v_inf = Mach_inf * sqrt(gamma_) * sin(alpha_);
        }
        T rho_inf_{ 1.0f }, P_inf_{ 1.0f }, gamma_, Mach_inf_, velocity_u_inf, velocity_v_inf, alpha_;
    };

}
