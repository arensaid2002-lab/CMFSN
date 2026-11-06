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
		CellsData(int N_ghostlayers_, const MeshData_Structured<T>& meshdata) :	Ni_mesh(meshdata.Ni), Nj_mesh(meshdata.Nj), N_edges(meshdata.N_edges), N_cells_mesh(meshdata.N_Cells)
        {
            Ni_cells = meshdata.Ni - 1;
            Nj_cells = meshdata.Nj - 1;


            N_ghostlayers   = N_ghostlayers_;
            Nitotal         = meshdata.Ni + 2 * N_ghostlayers;
            Njtotal         = meshdata.Nj + 2 * N_ghostlayers;
            N_cells_total   = (Nitotal - 1) * (Njtotal - 1);

            Normals     = View3D<T>("Normals", N_cells_total, meshdata.N_edges, 3); // {nx, ny, norm}
            Cell_Areas  = View1D<T>("Cell_Areas", N_cells_total);
                
            Pressure = View1D<T>("Pressure", N_cells_total);
            SoundSpeed = View1D<T>("Soundspeed", N_cells_total);
            Lambdas = View2D<T>("Lambdas", N_cells_total, meshdata.N_edges);

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

            Flux_mass       = View2D<T>("Flux_mass", N_cells_total, meshdata.N_edges);
            Flux_Momentum_u = View2D<T>("Flux_Momentum_u", N_cells_total, meshdata.N_edges);
            Flux_Momentum_v = View2D<T>("Flux_Momentum_v", N_cells_total, meshdata.N_edges);
            Flux_E          = View2D<T>("Flux_E", N_cells_total, meshdata.N_edges);
            D               = View3D<T>("D", N_cells_total, meshdata.N_edges, 4);
        }
        View1D<T> Cell_Areas;
        View1D<T> velocity_u, velocity_v, Pressure, SoundSpeed;
        View1D<T> rho_old, Momentum_u_old, Momentum_v_old, E_old;
        View1D<T> rho_new, Momentum_u_new, Momentum_v_new, E_new;
        View2D<T> Flux_mass, Flux_Momentum_u, Flux_Momentum_v, Flux_E, Lambdas;
        View3D<T> Normals, D; // D(cell, face, conservative component)
        int Ni_cells, Nj_cells;
        int Ni_mesh, Nj_mesh, N_edges, N_cells_mesh;
        int N_cells_total, Nitotal, Njtotal, N_ghostlayers;
    };

    template <typename T = float>
    struct initialization_data
    {
        initialization_data(T Mach_inf, T gamma, T alpha)
        {
            gamma_ = gamma;
            c_inf = sqrt(gamma);
            Mach_inf_ = Mach_inf;
            alpha_ = alpha; // Radians
            velocity_u_inf = Mach_inf * sqrt(gamma_) * cos(alpha_);
            velocity_v_inf = Mach_inf * sqrt(gamma_) * sin(alpha_);
        }
        T rho_inf{ 1.0f }, P_inf{ 1.0f }, c_inf, gamma_, Mach_inf_, velocity_u_inf, velocity_v_inf, alpha_;
    };

}
