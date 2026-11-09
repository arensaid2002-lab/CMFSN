#pragma once

#include "MeshFormat.h"
#include "NumOps/NumFunctions.h"

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
            Time_Step = View1D<T>("Time_Step", N_cells_total);
                
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

            rho_0           = View1D<T>("rho_0", N_cells_total);
            Momentum_u_0    = View1D<T>("Momentum_u_0", N_cells_total);
            Momentum_v_0    = View1D<T>("Momentum_v_0", N_cells_total);
            Energy_0        = View1D<T>("Energy_0", N_cells_total);
            D_history       = View2D<T>("D_history", N_cells_total, 4);

            Flux_mass       = View2D<T>("Flux_mass", N_cells_total, meshdata.N_edges);
            Flux_Momentum_u = View2D<T>("Flux_Momentum_u", N_cells_total, meshdata.N_edges);
            Flux_Momentum_v = View2D<T>("Flux_Momentum_v", N_cells_total, meshdata.N_edges);
            Flux_E          = View2D<T>("Flux_E", N_cells_total, meshdata.N_edges);


        }
        View1D<T> Cell_Areas, Time_Step;
        View1D<T> velocity_u, velocity_v, Pressure, SoundSpeed;
        View1D<T> rho_old, Momentum_u_old, Momentum_v_old, E_old;
        View1D<T> rho_new, Momentum_u_new, Momentum_v_new, E_new;
        View1D<T> rho_0, Momentum_u_0, Momentum_v_0, Energy_0;
        View2D<T> D_history; // Retained RK contribution per cell, not per face.
        View2D<T> Flux_mass, Flux_Momentum_u, Flux_Momentum_v, Flux_E, Lambdas;
        View3D<T> Normals; 
        int Ni_cells, Nj_cells;
        int Ni_mesh, Nj_mesh, N_edges, N_cells_mesh;
        int N_cells_total, Nitotal, Njtotal, N_ghostlayers;
    };

    template <typename T = float>
    struct initialization_data
    {
        initialization_data(T Mach_inf, T gamma, T alpha) : gamma(gamma), Mach_inf(Mach_inf), alpha(alpha)
        {
            c_inf = sqrt(gamma); // alpha is in radians.
            velocity_u_inf = Mach_inf * c_inf * cos(alpha);
            velocity_v_inf = Mach_inf * c_inf * sin(alpha);

            Momentum_u_inf = rho_inf * velocity_u_inf;
            Momentum_v_inf = rho_inf * velocity_v_inf;
            E_inf = NumOps::eval_Energy<T>(rho_inf, velocity_u_inf, velocity_v_inf, P_inf, gamma);
        }
        T rho_inf{ 1.0f }, P_inf{ 1.0f }, c_inf, gamma, Mach_inf, velocity_u_inf, velocity_v_inf, alpha;
        T Momentum_u_inf, Momentum_v_inf, E_inf; // E_inf is total energy per unit volume.
    };

}
