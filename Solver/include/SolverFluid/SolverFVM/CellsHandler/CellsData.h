#pragma once

#include "SolverTypes.h"

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
        static constexpr int N_ghost_layers = 3;

        explicit CellsData(const MeshData_Structured<T>& meshdata) :	Ni_nodes(meshdata.Ni), Nj_nodes(meshdata.Nj), N_faces(meshdata.N_faces), N_physical_cells(meshdata.N_Cells)
        {
            Ni_physical_cells = meshdata.Ni - 1;
            Nj_physical_cells = meshdata.Nj - 1;


            Ni_total_cells         = Ni_physical_cells + 2 * N_ghost_layers;
            Nj_total_cells         = Nj_physical_cells + 2 * N_ghost_layers;
            N_total_cells   = (Ni_total_cells) * (Nj_total_cells);

            Normals     = View3D<T>("Normals", N_total_cells, meshdata.N_faces, 3); // {nx, ny, norm}
            Cell_Areas  = View1D<T>("Cell_Areas", N_total_cells);
			Cell_Centers = View2D<T>("Cell_Centers", N_total_cells, 2); // {x, y}
            Time_Step = View1D<T>("Time_Step", N_total_cells);
             
            
            Pressure = View1D<T>("Pressure", N_total_cells);
            SoundSpeed = View1D<T>("Soundspeed", N_total_cells);
            Lambdas = View2D<T>("Lambdas", N_total_cells, meshdata.N_faces);

            velocity_u      = View1D<T>("velocity_u", N_total_cells);
            velocity_v      = View1D<T>("velocity_v", N_total_cells);
            rho_old         = View1D<T>("rho_old", N_total_cells);
            Momentum_u_old  = View1D<T>("Momentum_u_old", N_total_cells);
            Momentum_v_old  = View1D<T>("Momentum_v_old", N_total_cells);
            E_old           = View1D<T>("E_old", N_total_cells);

            rho_new         = View1D<T>("rho_new", N_total_cells);
            Momentum_u_new  = View1D<T>("Momentum_u_new", N_total_cells);
            Momentum_v_new  = View1D<T>("Momentum_v_new", N_total_cells);
            E_new           = View1D<T>("E_new", N_total_cells);

            rho_0           = View1D<T>("rho_0", N_total_cells);
            Momentum_u_0    = View1D<T>("Momentum_u_0", N_total_cells);
            Momentum_v_0    = View1D<T>("Momentum_v_0", N_total_cells);
            Energy_0        = View1D<T>("Energy_0", N_total_cells);
            D_history       = View2D<T>("D_history", N_total_cells, 4);
            DJacobien        = View3D<T>("dxy_local", N_total_cells, meshdata.N_faces, 4);
            

            Flux_mass       = View2D<T>("Flux_mass", N_total_cells, meshdata.N_faces);
            Flux_Momentum_u = View2D<T>("Flux_Momentum_u", N_total_cells, meshdata.N_faces);
            Flux_Momentum_v = View2D<T>("Flux_Momentum_v", N_total_cells, meshdata.N_faces);
            Flux_E          = View2D<T>("Flux_E", N_total_cells, meshdata.N_faces);


        }
        View1D<T> Cell_Areas, Time_Step;
        View1D<T> velocity_u, velocity_v, Pressure, SoundSpeed;
        View1D<T> rho_old, Momentum_u_old, Momentum_v_old, E_old;
        View1D<T> rho_new, Momentum_u_new, Momentum_v_new, E_new;
        View1D<T> rho_0, Momentum_u_0, Momentum_v_0, Energy_0;
        View2D<T> D_history; // Retained RK contribution per cell, not per face.
        View2D<T> Cell_Centers, Flux_mass, Flux_Momentum_u, Flux_Momentum_v, Flux_E, Lambdas;
        View3D<T> Normals, DJacobien; 
        int Ni_physical_cells, Nj_physical_cells;
        int Ni_nodes, Nj_nodes, N_faces, N_physical_cells;
        int N_total_cells, Ni_total_cells, Nj_total_cells;
    };

    template <typename T = float>
    struct initialization_data
    {
        initialization_data(T Mach_inf, T gamma, T alpha, T Re, T Pr, BoundaryCondition wall_condition_ = BoundaryCondition::AdiabaticWall, T wall_temperature_ = T{1}) : gamma(gamma), Mach_inf(Mach_inf), alpha(alpha), Re_inf(Re), Pr_inf(Pr), wall_condition(wall_condition_), wall_temperature(wall_temperature_)
        {
            c_inf = sqrt(gamma); // alpha is in radians.
            velocity_u_inf = Mach_inf * c_inf * cos(alpha);
            velocity_v_inf = Mach_inf * c_inf * sin(alpha);
			mu_inf = Mach_inf * c_inf / Re_inf;
			k_inf = (gamma/(gamma - 1)) *mu_inf / Pr_inf;
            Momentum_u_inf = rho_inf * velocity_u_inf;
            Momentum_v_inf = rho_inf * velocity_v_inf;
            E_inf = NumOps::eval_Energy<T>(rho_inf, velocity_u_inf, velocity_v_inf, P_inf, gamma);
        }
        T rho_inf = T{ 1.0 }, P_inf = T{ 1.0 }, theta_inf = T{1.0}, mu_inf, k_inf,  c_inf, gamma, Mach_inf, velocity_u_inf, velocity_v_inf, alpha;
        T Momentum_u_inf, Momentum_v_inf, E_inf; // E_inf is total energy per unit volume.
        T Re_inf, Pr_inf;
        BoundaryCondition wall_condition;
        T wall_temperature; // Same nondimensional temperature convention as p/rho.

    };
}
