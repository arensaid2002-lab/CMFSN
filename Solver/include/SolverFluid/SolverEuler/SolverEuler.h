#pragma once

#include "CellsEuler/CellsEuler.h"

using namespace Mesh;
namespace Solver
{
    template <typename T = float>
    class SolverEuler 
    {
        public:
            SolverEuler(T Mach_inf, T gamma, T alpha, int N_gLayer, const MeshData_Euler<T>& Mesh_Euler);

            void init_cells();

        private:
            Cells_Euler<Kokkos::DefaultExecutionSpace, T> cells_euler;

			T rho_inf_{1.0f}, P_inf_{1.0f}, gamma_, Mach_inf_, velocity_u_inf, velocity_v_inf, alpha_;
    };
}
#include "SolverEuler.tpp"
