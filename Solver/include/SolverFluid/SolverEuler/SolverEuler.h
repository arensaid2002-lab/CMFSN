#pragma once

#include "CellsEuler/CellsEuler.h"

using namespace Mesh;
using ExecSpace = Kokkos::DefaultExecutionSpace;
namespace Solver
{
	enum class TraversalPolicy
    {
        Flat1D,
        MD2D
    };
    template <typename T = float, TraversalPolicy Policy = TraversalPolicy::MD2D >
    class SolverEuler 
    {
        public:
            SolverEuler(T Mach_inf, T gamma, T alpha, int N_gLayer, const MeshData_Euler<T>& Mesh_Euler);

            void init_cells();



            void eval_dissipation()

        private:
            Cells_Euler<ExecSpace, T> cells_euler;

			T rho_inf_{1.0f}, P_inf_{1.0f}, gamma_, Mach_inf_, velocity_u_inf, velocity_v_inf, alpha_;
    };
}
#include "SolverEuler.tpp"
