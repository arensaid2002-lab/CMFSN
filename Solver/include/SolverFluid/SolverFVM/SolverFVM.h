#pragma once

#include "CellsHandler/CellsHandler.h"

using namespace Mesh;
using ExecSpace = Kokkos::DefaultExecutionSpace;
namespace Solver
{
	enum class SolverTypes
    {
        Euler2D,
        RANS2D,
    };
    template <typename T = float, SolverTypes SType = SolverTypes::Euler2D >
    class SolverFVM 
    {
        public:
            SolverFVM(T Mach_inf, T gamma, T alpha, int N_gLayer, const MeshData_Structured<T>& Mesh_Euler);

            void init_cells();



            void eval_dissipation();

        private:
            CellsHandler<ExecSpace, T> cells_handler;

    };
}
#include "SolverFVM.tpp"
