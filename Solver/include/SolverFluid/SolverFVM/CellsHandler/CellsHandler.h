#pragma once

#include <limits>

#include "CellsHandler/CellsData.h"
#include "CellOps/CellOps.h"

namespace Solver
{
    template <typename ExecSpace = Kokkos::DefaultExecutionSpace, typename T = float, SolverTypes SType = SolverTypes::Euler2D>
    class CellsHandler
    {
        public:
            CellsHandler(T Mach_inf, T gamma, T alpha, int N_ghost_layers, const MeshData_Structured<T>& meshdata);
            CellsHandler(T Mach_inf, T gamma, T alpha, T Re, T Pr, int N_ghost_layers, const MeshData_Structured<T>& meshdata, BoundaryCondition wall_condition = BoundaryCondition::AdiabaticWall, T wall_temperature = T{1});
            ~CellsHandler() = default;

            void initialize_cells();
            void initialize_ghosts();

            CellsData<T>& get_CellsData() { return CellsData_; }
            const CellsData<T>& get_CellsData() const { return CellsData_; }
            const initialization_data<T>& get_InitializationData() const { return init_data_; }

        private:
            CellsData<T> CellsData_;
            initialization_data<T> init_data_;
            const MeshData_Structured<T>& MeshData_;
    };
}
#include "CellsHandler/CellsHandler.tpp"
