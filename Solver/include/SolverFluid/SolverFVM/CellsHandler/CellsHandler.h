#pragma once

#include "CellsHandler/CellsData.h"
#include "CellOps/CellOps.h"

namespace Solver
{
    template <typename ExecSpace = Kokkos::DefaultExecutionSpace, typename T = float>
    class CellsHandler
    {
        public:
            CellsHandler(T Mach_inf, T gamma, T alpha, int N_GhostLayers, const MeshData_Structured<T>& meshdata);
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
