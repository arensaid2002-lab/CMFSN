#pragma once

#include "CellsHandler/CellsData.h"
#include "CellOps/CellOps.h"

namespace Solver
{
        template <typename ExecSpace = Kokkos::DefaultExecutionSpace, typename T = float>
        class CellsHandler
        {
            public:
                CellsHandler(T Mach_inf, T gamma, T alpha,int N_GhostLayers, const MeshData_Structured<T>& MeshData);
                ~CellsHandler() = default;

                void initialize_cells();
                void initialize_ghosts(CellsData<T> cellsdata, const MeshData_Structured<T>& meshdata);

			    CellsData<T>& get_CellsData() const { return CellsData_; }
     
            private:
                initialization_data<T> init_data_;
                const MeshData_Structured<T>& MeshData_;
                CellsData<T> CellsData_;
        };
}
#include "CellsHandler/CellsHandler.tpp"
