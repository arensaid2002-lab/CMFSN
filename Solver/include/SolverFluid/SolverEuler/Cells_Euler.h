#pragma once

#include "CellsData.h"
#include "CellOps/CellOps.h"

namespace Solver
{
    template <typename ExecSpace = Kokkos::DefaultExecutionSpace>
    class Cells_Euler
    {
        public:
            Cells_Euler(int N_gCell_layers, const MeshData_Euler& Mesh_Euler);
            ~Cells_Euler() = default;

			void initialize_Cells_metadata(int N_gCell_layers);

            void initialize_gcells(CellsData cellsdata, const MeshData_Euler& meshdata);

			CellsData get_CellsData() const { return CellsData_; }
            //CellsData get_CellsData() const {return CellsData_;}
        private:
            
            const MeshData_Euler MeshData;
            CellsData CellsData_;
    };
}
#include "Cells_Euler.tpp"
