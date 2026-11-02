#pragma once

#include "CellsEuler/CellsData.h"
#include "CellOps/CellOps.h"

namespace Solver
{
    template <typename ExecSpace = Kokkos::DefaultExecutionSpace, typename T = float>
    class Cells_Euler
    {
        public:
            Cells_Euler(int N_gCell_layers, const MeshData_Euler<T>& Mesh_Euler);
            ~Cells_Euler() = default;

			void init_Cells_metadata(int N_gCell_layers);

            void init_gcells(CellsData<T> cellsdata, const MeshData_Euler<T>& meshdata);

			CellsData<T> get_CellsData() const { return CellsData_; }
            //CellsData<T> get_CellsData() const {return CellsData_;}
        private:
            
            const MeshData_Euler<T> MeshData;
            CellsData<T> CellsData_;
    };
}
#include "CellsEuler/CellsEuler.tpp"
