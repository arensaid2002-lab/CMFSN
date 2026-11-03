#pragma once

#include "CellsEuler/CellsData.h"
#include "CellOps/CellOps.h"

namespace Solver
{
    template <typename ExecSpace = Kokkos::DefaultExecutionSpace, typename T = float>
    class Cells_Euler
    {
        public:
            Cells_Euler(T Mach_inf, T gamma, T alpha,int N_GhostLayers, const MeshData_Structured<T>& MeshData);
            ~Cells_Euler() = default;

            void initialize_cells();
            void initialize_ghosts(CellsData<T> cellsdata, const MeshData_Structured<T>& meshdata);

			CellsData<T> get_CellsData() const { return CellsData_; }
     
        private:
            T rho_inf_{ 1.0f }, P_inf_{ 1.0f }, gamma_, Mach_inf_, velocity_u_inf, velocity_v_inf, alpha_;
            const MeshData_Structured<T>& MeshData_;
            CellsData<T> CellsData_;
    };
}
#include "CellsEuler/CellsEuler.tpp"
