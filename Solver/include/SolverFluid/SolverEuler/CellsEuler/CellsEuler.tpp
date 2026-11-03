#include "CellsEuler/CellsEuler.h"

namespace Solver
{
        template <typename ExecSpace, typename T>
        Cells_Euler<ExecSpace, T>::Cells_Euler(T Mach_inf, T gamma, T alpha, int N_GhostLayers, const MeshData_Structured<T>& MeshData) : MeshData_(MeshData), CellsData_(N_GhostLayers, MeshData)
        {
            gamma_ = gamma;
            Mach_inf_ = Mach_inf;
            alpha_ = alpha; // Radians
            velocity_u_inf = Mach_inf * sqrt(gamma_) * cos(alpha_);
            velocity_v_inf = Mach_inf * sqrt(gamma_) * sin(alpha_);


            CellOps::eval_GeometricalData<ExecSpace>(CellsData_, MeshData_);
			initialize_ghosts(CellsData_, MeshData_);
        }

        template<typename ExecSpace, typename T>
        void Cells_Euler<ExecSpace, T>::initialize_ghosts(CellsData<T> cellsdata, const MeshData_Structured<T>& meshdata)
        {
            CellOps::ApplyPeriodicity<ExecSpace>(cellsdata, meshdata);
        }
}
