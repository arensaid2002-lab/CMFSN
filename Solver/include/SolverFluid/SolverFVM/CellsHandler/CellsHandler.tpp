#include "CellsHandler/CellsHandler.h"

namespace Solver
{
        template <typename ExecSpace, typename T>
        CellsHandler<ExecSpace, T>::CellsHandler(T Mach_inf, T gamma, T alpha, int N_GhostLayers, const MeshData_Structured<T>& meshdata) : init_data_(Mach_inf, gamma, alpha), MeshData_(meshdata), CellsData_(N_GhostLayers, meshdata)
        {
            CellOps::eval_GeometricalData<ExecSpace>(CellsData_, MeshData_);
			initialize_cells();
			initialize_ghosts();
        }

        template<typename ExecSpace, typename T>
		void CellsHandler<ExecSpace, T>::initialize_cells()
		{
			CellOps::ApplyInitialConditions<ExecSpace>(CellsData_, init_data_);
		}

        template<typename ExecSpace, typename T>
        void CellsHandler<ExecSpace, T>::initialize_ghosts()
        {
            CellOps::eval_DerivedFields<ExecSpace>(CellsData_, init_data_);
            Kokkos::fence();
            CellOps::ApplyPeriodicity<ExecSpace>(CellsData_);
            Kokkos::fence();
            CellOps::eval_DerivedFields<ExecSpace>(CellsData_, init_data_);
            Kokkos::fence();
        }

        
}
