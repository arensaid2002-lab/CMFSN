#include "CellsEuler/CellsEuler.h"

namespace Solver
{
        template <typename ExecSpace, typename T>
        Cells_Euler<ExecSpace, T>::Cells_Euler(int N_gCell_layers, const MeshData_Euler<T>& Mesh_Euler) : MeshData(Mesh_Euler)
        {
            init_Cells_metadata(N_gCell_layers);

            CellOps::eval_geometrical_data<ExecSpace>(CellsData_, MeshData);
			init_gcells(CellsData_, MeshData);
        }

        template <typename ExecSpace, typename T>
        void Cells_Euler<ExecSpace, T>::init_Cells_metadata(int N_gCell_layers)
        {
            CellsData_.N_gCell_layers = N_gCell_layers;
			CellsData_.Nitotal = MeshData.Ni + 2 * CellsData_.N_gCell_layers;
			CellsData_.Njtotal = MeshData.Nj + 2 * CellsData_.N_gCell_layers;
            CellsData_.N_cells_total = (CellsData_.Nitotal-1) * (CellsData_.Njtotal-1);

            CellsData_.Normals = View3D<T>("Normals", CellsData_.N_cells_total, MeshData.N_edges, 2);
            CellsData_.Cell_Areas = View1D<T>("Cell_Areas", CellsData_.N_cells_total);

			CellsData_.velocity_u = View1D<T>("velocity_u", CellsData_.N_cells_total);
			CellsData_.velocity_v = View1D<T>("velocity_v", CellsData_.N_cells_total);
            CellsData_.rho_old = View1D<T>("rho_old", CellsData_.N_cells_total);
            CellsData_.Momentum_u_old = View1D<T>("Momentum_u_old", CellsData_.N_cells_total);
            CellsData_.Momentum_v_old = View1D<T>("Momentum_v_old", CellsData_.N_cells_total);
            CellsData_.E_old = View1D<T>("E_old", CellsData_.N_cells_total);

            CellsData_.rho_new = View1D<T>("rho_new", CellsData_.N_cells_total);
            CellsData_.Momentum_u_new = View1D<T>("Momentum_u_new", CellsData_.N_cells_total);
            CellsData_.Momentum_v_new = View1D<T>("Momentum_v_new", CellsData_.N_cells_total);
            CellsData_.E_new = View1D<T>("E_new", CellsData_.N_cells_total);

            CellsData_.Flux_mass = View2D<T>("Flux_mass", CellsData_.N_cells_total, MeshData.N_edges);
            CellsData_.Flux_Momentum_u = View2D<T>("Flux_Momentum_u", CellsData_.N_cells_total, MeshData.N_edges);
            CellsData_.Flux_Momentum_v = View2D<T>("Flux_Momentum_v", CellsData_.N_cells_total, MeshData.N_edges);
            CellsData_.Flux_E = View2D<T>("Flux_E", CellsData_.N_cells_total, MeshData.N_edges);
            CellsData_.D = View2D<T>("D", CellsData_.N_cells_total, MeshData.N_edges);

        }

        template<typename ExecSpace, typename T>
        void Cells_Euler<ExecSpace, T>::init_gcells(CellsData<T> cellsdata, const MeshData_Euler<T>& meshdata)
        {
            CellOps::ApplyPeriodicity<ExecSpace>(cellsdata, meshdata);
        }
}
