#include "Cells_Euler.h"

namespace Solver
{
        template <typename ExecSpace>
        Cells_Euler<ExecSpace>::Cells_Euler(int N_gCell_layers, const MeshData_Euler& Mesh_Euler) : MeshData(Mesh_Euler)
        {
            initialize_Cells_metadata(N_gCell_layers);

            CellOps::eval_geometrical_data<ExecSpace>(CellsData_, MeshData);
			initialize_gcells(CellsData_, MeshData);
        }

        template <typename ExecSpace>
        void Cells_Euler<ExecSpace>::initialize_Cells_metadata(int N_gCell_layers)
        {
            CellsData_.N_gCell_layers = N_gCell_layers;
			CellsData_.Nitotal = MeshData.Ni + 2 * CellsData_.N_gCell_layers;
			CellsData_.Njtotal = MeshData.Nj + 2 * CellsData_.N_gCell_layers;
            CellsData_.N_cells_total = (CellsData_.Nitotal-1) * (CellsData_.Njtotal-1);

            CellsData_.Normals = View3D("Normals", CellsData_.N_cells_total, MeshData.N_edges, 2);
            CellsData_.Cell_Areas = View1D("Cell_Areas", CellsData_.N_cells_total);

            CellsData_.rho_old = View1D("rho_old", CellsData_.N_cells_total);
            CellsData_.Momentum_x_old = View1D("Momentum_x_old", CellsData_.N_cells_total);
            CellsData_.Momentum_y_old = View1D("Momentum_y_old", CellsData_.N_cells_total);
            CellsData_.E_old = View1D("E_old", CellsData_.N_cells_total);

            CellsData_.rho_new = View1D("rho_new", CellsData_.N_cells_total);
            CellsData_.Momentum_x_new = View1D("Momentum_x_new", CellsData_.N_cells_total);
            CellsData_.Momentum_y_new = View1D("Momentum_y_new", CellsData_.N_cells_total);
            CellsData_.E_new = View1D("E_new", CellsData_.N_cells_total);

            CellsData_.Flux_mass = View2D("Flux_mass", CellsData_.N_cells_total, MeshData.N_edges);
            CellsData_.Flux_Momentum_x = View2D("Flux_Momentum_x", CellsData_.N_cells_total, MeshData.N_edges);
            CellsData_.Flux_Momentum_y = View2D("Flux_Momentum_y", CellsData_.N_cells_total, MeshData.N_edges);
            CellsData_.Flux_E = View2D("Flux_E", CellsData_.N_cells_total, MeshData.N_edges);
            CellsData_.D = View2D("D", CellsData_.N_cells_total, MeshData.N_edges);

        }

        template<typename ExecSpace>
        void Cells_Euler<ExecSpace>::initialize_gcells(CellsData cellsdata, const MeshData_Euler& meshdata)
        {
            CellOps::ApplyPeriodicity<ExecSpace>(cellsdata, meshdata);
        }
}
