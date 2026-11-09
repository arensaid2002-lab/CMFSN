#pragma once

#include "CellsHandler/CellsHandler.h"

namespace Solver
{
    template <typename ExecSpace, typename T>
    CellsHandler<ExecSpace, T>::CellsHandler(T Mach_inf, T gamma, T alpha, int N_GhostLayers, const MeshData_Structured<T>& meshdata) : CellsData_(N_GhostLayers, meshdata), init_data_(Mach_inf, gamma, alpha), MeshData_(meshdata)
    {
        if (CellsData_.N_edges != 4 || N_GhostLayers < 1 || CellsData_.Ni_cells < N_GhostLayers || CellsData_.Nj_cells < N_GhostLayers)
            throw std::invalid_argument("Cell initialization requires four faces and enough physical cells for the ghost layers");
        if (!(gamma > T{1}))
            throw std::invalid_argument("Euler initialization requires gamma greater than one");
        CellOps::eval_GeometricalData<ExecSpace>(CellsData_, MeshData_);
        ExecSpace().fence("Cell geometry ready");
        initialize_cells();
        initialize_ghosts();
    }

    template <typename ExecSpace, typename T>
    void CellsHandler<ExecSpace, T>::initialize_cells()
    {
        const auto cellsdata = CellsData_;
        CellOps::ApplyInitialConditions<ExecSpace>(cellsdata, init_data_);
        ExecSpace().fence("Initial conservative and derived fields ready");
        Kokkos::parallel_for("InitializeCellFluxes", Kokkos::RangePolicy<ExecSpace>(0, cellsdata.N_cells_mesh), KOKKOS_LAMBDA(int k_mesh)
        {
            const int i_cell = k_mesh % cellsdata.Ni_cells + cellsdata.N_ghostlayers;
            const int j_cell = k_mesh / cellsdata.Ni_cells + cellsdata.N_ghostlayers;
            const int k_cell = CellOps::Cell_index(i_cell, j_cell, cellsdata.Nitotal);
            for (int direction = 0; direction < cellsdata.N_edges; ++direction)
                CellOps::Update_DataField<T, CellOps::DataFields::Fluxes>(cellsdata, k_cell, direction);
        });
        ExecSpace().fence("Initial physical fluxes ready");
    }

    template <typename ExecSpace, typename T>
    void CellsHandler<ExecSpace, T>::initialize_ghosts()
    {
        CellOps::ApplyPeriodicity<ExecSpace>(CellsData_);
        CellOps::ApplyNonPenetrationConditions<ExecSpace>(CellsData_, init_data_);
        CellOps::ApplyRiemanInvariantConditions<ExecSpace>(CellsData_, init_data_);
        ExecSpace().fence("Initial ghost fields ready");
    }
}
