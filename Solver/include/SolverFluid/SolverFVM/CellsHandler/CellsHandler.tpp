#pragma once

#include <cmath>
#include <stdexcept>
#include "CellsHandler/CellsHandler.h"

namespace Solver
{
    // Euler2D ignores the Reynolds number: the stored conditions are inviscid (mu_inf = k_inf = 0).
    template <typename T>
    FlowConditions<T> inviscid(FlowConditions<T> flow)
    {
        flow.Re_inf = std::numeric_limits<T>::infinity();
        return flow;
    }

    template <typename ExecSpace, typename T, SolverTypes SType>
    CellsHandler<ExecSpace, T, SType>::CellsHandler(const FlowConditions<T>& flow, const MeshData_Structured<T>& meshdata) : CellsData_(meshdata), init_data_(SType == SolverTypes::Euler2D ? inviscid(flow) : flow), MeshData_(meshdata)
    {
        const auto positive = [](T value) { return value > T{0} && Kokkos::isfinite(value); };
        if (CellsData_.N_faces != 4 ||
            CellsData_.Ni_physical_cells < CellsData<T>::N_ghost_layers ||
            CellsData_.Nj_physical_cells < CellsData<T>::N_ghost_layers)
            throw std::invalid_argument("Cell initialization requires four faces and enough physical cells for the ghost layers");
        if (!(flow.gamma > T{1}) || !Kokkos::isfinite(flow.gamma))
            throw std::invalid_argument("Euler initialization requires gamma greater than one");
        if (!(flow.Mach_inf >= T{0}) || !Kokkos::isfinite(flow.Mach_inf) || !Kokkos::isfinite(flow.alpha))
            throw std::invalid_argument("Mach must be nonnegative and finite, and incidence finite");
        if constexpr (SType == SolverTypes::RANS2D)
            if (!positive(flow.Re_inf) || !positive(flow.Pr_inf))
                throw std::invalid_argument("RANS2D requires positive finite Reynolds and Prandtl numbers");
        if (flow.wall_condition != BoundaryCondition::AdiabaticWall && flow.wall_condition != BoundaryCondition::isothermalWall)
            throw std::invalid_argument("Unknown wall boundary condition");
        if (flow.wall_condition == BoundaryCondition::isothermalWall && !positive(flow.wall_temperature))
            throw std::invalid_argument("Isothermal wall temperature must be positive, finite and nondimensional");
        CellOps::eval_GeometricalData<ExecSpace>(CellsData_, MeshData_);
        ExecSpace().fence("Cell geometry ready");
        initialize_cells();
        initialize_ghosts();
    }

    template <typename ExecSpace, typename T, SolverTypes SType>
    void CellsHandler<ExecSpace, T, SType>::initialize_cells()
    {
        const auto cellsdata = CellsData_;
        CellOps::ApplyInitialConditions<ExecSpace>(cellsdata, init_data_);
        ExecSpace().fence("Initial conservative and derived fields ready");
        Kokkos::parallel_for("InitializeCellFluxes", Kokkos::RangePolicy<ExecSpace>(0, cellsdata.N_physical_cells), KOKKOS_LAMBDA(int k_mesh)
        {
            const int i_cell = k_mesh % cellsdata.Ni_physical_cells + cellsdata.N_ghost_layers;
            const int j_cell = k_mesh / cellsdata.Ni_physical_cells + cellsdata.N_ghost_layers;
            const int k_cell = CellOps::Cell_index(i_cell, j_cell, cellsdata.Ni_total_cells);
            for (int direction = 0; direction < cellsdata.N_faces; ++direction)
                CellOps::Update_DataField<T, CellOps::DataFields::Fluxes>(cellsdata, k_cell, direction);
        });
        ExecSpace().fence("Initial physical fluxes ready");
    }

    template <typename ExecSpace, typename T, SolverTypes SType>
    void CellsHandler<ExecSpace, T, SType>::initialize_ghosts()
    {
        CellOps::ApplyBoundaryConditions<ExecSpace, T, SType>(CellsData_, init_data_);
        ExecSpace().fence("Initial ghost fields ready");
    }
}
