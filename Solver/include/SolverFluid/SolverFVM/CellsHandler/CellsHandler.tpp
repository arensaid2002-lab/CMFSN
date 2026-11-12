#pragma once

#include "CellsHandler/CellsHandler.h"

namespace Solver
{
    template <typename ExecSpace, typename T, SolverTypes SType>
    CellsHandler<ExecSpace, T, SType>::CellsHandler(T Mach_inf, T gamma, T alpha, const MeshData_Structured<T>& meshdata) : CellsHandler(Mach_inf, gamma, alpha, std::numeric_limits<T>::infinity(), T{1}, meshdata)
    {
        static_assert(SType == SolverTypes::Euler2D, "RANS2D requires explicit Reynolds and Prandtl numbers");
    }

    template <typename ExecSpace, typename T, SolverTypes SType>
    CellsHandler<ExecSpace, T, SType>::CellsHandler(T Mach_inf, T gamma, T alpha, T Re, T Pr, const MeshData_Structured<T>& meshdata, BoundaryCondition wall_condition, T wall_temperature) : CellsData_(meshdata), init_data_(Mach_inf, gamma, alpha, Re, Pr, wall_condition, wall_temperature), MeshData_(meshdata)
    {
        if (CellsData_.N_faces != 4 ||
            CellsData_.Ni_physical_cells < CellsData<T>::N_ghost_layers ||
            CellsData_.Nj_physical_cells < CellsData<T>::N_ghost_layers)
            throw std::invalid_argument("Cell initialization requires four faces and enough physical cells for the ghost layers");
        if (!(gamma > T{1}) || !std::isfinite(gamma))
            throw std::invalid_argument("Euler initialization requires gamma greater than one");
        if constexpr (SType == SolverTypes::RANS2D)
            if (!(Re > T{0}) || !std::isfinite(Re) || !(Pr > T{0}) || !std::isfinite(Pr))
                throw std::invalid_argument("RANS2D requires positive finite Reynolds and Prandtl numbers");
        if (!(Mach_inf >= T{0}) || !std::isfinite(Mach_inf) || !std::isfinite(alpha))
            throw std::invalid_argument("Mach must be nonnegative and finite, and incidence finite");
        if (wall_condition != BoundaryCondition::AdiabaticWall && wall_condition != BoundaryCondition::isothermalWall)
            throw std::invalid_argument("Unknown wall boundary condition");
        if (wall_condition == BoundaryCondition::isothermalWall && (!(wall_temperature > T{0}) || !std::isfinite(wall_temperature)))
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
