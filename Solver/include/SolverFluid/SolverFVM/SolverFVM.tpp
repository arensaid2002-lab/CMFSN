#pragma once

#include <cmath>
#include <chrono>
#include <iomanip>
#include <limits>
#include "SolverFVM.h"

namespace Solver
{
    template <typename T, SolverTypes SType, TimeEvolutionScheme TEScheme, RungeKuttaScheme RKScheme>
    SolverFVM<T, SType, TEScheme, RKScheme>::SolverFVM(T Mach_inf, T gamma, T alpha, const MeshData_Structured<T>& Mesh_Euler) : SolverFVM(Mach_inf, gamma, alpha, std::numeric_limits<T>::infinity(), T{1}, Mesh_Euler)
    {
        static_assert(SType == SolverTypes::Euler2D, "RANS2D requires explicit Reynolds and Prandtl numbers");
    }

    template <typename T, SolverTypes SType, TimeEvolutionScheme TEScheme, RungeKuttaScheme RKScheme>
    SolverFVM<T, SType, TEScheme, RKScheme>::SolverFVM(T Mach_inf, T gamma, T alpha, T Re, T Pr, const MeshData_Structured<T>& Mesh_Euler, BoundaryCondition wall_condition, T wall_temperature) : cells_handler(Mach_inf, gamma, alpha, Re, Pr, Mesh_Euler, wall_condition, wall_temperature)
    {
    }

    
    template <typename T, SolverTypes SType, TimeEvolutionScheme TEScheme, RungeKuttaScheme RKScheme>
    void SolverFVM<T, SType, TEScheme, RKScheme>::run(int iterations, T CFL, int report_interval, std::ostream* history, StateObserver observer, int snapshot_interval, double final_time)
    {
        if (iterations < 0 || !(CFL > T{0}) || !std::isfinite(CFL) || report_interval < 1 || snapshot_interval < 1)
            throw std::invalid_argument("Solver run requires nonnegative iterations, positive finite CFL, and positive reporting intervals");
        if constexpr (TEScheme == TimeEvolutionScheme::UnSteady)
            if (std::isnan(final_time) || final_time < state_info.physical_time)
                throw std::invalid_argument("Final physical time must not precede the current state");
        const auto start = std::chrono::steady_clock::now();
        auto previous_report = start;
        int previous_iteration = state_info.iteration;
        int last_reported_iteration = state_info.iteration;
        const auto initial = eval_diagnostics(cells_handler.get_CellsData());
        const double initial_rms = initial.residual_rms;
        if (history)
            *history << "iteration,rms,relative_rms,elapsed_seconds,block_iterations,block_seconds,min_density,min_pressure,invalid_cells,physical_time,time_step\n";
        const auto report = [&](const Diagnostics& diagnostics)
        {
            const auto now = std::chrono::steady_clock::now();
            const double elapsed = std::chrono::duration<double>(now - start).count();
            const double block_seconds = std::chrono::duration<double>(now - previous_report).count();
            const double rms = diagnostics.residual_rms;
            const double relative_rms = initial_rms > 0 ? rms / initial_rms : (rms == 0 ? 0 : std::numeric_limits<double>::quiet_NaN());
            const int block_iterations = state_info.iteration - previous_iteration;
            std::cout << "Temps = " << std::fixed << std::setprecision(3) << elapsed << " s | RMS = " << std::scientific << std::setprecision(6) << rms << std::defaultfloat << '\n' << std::flush;
            if (history)
            {
                *history << std::setprecision(17) << state_info.iteration << ',' << rms << ',' << relative_rms << ',' << elapsed << ',' << block_iterations << ',' << block_seconds << ',' << diagnostics.min_density << ',' << diagnostics.min_pressure << ',' << diagnostics.invalid_cells << ',' << state_info.physical_time << ',' << state_info.time_step << '\n' << std::flush;
                if (!*history)
                    throw std::runtime_error("Residual history write failed");
            }
            previous_report = now;
            previous_iteration = state_info.iteration;
            last_reported_iteration = state_info.iteration;
            if (diagnostics.invalid_cells != 0 || !std::isfinite(rms))
                throw std::runtime_error("Invalid physical state or residual; calculation stopped");
        };
        report(initial);
        if (observer && state_info.iteration == 0)
            observer(state_info, cells_handler.get_CellsData());
        for (int iteration = 0; iteration < iterations; ++iteration)
        {
            T step_limit = std::numeric_limits<T>::max();
            if constexpr (TEScheme == TimeEvolutionScheme::UnSteady)
            {
                if (state_info.physical_time >= final_time)
                    break;
                if (std::isfinite(final_time))
                {
                    const double remaining = final_time - state_info.physical_time;
                    step_limit = T(remaining);
                    if (double(step_limit) > remaining)
                        step_limit = std::nextafter(step_limit, T{0});
                    if (!(step_limit > T{0}) || !std::isfinite(step_limit))
                        throw std::runtime_error("Remaining physical duration is outside the scalar precision range");
                }
            }
            if constexpr (TEScheme == TimeEvolutionScheme::Steady)
            {
                const auto diagnostics = eval_diagnostics(cells_handler.get_CellsData());
                if (diagnostics.invalid_cells != 0 || !std::isfinite(diagnostics.residual_rms))
                    throw std::runtime_error("Invalid physical state or residual; calculation stopped");
                if (diagnostics.residual_rms <= 1e-6)
                    break;
            }
            T time_step;
            if constexpr (RKScheme == RungeKuttaScheme::CarpenterKennedy)
                time_step = CellOps::RungeKuttaCarpenterKennedy5<Kokkos::DefaultExecutionSpace, T, SType, TEScheme>(cells_handler.get_CellsData(), cells_handler.get_InitializationData(), CFL, step_limit);
            else
                time_step = CellOps::RungeKutta5<Kokkos::DefaultExecutionSpace, T, SType, TEScheme>(cells_handler.get_CellsData(), cells_handler.get_InitializationData(), CFL, step_limit);
            ++state_info.iteration;
            if constexpr (TEScheme == TimeEvolutionScheme::UnSteady)
            {
                if (!(time_step > T{0}) || !std::isfinite(time_step))
                    throw std::runtime_error("Invalid physical time step");
                state_info.time_step = time_step;
                state_info.physical_time += double(time_step);
            }
            const bool reached_final_time = TEScheme == TimeEvolutionScheme::UnSteady && state_info.physical_time >= final_time;
            if (observer && ((iteration + 1) % snapshot_interval == 0 || iteration == iterations - 1 || reached_final_time))
                observer(state_info, cells_handler.get_CellsData());
            if ((iteration + 1) % report_interval == 0 || iteration == iterations - 1 || reached_final_time)
                report(eval_diagnostics(cells_handler.get_CellsData()));
        }
        if (last_reported_iteration != state_info.iteration)
            report(eval_diagnostics(cells_handler.get_CellsData()));
    }

    template <typename T, SolverTypes SType, TimeEvolutionScheme TEScheme, RungeKuttaScheme RKScheme>
    void SolverFVM<T, SType, TEScheme, RKScheme>::run_until(double final_time, T CFL, int report_interval, std::ostream* history, StateObserver observer, int snapshot_interval, int max_iterations)
    {
        if constexpr (TEScheme != TimeEvolutionScheme::UnSteady)
            throw std::logic_error("Physical-time stopping requires the UnSteady solver mode");
        if (!std::isfinite(final_time) || final_time < state_info.physical_time || max_iterations < 1)
            throw std::invalid_argument("A finite final time and a positive step limit are required");
        run(max_iterations, CFL, report_interval, history, observer, snapshot_interval, final_time);
        if (state_info.physical_time < final_time)
            throw std::runtime_error("Maximum number of physical steps reached before the requested final time");
    }
    template <typename T, SolverTypes SType, TimeEvolutionScheme TEScheme, RungeKuttaScheme RKScheme>
    typename SolverFVM<T, SType, TEScheme, RKScheme>::Diagnostics SolverFVM<T, SType, TEScheme, RKScheme>::eval_diagnostics(const CellsData<T>& cellsdata)
    {
        if (cellsdata.N_physical_cells < 1)
            throw std::invalid_argument("Diagnostics require physical cells");
        Diagnostics result;
        Kokkos::parallel_reduce("Euler2D_StateDiagnostics", Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, cellsdata.N_physical_cells), KOKKOS_LAMBDA(int k_mesh, T& min_density, T& min_pressure, int& invalid)
        {
            const int i_cell = k_mesh % cellsdata.Ni_physical_cells + cellsdata.N_ghost_layers;
            const int j_cell = k_mesh / cellsdata.Ni_physical_cells + cellsdata.N_ghost_layers;
            const int k_cell = CellOps::Cell_index(i_cell, j_cell, cellsdata.Ni_total_cells);
            const T density = cellsdata.rho_old(k_cell), pressure = cellsdata.Pressure(k_cell);
            if (Kokkos::isfinite(density) && density < min_density)
                min_density = density;
            if (Kokkos::isfinite(pressure) && pressure < min_pressure)
                min_pressure = pressure;
            if (!(density > T{0}) || !(pressure > T{0}) || !Kokkos::isfinite(density) || !Kokkos::isfinite(pressure) || !Kokkos::isfinite(cellsdata.Momentum_u_old(k_cell)) || !Kokkos::isfinite(cellsdata.Momentum_v_old(k_cell)) || !Kokkos::isfinite(cellsdata.E_old(k_cell)))
                ++invalid;
        }, Kokkos::Min<T>(result.min_density), Kokkos::Min<T>(result.min_pressure), result.invalid_cells);
        if (result.invalid_cells == 0)
            result.residual_rms = eval_residual_rms(cellsdata);
        return result;
    }
    template <typename T, SolverTypes SType, TimeEvolutionScheme TEScheme, RungeKuttaScheme RKScheme>
    double SolverFVM<T, SType, TEScheme, RKScheme>::eval_residual_rms(const CellsData<T>& cellsdata)
    {
		const auto init_data = cells_handler.get_InitializationData();
		double residual_squared = 0;
		Kokkos::parallel_reduce("Euler2D_ResidualDiagnostic", Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, cellsdata.N_physical_cells), KOKKOS_LAMBDA(int k_mesh, double& sum)
		{
			const int i_cell = k_mesh % cellsdata.Ni_physical_cells + cellsdata.N_ghost_layers;
			const int j_cell = k_mesh / cellsdata.Ni_physical_cells + cellsdata.N_ghost_layers;
			const int k_cell = CellOps::Cell_index(i_cell, j_cell, cellsdata.Ni_total_cells);
			const auto fluxes_dissipation = CellOps::eval_All_FluxesDissipation<T, SType>(cellsdata, init_data, k_cell);
			for (int component = 0; component < 4; ++component)
			{
				const double residual = double(fluxes_dissipation.first[component] - fluxes_dissipation.second[component]) / double(cellsdata.Cell_Areas(k_cell));
				sum += residual * residual;
			}
		}, residual_squared);
		return std::sqrt(residual_squared / (4.0 * cellsdata.N_physical_cells));
    }
}
