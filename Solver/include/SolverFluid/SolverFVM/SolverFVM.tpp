#pragma once

#include <cmath>
#include <chrono>
#include <iomanip>
#include <limits>
#include "SolverFVM.h"

namespace Solver
{
    template <typename T, SolverTypes SType, TimeEvolutionScheme TEScheme, RungeKuttaScheme RKScheme>
    void SolverFVM<T, SType, TEScheme, RKScheme>::run(int iterations, T tol, T CFL, int report_interval, int history_interval, std::ostream* history, StateObserver observer, int snapshot_interval, double final_time)
    {
        if (iterations < 0 || !(CFL > T{0}) || !Kokkos::isfinite(CFL) || report_interval < 1 || history_interval < 1 || snapshot_interval < 1)
            throw std::invalid_argument("Solver run requires nonnegative iterations, positive finite CFL, and positive reporting intervals");
        // tol = 0 disables the residual stop; a NaN or negative value would disable it silently.
        if (!(tol >= T{0}) || !Kokkos::isfinite(tol))
            throw std::invalid_argument("Solver run requires a nonnegative finite residual tolerance");
        if constexpr (TEScheme == TimeEvolutionScheme::UnSteady)
            if (Kokkos::isnan(final_time) || final_time < state_info.physical_time)
                throw std::invalid_argument("Final physical time must not precede the current state");
        const auto start = std::chrono::steady_clock::now();
        auto previous_record = start;
        int previous_iteration = state_info.iteration;
        int last_console_iteration = -1, last_history_iteration = -1;
        Diagnostics current = eval_diagnostics(cells_handler.get_CellsData());
        const double initial_rms = current.residual_rms;
        if (history)
            *history << "iteration,rms,rms_rho,relative_rms,elapsed_seconds,block_iterations,block_seconds,min_density,min_pressure,invalid_cells,physical_time,time_step\n";
        // Console line every report_interval iterations, history row every history_interval iterations.
        const auto report = [&](const Diagnostics& diagnostics, bool to_console, bool to_history)
        {
            const auto now = std::chrono::steady_clock::now();
            const double elapsed = std::chrono::duration<double>(now - start).count();
            const double rms = diagnostics.residual_rms;
            if (to_console)
            {
                std::cout << "Iteration " << state_info.iteration << " | Temps = " << std::fixed << std::setprecision(3) << elapsed << " s | RMS = " << std::scientific << std::setprecision(6) << rms << " | RMS rho = " << diagnostics.residual_rms_rho << std::defaultfloat << '\n' << std::flush;
                last_console_iteration = state_info.iteration;
            }
            if (to_history && history)
            {
                const double block_seconds = std::chrono::duration<double>(now - previous_record).count();
                const double relative_rms = initial_rms > 0 ? rms / initial_rms : (rms == 0 ? 0 : std::numeric_limits<double>::quiet_NaN());
                const int block_iterations = state_info.iteration - previous_iteration;
                *history << std::setprecision(17) << state_info.iteration << ',' << rms << ',' << diagnostics.residual_rms_rho << ',' << relative_rms << ',' << elapsed << ',' << block_iterations << ',' << block_seconds << ',' << diagnostics.min_density << ',' << diagnostics.min_pressure << ',' << diagnostics.invalid_cells << ',' << state_info.physical_time << ',' << state_info.time_step << '\n' << std::flush;
                if (!*history)
                    throw std::runtime_error("Residual history write failed");
                previous_record = now;
                previous_iteration = state_info.iteration;
                last_history_iteration = state_info.iteration;
            }
            if (diagnostics.invalid_cells != 0 || !Kokkos::isfinite(rms) || !Kokkos::isfinite(diagnostics.residual_rms_rho))
                throw std::runtime_error("Invalid physical state or residual; calculation stopped");
        };
        report(current, true, history != nullptr);
        if (observer && state_info.iteration == 0)
            observer(state_info, cells_handler.get_CellsData());
        for (int iteration = 0; iteration < iterations; ++iteration)
        {
            T step_limit = std::numeric_limits<T>::max();
            if constexpr (TEScheme == TimeEvolutionScheme::UnSteady)
            {
                if (state_info.physical_time >= final_time)
                    break;
                if (Kokkos::isfinite(final_time))
                {
                    const double remaining = final_time - state_info.physical_time;
                    step_limit = T(remaining);
                    if (double(step_limit) > remaining)
                        step_limit = Kokkos::nextafter(step_limit, T{0});
                    if (!(step_limit > T{0}) || !Kokkos::isfinite(step_limit))
                        throw std::runtime_error("Remaining physical duration is outside the scalar precision range");
                }
            }
            // Steady stop: both the global and the density RMS residuals must reach tol.
            if constexpr (TEScheme == TimeEvolutionScheme::Steady)
                if (current.residual_rms <= double(tol) && current.residual_rms_rho <= double(tol))
                    break;
            T time_step;
            if constexpr (RKScheme == RungeKuttaScheme::CarpenterKennedy)
                time_step = CellOps::RungeKuttaCarpenterKennedy5<Kokkos::DefaultExecutionSpace, T, SType, TEScheme>(cells_handler.get_CellsData(), cells_handler.get_InitializationData(), CFL, step_limit);
            else
                time_step = CellOps::RungeKutta5<Kokkos::DefaultExecutionSpace, T, SType, TEScheme>(cells_handler.get_CellsData(), cells_handler.get_InitializationData(), CFL, step_limit);
            ++state_info.iteration;
            if constexpr (TEScheme == TimeEvolutionScheme::UnSteady)
            {
                if (!(time_step > T{0}) || !Kokkos::isfinite(time_step))
                    throw std::runtime_error("Invalid physical time step");
                state_info.time_step = time_step;
                state_info.physical_time += double(time_step);
            }
            const bool reached_final_time = TEScheme == TimeEvolutionScheme::UnSteady && state_info.physical_time >= final_time;
            const bool last = iteration == iterations - 1 || reached_final_time;
            if (observer && ((iteration + 1) % snapshot_interval == 0 || last))
                observer(state_info, cells_handler.get_CellsData());
            const bool to_console = (iteration + 1) % report_interval == 0 || last;
            const bool to_history = history && ((iteration + 1) % history_interval == 0 || last);
            if constexpr (TEScheme == TimeEvolutionScheme::Steady)
            {
                // Evaluated once per iteration: it feeds the history and the next stopping test.
                current = eval_diagnostics(cells_handler.get_CellsData());
                if (current.invalid_cells != 0 || !Kokkos::isfinite(current.residual_rms) || !Kokkos::isfinite(current.residual_rms_rho))
                    throw std::runtime_error("Invalid physical state or residual; calculation stopped");
                if (to_console || to_history)
                    report(current, to_console, to_history);
            }
            else if (to_console || to_history)
            {
                current = eval_diagnostics(cells_handler.get_CellsData());
                report(current, to_console, to_history);
            }
        }
        // The last state is always reported, on the console and in the history.
        const bool console_missing = last_console_iteration != state_info.iteration;
        const bool history_missing = history && last_history_iteration != state_info.iteration;
        if (console_missing || history_missing)
        {
            if constexpr (TEScheme == TimeEvolutionScheme::UnSteady)
                current = eval_diagnostics(cells_handler.get_CellsData());
            report(current, console_missing, history_missing);
        }
    }

    template <typename T, SolverTypes SType, TimeEvolutionScheme TEScheme, RungeKuttaScheme RKScheme>
    void SolverFVM<T, SType, TEScheme, RKScheme>::run_until(double final_time, T CFL, int report_interval, int history_interval, std::ostream* history, StateObserver observer, int snapshot_interval, int max_iterations)
    {
        if constexpr (TEScheme != TimeEvolutionScheme::UnSteady)
            throw std::logic_error("Physical-time stopping requires the UnSteady solver mode");
        if (!Kokkos::isfinite(final_time) || final_time < state_info.physical_time || max_iterations < 1)
            throw std::invalid_argument("A finite final time and a positive step limit are required");
        run(max_iterations, T{0}, CFL, report_interval, history_interval, history, observer, snapshot_interval, final_time); // No residual stop in UnSteady.
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
        {
            const auto residuals = eval_residuals(cellsdata);
            result.residual_rms = residuals.rms;
            result.residual_rms_rho = residuals.rms_rho;
        }
        return result;
    }
    // RMS over the physical cells of R/A = -dW/dt (local time step): all four components, and the density alone (dRHO/dT).
    template <typename T, SolverTypes SType, TimeEvolutionScheme TEScheme, RungeKuttaScheme RKScheme>
    typename SolverFVM<T, SType, TEScheme, RKScheme>::Residuals SolverFVM<T, SType, TEScheme, RKScheme>::eval_residuals(const CellsData<T>& cellsdata)
    {
        const auto init_data = cells_handler.get_InitializationData();
        double residual_squared = 0, density_squared = 0;
        Kokkos::parallel_reduce("Euler2D_ResidualDiagnostic", Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, cellsdata.N_physical_cells), KOKKOS_LAMBDA(int k_mesh, double& sum, double& sum_rho)
        {
            const int i_cell = k_mesh % cellsdata.Ni_physical_cells + cellsdata.N_ghost_layers;
            const int j_cell = k_mesh / cellsdata.Ni_physical_cells + cellsdata.N_ghost_layers;
            const int k_cell = CellOps::Cell_index(i_cell, j_cell, cellsdata.Ni_total_cells);
            const auto fluxes_dissipation = CellOps::eval_All_FluxesDissipation<T, SType>(cellsdata, init_data, k_cell);
            for (int component = 0; component < 4; ++component)
            {
                const double residual = double(fluxes_dissipation.first[component] - fluxes_dissipation.second[component]) / double(cellsdata.Cell_Areas(k_cell));
                sum += residual * residual;
                if (component == 0)
                    sum_rho += residual * residual;
            }
        }, residual_squared, density_squared);
        return { Kokkos::sqrt(residual_squared / (4.0 * cellsdata.N_physical_cells)), Kokkos::sqrt(density_squared / double(cellsdata.N_physical_cells)) };
    }
}
