#pragma once

#include <ostream>
#include <functional>
#include <limits>
#include "SolverTypes.h"
#include "CellsHandler/CellsHandler.h"

namespace Solver
{
    template <typename T = float, SolverTypes SType = SolverTypes::Euler2D, TimeEvolutionScheme TEScheme = TimeEvolutionScheme::Steady, RungeKuttaScheme RKScheme = RungeKuttaScheme::Unspecified>
    class SolverFVM
    {
        static_assert(SType == SolverTypes::Euler2D || SType == SolverTypes::RANS2D, "Unsupported solver type");
        static_assert(TEScheme != TimeEvolutionScheme::UnSteady || RKScheme != RungeKuttaScheme::Unspecified, "UnSteady mode requires an explicit RungeKuttaScheme");
        static_assert(TEScheme == TimeEvolutionScheme::UnSteady || RKScheme != RungeKuttaScheme::CarpenterKennedy, "Carpenter-Kennedy is only available in UnSteady mode");

        public:
            struct StateInfo { int iteration = 0; double physical_time = 0; T time_step = T{0}; };
            using StateObserver = std::function<void(const StateInfo&, const CellsData<T>&)>;
            const StateInfo& get_StateInfo() const { return state_info; }
            SolverFVM(const FlowConditions<T>& flow, const MeshData_Structured<T>& mesh) : cells_handler(flow, mesh) {}
            // Steady: stops once both the global and the density RMS residuals reach tol (0 disables the stop).
            // report_interval paces the console, history_interval the rows of the residual history.
            void run(int iterations, T tol, T CFL = T{1}, int report_interval = 1000, int history_interval = 1, std::ostream* history = nullptr, StateObserver observer = {}, int snapshot_interval = 1, double final_time = std::numeric_limits<double>::infinity());
            void run_until(double final_time, T CFL = T{1}, int report_interval = 1000, int history_interval = 1, std::ostream* history = nullptr, StateObserver observer = {}, int snapshot_interval = 1, int max_iterations = 200000);
            const CellsData<T>& get_CellsData() const { return cells_handler.get_CellsData(); }
            const initialization_data<T>& get_InitializationData() const { return cells_handler.get_InitializationData(); }
            // residual_rms: all four components; residual_rms_rho: density only (dRHO/dT).
            struct Diagnostics { T min_density{}, min_pressure{}; int invalid_cells{}; double residual_rms = std::numeric_limits<double>::quiet_NaN(), residual_rms_rho = std::numeric_limits<double>::quiet_NaN(); };
            struct Residuals { double rms, rms_rho; };
            Diagnostics eval_diagnostics(const CellsData<T>& cellsdata);
            Residuals eval_residuals(const CellsData<T>& cellsdata);
        private:
            StateInfo state_info;
            CellsHandler<Kokkos::DefaultExecutionSpace, T, SType> cells_handler;
    };
}
#include "SolverFVM.tpp"
