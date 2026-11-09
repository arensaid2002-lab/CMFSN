#pragma once

#include <ostream>
#include <functional>
#include <limits>
#include "SolverTypes.h"
#include "CellsHandler/CellsHandler.h"

namespace Solver
{
    template <typename T = float, SolverTypes SType = SolverTypes::Euler2D, TimeEvolutionScheme TEScheme = TimeEvolutionScheme::Steady>
    class SolverFVM
    {
        static_assert(SType == SolverTypes::Euler2D, "Only the Euler2D solver is implemented");

        public:
            struct StateInfo { int iteration = 0; double physical_time = 0; T time_step = T{0}; };
            using StateObserver = std::function<void(const StateInfo&, const CellsData<T>&)>;
            const StateInfo& get_StateInfo() const { return state_info; }
            SolverFVM(T Mach_inf, T gamma, T alpha, int N_gLayer, const MeshData_Structured<T>& Mesh_Euler);
            void run(int iterations, T CFL = T{1}, int report_interval = 1000, std::ostream* history = nullptr, StateObserver observer = {}, int snapshot_interval = 1, double final_time = std::numeric_limits<double>::infinity());
            void run_until(double final_time, T CFL = T{1}, int report_interval = 1000, std::ostream* history = nullptr, StateObserver observer = {}, int snapshot_interval = 1, int max_iterations = 200000);
            const CellsData<T>& get_CellsData() const { return cells_handler.get_CellsData(); }
            const initialization_data<T>& get_InitializationData() const { return cells_handler.get_InitializationData(); }
            struct Diagnostics { T min_density{}, min_pressure{}; int invalid_cells{}; double residual_rms = std::numeric_limits<double>::quiet_NaN(); };
            Diagnostics eval_diagnostics(const CellsData<T>& cellsdata);
            double eval_residual_rms(const CellsData<T>& cellsdata);    
        private:
            StateInfo state_info;
            CellsHandler<Kokkos::DefaultExecutionSpace, T> cells_handler;
    };
}
#include "SolverFVM.tpp"
