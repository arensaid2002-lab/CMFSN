#pragma once

#include <filesystem>
#include <fstream>
#include <locale>
#include <stdexcept>
#include "MeshHandler.h"
#include "FileIO.h"
#include "SolverFluid/SolverFVM/SolverFVM.h"

namespace CMFSN
{
    template<typename T>
    struct SimulationParameters
    {
        Solver::FlowConditions<T> flow; // Mach, gamma, incidence, Reynolds, Prandtl and wall condition.
        T CFL = T{1};
        int iterations = 100000;  // Steady: number of iterations. UnSteady: maximum number of steps.
        double final_time = 1.;   // UnSteady only.
        std::filesystem::path mesh_path;
        std::filesystem::path output_root = "results";
        int report_interval = 1000, snapshot_interval = 1;
    };

    template<typename T, Solver::SolverTypes SType, Solver::TimeEvolutionScheme TimeScheme, Solver::RungeKuttaScheme RKScheme>
    void RunSimulation(const SimulationParameters<T>& parameters)
    {
        Mesh::MeshHandler<T> mesh_handler;
        auto& meshdata = mesh_handler.get_Mesh_Euler();
        IO::Import_Mesh<T, Mesh::MeshFileFormat::Plot3D>(parameters.mesh_path.string(), meshdata);
        Solver::SolverFVM<T, SType, TimeScheme, RKScheme> solver(parameters.flow, meshdata);

        const auto output_dir = IO::Create_RunDirectory(parameters.output_root);
        std::ofstream history(output_dir / "residuals.csv");
        if (!history)
            throw std::runtime_error("Could not open residual history file");
        history.imbue(std::locale::classic());

        if constexpr (TimeScheme == Solver::TimeEvolutionScheme::UnSteady)
        {
            IO::SolutionSeriesWriter series(output_dir, meshdata, solver.get_InitializationData());
            const auto write_snapshot = [&](const auto& info, const auto& cellsdata) { series.write(cellsdata, info.iteration, info.physical_time, info.time_step); };
            solver.run_until(parameters.final_time, parameters.CFL, parameters.report_interval, &history, write_snapshot, parameters.snapshot_interval, parameters.iterations);
        }
        else
        {
            solver.run(parameters.iterations, parameters.CFL, parameters.report_interval, &history);
            IO::Export_Solution<T>((output_dir / "solution.vts").string(), meshdata, solver.get_CellsData(), solver.get_InitializationData(), solver.get_StateInfo().iteration);
        }
    }
}
