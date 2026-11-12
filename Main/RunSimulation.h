#pragma once

#include <filesystem>
#include <fstream>
#include <limits>
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
        T Mach_inf = T{2}, gamma = T(1.4), alpha = T{0}; // alpha is in radians.
        int N_ghost_layers = 3;
        Solver::BoundaryCondition wall_condition = Solver::BoundaryCondition::AdiabaticWall;
        T wall_temperature = T{1}; // Nondimensional temperature, consistent with p/rho.
        T CFL = T(3.9);
        double final_time = 1.;
        int iterations = 100000;
        T Re_inf = std::numeric_limits<T>::infinity(), Pr_inf = T{1};
        std::filesystem::path mesh_path = std::filesystem::path(__FILE__).parent_path().parent_path() / "NACA_0012_mesh_PLOT3D_format" / "x.1";
        std::filesystem::path output_root = "results";
        int report_interval = 1000, snapshot_interval = 1;
    };

    template<typename T, Solver::SolverTypes SType, Solver::TimeEvolutionScheme TimeScheme, Solver::RungeKuttaScheme RKScheme>
    int RunSimulation(const SimulationParameters<T>& parameters)
    {
        Mesh::MeshHandler<T> mesh_handler;
        auto& meshdata = mesh_handler.get_Mesh_Euler();
        IO::Import_Mesh<T, Mesh::MeshFileFormat::Plot3D>(parameters.mesh_path.string(), meshdata);
        Solver::SolverFVM<T, SType, TimeScheme, RKScheme> solver(parameters.Mach_inf, parameters.gamma, parameters.alpha, parameters.Re_inf, parameters.Pr_inf, parameters.N_ghost_layers, meshdata, parameters.wall_condition, parameters.wall_temperature);
        const auto output_dir = IO::Create_RunDirectory(parameters.output_root);
        std::ofstream history(output_dir / "residuals.csv");
        if (!history)
            throw std::runtime_error("Could not open residual history file");
        history.imbue(std::locale::classic());
        if constexpr (TimeScheme == Solver::TimeEvolutionScheme::UnSteady)
        {
            IO::SolutionSeriesWriter series(output_dir, meshdata, solver.get_InitializationData());
            solver.run_until(parameters.final_time, parameters.CFL, parameters.report_interval, &history, [&](const auto& info, const auto& cellsdata) { series.write(cellsdata, info.iteration, info.physical_time, info.time_step); }, parameters.snapshot_interval, parameters.iterations);
        }
        else
        {
            solver.run(parameters.iterations, parameters.CFL, parameters.report_interval, &history);
            IO::Export_Solution<T>((output_dir / "solution.vts").string(), meshdata, solver.get_CellsData(), solver.get_InitializationData(), solver.get_StateInfo().iteration);
        }
        return 0;
    }
}
