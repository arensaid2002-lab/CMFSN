#pragma once

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <locale>
#include <sstream>
#include <stdexcept>
#include <vector>
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
        T tol = T(1e-6);          // Steady only: stops when the density RMS residual reaches it.
        std::filesystem::path mesh_path;
        std::filesystem::path output_root = "results";
        std::filesystem::path output_directory; // Empty: a new output_root/run_NNNN. Otherwise written here (created if needed).
        int report_interval = 1000, snapshot_interval = 1;
        int history_interval = 1; // Rows of residuals.csv (console lines follow report_interval).
    };

    template<typename T, Solver::SolverTypes SType, Solver::TimeEvolutionScheme TimeScheme, Solver::RungeKuttaScheme RKScheme>
    void RunSimulation(const SimulationParameters<T>& parameters)
    {
        Mesh::MeshHandler<T> mesh_handler;
        auto& meshdata = mesh_handler.get_Mesh_Euler();
        IO::Import_Mesh<T, Mesh::MeshFileFormat::Plot3D>(parameters.mesh_path.string(), meshdata);
        Solver::SolverFVM<T, SType, TimeScheme, RKScheme> solver(parameters.flow, meshdata);

        std::filesystem::path output_dir = parameters.output_directory;
        if (output_dir.empty())
            output_dir = IO::Create_RunDirectory(parameters.output_root);
        else
            std::filesystem::create_directories(output_dir);
        std::ofstream history(output_dir / "residuals.csv");
        if (!history)
            throw std::runtime_error("Could not open residual history file");
        history.imbue(std::locale::classic());

        if constexpr (TimeScheme == Solver::TimeEvolutionScheme::UnSteady)
        {
            IO::SolutionSeriesWriter series(output_dir, meshdata, solver.get_InitializationData());
            const auto write_snapshot = [&](const auto& info, const auto& cellsdata) { series.write(cellsdata, info.iteration, info.physical_time, info.time_step); };
            solver.run_until(parameters.final_time, parameters.CFL, parameters.report_interval, parameters.history_interval, &history, write_snapshot, parameters.snapshot_interval, parameters.iterations);
        }
        else
        {
            solver.run(parameters.iterations, parameters.tol, parameters.CFL, parameters.report_interval, parameters.history_interval, &history);
            IO::Export_Solution<T>((output_dir / "solution.vts").string(), meshdata, solver.get_CellsData(), solver.get_InitializationData(), solver.get_StateInfo().iteration);
        }
    }

    // Grid study: the same conditions on every mesh, each written to output_root/study_NNNN/<mesh file name>.
    // A failing mesh (divergence, unreadable file) is reported and skipped; study.csv lists the status of every mesh.
    // Returns the number of failed meshes.
    template<typename T, Solver::SolverTypes SType, Solver::TimeEvolutionScheme TimeScheme, Solver::RungeKuttaScheme RKScheme>
    int RunStudy(SimulationParameters<T> parameters, const std::vector<std::filesystem::path>& meshes)
    {
        if (meshes.empty())
            throw std::invalid_argument("A grid study needs at least one mesh");
        for (std::size_t m = 0; m < meshes.size(); ++m)
            for (std::size_t n = m + 1; n < meshes.size(); ++n)
                if (meshes[m].filename() == meshes[n].filename())
                    throw std::invalid_argument("Two study meshes share the file name " + meshes[m].filename().string());
        const auto study_dir = IO::Create_RunDirectory(parameters.output_root, "study_");
        std::ofstream summary(study_dir / "study.csv");
        if (!summary)
            throw std::runtime_error("Could not open the study summary file");
        summary.imbue(std::locale::classic());
        summary << "mesh,directory,status,message\n";
        int failures = 0;
        for (const auto& mesh : meshes)
        {
            parameters.mesh_path = mesh;
            parameters.output_directory = study_dir / mesh.filename();
            std::cout << "=== " << mesh.filename().string() << " -> " << parameters.output_directory.string() << " ===\n" << std::flush;
            std::string status = "ok", message;
            try
            {
                RunSimulation<T, SType, TimeScheme, RKScheme>(parameters);
            }
            catch (const std::exception& error)
            {
                ++failures;
                status = "failed";
                message = error.what();
                std::cerr << "Mesh " << mesh.filename().string() << " failed: " << message << '\n';
            }
            const auto csv_field = [](std::string text) // Quoted, with doubled quotes: paths and messages may hold commas.
            {
                for (std::size_t at = text.find('"'); at != std::string::npos; at = text.find('"', at + 2))
                    text.insert(at, 1, '"');
                return '"' + text + '"';
            };
            summary << csv_field(mesh.string()) << ',' << csv_field(mesh.filename().string()) << ',' << status << ',' << csv_field(message) << '\n' << std::flush;
        }
        return failures;
    }
}
