#include <filesystem>
#include <fstream>
#include <iostream>
#include <locale>
#include <limits>
#include <stdexcept>
#include "MeshHandler.h"
#include "FileIO.h"
#include "SolverFluid/SolverFVM/SolverFVM.h"

int main(int argc, char* argv[])
{
    Kokkos::ScopeGuard guard(argc, argv);
    try
    {
        using T = double;
        constexpr auto SType = Solver::SolverTypes::Euler2D;
        constexpr auto TimeScheme = Solver::TimeEvolutionScheme::UnSteady;
        constexpr auto RKScheme = Solver::RungeKuttaScheme::CarpenterKennedy;
        const T Mach_inf = T(2), gamma = T(1.4), alpha = T{0}; // alpha is in radians.
        const int N_ghostlayers = 3;
        const double final_time = argc > 2 ? std::stod(argv[2]) : 1;
        const int iterations = TimeScheme == Solver::TimeEvolutionScheme::Steady && argc > 2 ? std::stoi(argv[2]) : 100000;
        const T CFL = argc > 3 ? T(std::stod(argv[3])) : T(1);
        const auto default_path = std::filesystem::path(__FILE__).parent_path().parent_path() / "NACA_0012_mesh_PLOT3D_format" / "x.3";
        const std::string path = argc > 1 ? argv[1] : default_path.string();
        const std::filesystem::path output_root = argc > 4 ? argv[4] : "results";
        const int report_interval = argc > 5 ? std::stoi(argv[5]) : 1000;

        const int snapshot_interval = argc > 6 ? std::stoi(argv[6]) : 1;
        // RANS2D takes Re and Pr after snapshot_interval on the command line.
        if constexpr (SType == Solver::SolverTypes::RANS2D)
            if (argc <= 8)
                throw std::invalid_argument("RANS2D requires Reynolds and Prandtl after snapshot_interval");
        const T Re_inf = argc > 7 ? T(std::stod(argv[7])) : std::numeric_limits<T>::infinity();
        const T Pr_inf = argc > 8 ? T(std::stod(argv[8])) : T{1};

        Mesh::MeshHandler<T> mesh_handler;
        auto& meshdata = mesh_handler.get_Mesh_Euler();
        IO::Import_Mesh<T, Mesh::MeshFileFormat::Plot3D>(path, meshdata);
        Solver::SolverFVM<T, SType, TimeScheme, RKScheme> solver(Mach_inf, gamma, alpha, Re_inf, Pr_inf, N_ghostlayers, meshdata);

        const auto output_dir = IO::Create_RunDirectory(output_root);
        std::ofstream history(output_dir / "residuals.csv");
        if (!history)
            throw std::runtime_error("Could not open residual history file");
        history.imbue(std::locale::classic());
        if constexpr (TimeScheme == Solver::TimeEvolutionScheme::UnSteady)
        {
            IO::SolutionSeriesWriter series(output_dir, meshdata, solver.get_InitializationData());
            solver.run_until(final_time, CFL, report_interval, &history, [&](const auto& info, const auto& cellsdata) { series.write(cellsdata, info.iteration, info.physical_time, info.time_step); }, snapshot_interval, iterations);
        }
        else
        {
            solver.run(iterations, CFL, report_interval, &history);
            IO::Export_Solution<T>((output_dir / "solution.vts").string(), meshdata, solver.get_CellsData(), solver.get_InitializationData(), solver.get_StateInfo().iteration);
        }
        history.close();
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "Solver test failed: " << error.what() << '\n';
        return 1;
    }
}
