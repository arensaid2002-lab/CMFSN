#include <iostream>
#include "RunSimulation.h"

int main(int argc, char* argv[])
{
    using T = double;
    constexpr auto SType = Solver::SolverTypes::Euler2D;
    constexpr auto TimeScheme = Solver::TimeEvolutionScheme::Steady;
    constexpr auto RKScheme = Solver::RungeKuttaScheme::Jameson;

    CMFSN::SimulationParameters<T> parameters;
    parameters.flow.Mach_inf = T{0.5};
    parameters.flow.gamma = T(1.4);
    parameters.flow.alpha = T{0}; // Radians.
    parameters.flow.Re_inf = T{5000}; // Ignored by Euler2D.
    parameters.flow.Pr_inf = T(0.72);
    parameters.flow.wall_condition = Solver::BoundaryCondition::AdiabaticWall;
    parameters.CFL = T(3);
    parameters.iterations = 100000;
    parameters.tol = T(1e-16); // Steady only: stops when the global and the density RMS residuals both reach it.
    parameters.final_time = 1.; // UnSteady only.
    parameters.output_root = "results";
    // Grid study (Vassberg & Jameson family): x.7 = 32x32 cells, x.6 = 64, x.5 = 128, x.4 = 256, x.3 = 512, x.2 = 1024, x.1 = 2048.
    // Each mesh is written to results/study_NNNN/<mesh file name>; study.csv lists the status of every mesh.
    const auto mesh_directory = std::filesystem::path(__FILE__).parent_path().parent_path() / "NACA_0012_mesh_PLOT3D_format";
    const std::vector<std::filesystem::path> meshes = { mesh_directory / "x.7", mesh_directory / "x.6", mesh_directory / "x.5", mesh_directory / "x.4" };
    parameters.report_interval = 1000;
    parameters.history_interval = 1; // residuals.csv rows (rms, rms_rho, ...).
    parameters.snapshot_interval = 1;

    Kokkos::ScopeGuard guard(argc, argv);
    try
    {
        const int failures = CMFSN::RunStudy<T, SType, TimeScheme, RKScheme>(parameters, meshes);
        if (failures != 0)
            std::cerr << failures << " mesh(es) failed, see study.csv\n";
        return failures == 0 ? 0 : 1;
    }
    catch (const std::exception& error)
    {
        std::cerr << "Solver failed: " << error.what() << '\n';
        return 1;
    }
}
