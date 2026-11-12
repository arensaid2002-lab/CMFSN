#include <iostream>
#include "RunSimulation.h"

int main(int argc, char* argv[])
{
    Kokkos::ScopeGuard guard(argc, argv);
    try
    {
        using T = double;
        constexpr auto SType = Solver::SolverTypes::RANS2D;
        constexpr auto TimeScheme = Solver::TimeEvolutionScheme::Steady;
        constexpr auto RKScheme = Solver::RungeKuttaScheme::Jameson;
        CMFSN::SimulationParameters<T> parameters;
        parameters.Mach_inf = T{2};
        parameters.gamma = T(1.4);
        parameters.alpha = T{0}; // Radians.
        parameters.wall_condition = Solver::BoundaryCondition::AdiabaticWall;
        parameters.wall_temperature = T{1}; // Nondimensional p/rho.
        parameters.CFL = T(1.5);
        parameters.final_time = 1.;
        parameters.iterations = 100000;
        parameters.Re_inf = T{5000 }; // Set a finite positive Reynolds number for laminar flow.
        parameters.Pr_inf = T{0.72};
        parameters.mesh_path = std::filesystem::path(__FILE__).parent_path().parent_path() / "NACA_0012_mesh_PLOT3D_format" / "x.4";
        parameters.output_root = "results";
        parameters.report_interval = 1000;
        parameters.snapshot_interval = 1;
        return CMFSN::RunSimulation<T, SType, TimeScheme, RKScheme>(parameters);
    }
    catch (const std::exception& error)
    {
        std::cerr << "Solver failed: " << error.what() << '\n';
        return 1;
    }
}
