#include <iostream>
#include "RunSimulation.h"

int main(int argc, char* argv[])
{
    using T = double;
    constexpr auto SType = Solver::SolverTypes::RANS2D;
    constexpr auto TimeScheme = Solver::TimeEvolutionScheme::Steady;
    constexpr auto RKScheme = Solver::RungeKuttaScheme::Jameson;

    CMFSN::SimulationParameters<T> parameters;
    parameters.flow.Mach_inf = T{0.5};
    parameters.flow.gamma = T(1.4);
    parameters.flow.alpha = T{0}; // Radians.
    parameters.flow.Re_inf = T{5000}; // Ignored by Euler2D.
    parameters.flow.Pr_inf = T(0.72);
    parameters.flow.wall_condition = Solver::BoundaryCondition::AdiabaticWall;
    parameters.flow.wall_temperature = T{1}; // Nondimensional p/rho, isothermal wall only.
    parameters.CFL = T(3.9);
    parameters.iterations = 100000;
    parameters.final_time = 1.; // UnSteady only.
    parameters.mesh_path = std::filesystem::path(__FILE__).parent_path().parent_path() / "NACA_0012_mesh_PLOT3D_format" / "x.4";
    parameters.output_root = "results";
    parameters.report_interval = 1000;
    parameters.snapshot_interval = 1;

    Kokkos::ScopeGuard guard(argc, argv);
    try
    {
        CMFSN::RunSimulation<T, SType, TimeScheme, RKScheme>(parameters);
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "Solver failed: " << error.what() << '\n';
        return 1;
    }
}
