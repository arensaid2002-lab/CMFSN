#pragma once

namespace Solver
{
    enum class SolverTypes
    {
        Euler2D,
        RANS2D,
    };

    enum class TimeEvolutionScheme
    {
        Steady,
        UnSteady,
    };

    enum class RungeKuttaScheme
    {
        Unspecified,
        Jameson,
        CarpenterKennedy,
    };
    enum class BoundaryCondition
    {
		AdiabaticWall,
		isothermalWall,
    };
}
