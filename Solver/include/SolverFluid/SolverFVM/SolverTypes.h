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
}
