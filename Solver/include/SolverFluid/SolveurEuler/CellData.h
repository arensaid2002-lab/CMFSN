#include <Kokkos_Core.hpp>
namespace Solver
{
    struct Cell_Euler
    {
        
        int i{}, j{};

        Kokkos::Array<float, 4> W; //{}
        Kokkos::Array<Kokkos::Array<float, 4>, 4> F_c;
        Kokkos::Array<Kokkos::Array<float, 2>, 4> Normals;

        KOKKOS_INLINE_FUNCTION
        float& rho() {return W[0];}

        KOKKOS_INLINE_FUNCTION
        float& Momentum_x() {return W[1];}

        KOKKOS_INLINE_FUNCTION
        float& Momentum_y() {return W[2];}

        KOKKOS_INLINE_FUNCTION
        float& E() {return W[3];}

    };
}