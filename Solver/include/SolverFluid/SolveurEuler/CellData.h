#include <Kokkos_Core.hpp>
namespace Solver
{
    struct Cell_Euler
    {
        int i{}, j{};
        
        Kokkos::Array<float, 4> W;
         

    };
}