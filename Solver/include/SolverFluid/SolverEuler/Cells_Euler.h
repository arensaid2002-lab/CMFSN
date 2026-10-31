#pragma once

#include <Kokkos_Core.hpp>
#include <MeshFormat.h>

using namespace std;
using namespace Mesh;

using View1D = Kokkos::View<float*, Kokkos::LayoutLeft>;
using View2D = Kokkos::View<float**, Kokkos::LayoutLeft>;
using View3D = Kokkos::View<float***, Kokkos::LayoutLeft>;

namespace Solver
{
    struct CellsData
    {
        CellsData(){}
        View1D Cell_Areas;
        View1D rho_old, Momentum_x_old, Momentum_y_old, E_old;
        View1D rho_new, Momentum_x_new, Momentum_y_new, E_new;
        View2D Flux_mass, Flux_Momentum_x, Flux_Momentum_y, Flux_E, D;
        View3D Normals;
        
    };

    template <typename ExecSpace = Kokkos::DefaultExecutionSpace>
    class Cells_Euler
    {
        public:


            Cells_Euler(const MeshData_Euler& Mesh_Euler);
            ~Cells_Euler() = default;

            init
            KOKKOS_FUNCTION
            static float eval_Cell_Area(int i, int j);

            KOKKOS_FUNCTION
            static Kokkos::pair<float, float> eval_Cell_Normal(int i, int j, int direction);
            KOKKOS_FUNCTION
            static int Node_index(int i, int j) {return i + MeshData.Ni*j;}
            KOKKOS_FUNCTION
            static int Cell_index(int i, int j) {return i + (MeshData.Ni-1)*j;}

            //CellsData get_CellsData() const {return CellsData_;}
        private:
            const MeshData_Euler& MeshData;
            CellsData CellsData_;
    };
}
#include "Cells_Euler.tpp"