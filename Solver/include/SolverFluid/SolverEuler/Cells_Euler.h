#include <Kokkos_Core.hpp>
#include <Eigen/Dense>
#include <iostream>
#include <MeshFormat.h>

using namespace std;
using namespace Mesh;

namespace Solver
{

    class Cells_Euler
    {
        public:


            Cells_Euler(const MeshData_Euler& Mesh_Euler);

        private:

            float eval_Cell_Area(int i, int j);

            std::pair<float, float> eval_Cell_Normal(int i, int j, int direction);

            void init_State();

        public:
            const MeshData_Euler& MeshData;

            int Node_index(int i, int j) {return i + MeshData.Ni*j;}
            int Cell_index(int i, int j) {return i + (MeshData.Ni-1)*j;}
            float& rho(int Cell_k) {return rho_old(Cell_k);}
            float& Momentum_x(int Cell_k) {return Momentum_x_old(Cell_k);}
            float& Momentum_y(int Cell_k) {return Momentum_y_old(Cell_k);}
            float& E(int Cell_k) {return E_old(Cell_k);}

            float& F_mass(int cell_k, int direction) {return Flux_mass(cell_k, direction);}
            float& F_Momentum_x(int cell_k, int direction) {return Flux_Momentum_x(cell_k, direction);}
            float& F_Momentum_y(int cell_k, int direction) {return Flux_Momentum_y(cell_k, direction);}
            float& F_E(int cell_k, int direction) {return Flux_E(cell_k, direction);}
            Kokkos::View<float***, Kokkos::LayoutLeft> Normals;
        private:
            Kokkos::View<float*, Kokkos::LayoutLeft> Cell_Areas;
            Kokkos::View<float*, Kokkos::LayoutLeft> rho_old, Momentum_x_old, Momentum_y_old, E_old;
            Kokkos::View<float*, Kokkos::LayoutLeft> rho_new, Momentum_x_new, Momentum_y_new, E_new;
            
            
            Kokkos::View<float**, Kokkos::LayoutLeft> Flux_mass, Flux_Momentum_x, Flux_Momentum_y, Flux_E, D;
            





 
    };
}