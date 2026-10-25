#include <Kokkos_Core.hpp>
#include <Eigen/Dense>
#include <iostream>
#include <MeshFormat.h>
#include <SolverEuler/Cells_Euler.h>

using namespace std;
using namespace Mesh;
namespace Solver
{

        Cells_Euler::Cells_Euler(int i, int j, MeshData_Euler& MeshData)
        {
                Ni = MeshData.Ni;
                Nj = MeshData.Nj;
                N_nodes = MeshData.N_nodes;
                N_Cells = (Ni - 1) * (Nj - 1);

                eval_Normals(i, j, MeshData);
                init_State();
        }



        KOKKOS_INLINE_FUNCTION
        void Cells_Euler::eval_Cell_normals(int i, int j)
        {
            int Index = get_Cell_index(int i, int j)
            {

            }
            

            for(int direction = 0; direction < 4; direction++ )
            {
                switch (direction)
                {
                case 0:
                    Normals.col(direction)[0] = MeshData.X[Position_index + 1] - 
                    break;
                
                default:
                    break;
                }
                Normals.col()
            }

        }


        KOKKOS_INLINE_FUNCTION
        float& Cell_Euler::rho() {return W[0];}

        KOKKOS_INLINE_FUNCTION
        float& Momentum_x() {return W[1];}

        KOKKOS_INLINE_FUNCTION
        float& Momentum_y() {return W[2];}

        KOKKOS_INLINE_FUNCTION
        float& E() {return W[3];}

        KOKKOS_INLINE_FUNCTION
        auto& Flux_rho(direction d) {return F_c.col(static_cast<Eigen::Index>(d))[0];}

        KOKKOS_INLINE_FUNCTION
        auto& Flux_Momentum_x(direction d) {return F_c.col(static_cast<Eigen::Index>(d))[1];}

        KOKKOS_INLINE_FUNCTION
        auto& Flux_Momentum_y(direction d) {return F_c.col(static_cast<Eigen::Index>(d))[2];}

        KOKKOS_INLINE_FUNCTION
        auto& Flux_E(direction d) {return F_c.col(static_cast<Eigen::Index>(d))[3];}

        KOKKOS_INLINE_FUNCTION
        int get_Position_index(int i, int j)
        {
            return Position_index = i + j*MeshData.Ni;
        }


}