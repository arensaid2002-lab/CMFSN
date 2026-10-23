#include <Kokkos_Core.hpp>
#include <Eigen/Dense>
#include <iostream>
#include <MeshFormat.h>

using namespace std;
using namespace Mesh;
namespace Solver
{

    class Cell_Euler
    {
        public:


            Cell_Euler(int i, int j, MeshData_Euler& MeshData);

            enum class direction : int
            {
                Down = 0, Right = 1, Up = 2, Left = 3
            };

            int& Mesh_Ni, Mesh_Nj; 

            float& rho();

            float& Momentum_x();

            float& Momentum_y();

            float& E();

            float& Flux_rho(direction d);

            float& Flux_Momentum_x(direction d);

            float& Flux_Momentum_y(direction d);

            float& Flux_E(direction d);

            int get_Position_index(int i, int j);

            void eval_Normals(int i, int j, MeshData_Euler MeshData);

            void init_State()
            {
                W << 1.0, 0.0, 0.0, 1.0;
                F_c.setZero();
                D.setZero();
            }
        private:
            Eigen::Vector4f W;      //Conservative Matrix quantities    
            Eigen::Matrix4f F_c, D; //Flux and Dissipation 
            Eigen::Matrix4f Normals;//Normals 




 
    };
}