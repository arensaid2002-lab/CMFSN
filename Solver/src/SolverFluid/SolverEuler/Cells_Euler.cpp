#include <Kokkos_Core.hpp>
#include <Eigen/Dense>
#include <iostream>
#include <MeshFormat.h>
#include <SolverEuler/Cells_Euler.h>

using namespace std;
using namespace Mesh;
namespace Solver
{

        Cells_Euler::Cells_Euler(const MeshData_Euler& Mesh_Euler): MeshData(Mesh_Euler)
        {
                Normals = Kokkos::View<float***, Kokkos::LayoutLeft >("Normals", MeshData.N_Cells, 4, 2);
                Cell_Areas = Kokkos::View<float*, Kokkos::LayoutLeft >("Cell_Areas", MeshData.N_Cells);

                rho_old = Kokkos::View<float*, Kokkos::LayoutLeft >("rho_old", MeshData.N_Cells);
                Momentum_x_old = Kokkos::View<float*, Kokkos::LayoutLeft >("Momentum_x_old", MeshData.N_Cells);
                Momentum_y_old = Kokkos::View<float*, Kokkos::LayoutLeft >("Momentum_y_old", MeshData.N_Cells);
                E_old = Kokkos::View<float*, Kokkos::LayoutLeft >("E_old", MeshData.N_Cells);

                rho_new = Kokkos::View<float*, Kokkos::LayoutLeft >("rho_new", MeshData.N_Cells);
                Momentum_x_new = Kokkos::View<float*, Kokkos::LayoutLeft >("Momentum_x_new", MeshData.N_Cells);
                Momentum_y_new = Kokkos::View<float*, Kokkos::LayoutLeft >("Momentum_y_new", MeshData.N_Cells);
                E_new = Kokkos::View<float*, Kokkos::LayoutLeft >("E_new", MeshData.N_Cells);

                Flux_mass = Kokkos::View<float**, Kokkos::LayoutLeft >("Flux_mass", MeshData.N_Cells, 4);
                Flux_Momentum_x = Kokkos::View<float**, Kokkos::LayoutLeft >("Flux_Momentum_x", MeshData.N_Cells, 4);
                Flux_Momentum_y = Kokkos::View<float**, Kokkos::LayoutLeft >("Flux_Momentum_y", MeshData.N_Cells, 4);
                Flux_E = Kokkos::View<float**, Kokkos::LayoutLeft >("Flux_E", MeshData.N_Cells, 4);
                D = Kokkos::View<float**, Kokkos::LayoutLeft >("D", MeshData.N_Cells, 4);
                
                Kokkos::parallel_for("eval_Cell_Areas_And_Normals", Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, MeshData.N_Cells), KOKKOS_LAMBDA(const int Cell_k)
                {
                    int i = Cell_k % (MeshData.Ni-1);
                    int j = Cell_k / (MeshData.Ni-1);
                    Cell_Areas(Cell_k) = eval_Cell_Area(i, j);
                    for(int direction = 0; direction < 4; direction++)
                    {
                        Normals(Cell_k, direction, 0) = Cells_Euler::eval_Cell_Normal(i, j, direction).first;
                        Normals(Cell_k, direction, 1) = Cells_Euler::eval_Cell_Normal(i, j, direction).second;
                    }
                });
                init_State();
        }



        KOKKOS_INLINE_FUNCTION
        std::pair<float, float> Cells_Euler::eval_Cell_Normal(int i, int j, int direction)
        {
            float delta_x, delta_y;
            int CellIndex = Cell_index(i, j);
            int NodeIndex = Node_index(i, j);
                        
            switch (direction)
            {
                case 0: // Bottom
                    delta_x = MeshData.X(Node_index(i + 1,j)) - MeshData.X(Node_index(i,j));
                    delta_y = MeshData.Y(Node_index(i + 1,j)) - MeshData.Y(Node_index(i,j));

                    break;
                
                case 1: // Right
                    delta_x = MeshData.X(Node_index(i+1,j + 1)) - MeshData.X(Node_index(i+1,j));
                    delta_y = MeshData.Y(Node_index(i+1,j + 1)) - MeshData.Y(Node_index(i+1,j));
                    
                    break;
                case 2: // Top
                    delta_x = MeshData.X(Node_index(i,j + 1)) - MeshData.X(Node_index(i + 1,j + 1));
                    delta_y = MeshData.Y(Node_index(i,j + 1)) - MeshData.Y(Node_index(i + 1,j + 1));
                    
                    break;
                case 3: // Left
                    delta_x = MeshData.X(Node_index(i,j)) - MeshData.X(Node_index(i,j + 1));
                    delta_y = MeshData.Y(Node_index(i,j)) - MeshData.Y(Node_index(i,j + 1));
                    
                    break;
            }
            return {delta_y, - delta_x};
        }

}