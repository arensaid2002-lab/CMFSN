#include <Kokkos_Core.hpp>
#include "MeshFormat.h"
#include "SolverFluid/SolverEuler/Cells_Euler.h"

using namespace std;
using namespace Mesh;

using View1D = Kokkos::View<float*, Kokkos::LayoutLeft>;
using View2D = Kokkos::View<float**, Kokkos::LayoutLeft>;
using View3D = Kokkos::View<float***, Kokkos::LayoutLeft>;

namespace Solver
{
        template <typename ExecSpace>
        Cells_Euler<ExecSpace>::Cells_Euler(const MeshData_Euler& Mesh_Euler) : MeshData(Mesh_Euler)
        {
                Normals = View3D("Normals", MeshData.N_Cells, MeshData.N_edges, 2);
                Cell_Areas = View1D("Cell_Areas", MeshData.N_Cells);

                rho_old = View1D("rho_old", MeshData.N_Cells);
                Momentum_x_old = View1D("Momentum_x_old", MeshData.N_Cells);
                Momentum_y_old = View1D("Momentum_y_old", MeshData.N_Cells);
                E_old = View1D("E_old", MeshData.N_Cells);

                rho_new = View1D("rho_new", MeshData.N_Cells);
                Momentum_x_new = View1D("Momentum_x_new", MeshData.N_Cells);
                Momentum_y_new = View1D("Momentum_y_new", MeshData.N_Cells);
                E_new = View1D("E_new", MeshData.N_Cells);

                Flux_mass = View2D("Flux_mass", MeshData.N_Cells, 4);
                Flux_Momentum_x = View2D("Flux_Momentum_x", MeshData.N_Cells, 4);
                Flux_Momentum_y = View2D("Flux_Momentum_y", MeshData.N_Cells, 4);
                Flux_E = View2D("Flux_E", MeshData.N_Cells, 4);
                D = View2D("D", MeshData.N_Cells, 4);

        }



        KOKKOS_FUNCTION
        static Kokkos::pair<float, float> Cells_Euler::eval_Cell_Normal(int i, int j, int direction)
        {
            float delta_x, delta_y;
            int CellIndex = Cells_Euler::Cell_index(i, j);
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