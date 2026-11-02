#include "Cells_Euler.h"


namespace Solver
{
    class SolverEuler 
    {
        public:
			SolverEuler(const Mesh::MeshData_Euler& Mesh_Euler);

			void initialize_SolverEuler()
			{
				// Initialize any additional data or parameters specific to the Euler solver here
			}

        private:
            Cells_Euler<> cells_euler;

    };
}