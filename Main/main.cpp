#include <iostream>
#include <filesystem>
#include "MeshHandler.h"
#include "FileIO.h"
#include "SolverFluid/SolverFVM/CellsHandler/CellsHandler.h"

using namespace Mesh;
using namespace std;

int main(int argc, char* argv[])
{
    Kokkos::ScopeGuard guard(argc, argv);
    const auto start = chrono::steady_clock::now();
    using T = float;
    MeshHandler<T> meshHandler;
    MeshData_Structured<T>& MeshData = meshHandler.get_Mesh_Euler();

    string path = "C:/Users/wali2/source/repos/arensaid2002-lab/CMFSN/NACA_0012_mesh_PLOT3D_format/x.1";

    IO::Import_Mesh<T, MeshFileFormat::Plot3D>(path, MeshData);
    
    Solver::CellsHandler<Kokkos::DefaultExecutionSpace, T> cells_handler(1,1,2,1 ,MeshData);
	auto cellsData = cells_handler.get_CellsData();
    const auto Normals = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, cellsData.Normals);
    for (int direction = 0; direction < 4; direction++)
    {
        cout << Normals(4104, direction, 0) << " " << Normals(4104, direction, 1) << endl;
    }
    
    const auto end = chrono::steady_clock::now();
    const auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);
    cout << "Execution time: " << duration.count() << " ms" << endl;

    return 0;
}
