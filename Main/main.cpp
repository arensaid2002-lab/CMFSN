#include <iostream>
#include <filesystem>
#include "MeshHandler.h"
#include "Errors.h"
#include "FileIO.h"

using namespace Mesh;
using namespace std;

int main(int argc, char* argv[])
{
    Kokkos::ScopeGuard guard(argc, argv);

    MeshHandler meshHandler;
    MeshData_Euler& MeshData = meshHandler.get_Mesh_Euler();
    
    string path = "../../../../NACA_0012_mesh_PLOT3D_format/x.1";
  
    IO::Import_Mesh(path, MeshData, MeshHandler::MeshFileFormat::Plot3D);
    
    for (int i = 0; i < MeshData.N_nodes; i++)
    {
        cout << MeshData.X(i)  << endl;
    }

    return 0;
}
