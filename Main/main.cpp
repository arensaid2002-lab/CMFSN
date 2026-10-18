#include <iostream>
#include <filesystem>
#include "MeshHandler.h"
#include "Errors.h"
#include "FileIO.h"

using namespace Mesh;
using namespace std;

int main(int, char**)
{
    MeshHandler meshHandler;
    MeshData_Euler& MeshData = meshHandler.get_Mesh_Euler();
    
    std::cout<< std::filesystem::current_path() << '\n'; 

    string path = "../../../../NACA_0012_mesh_PLOT3D_format/x.9";
    try
    {
        IO::Import_Mesh(path, MeshData, MeshHandler::MeshFileFormat::Plot3D);
    }
    catch (const Errors::FileOpenError& Error)
    {
        cerr << Error.what() << endl;
    }

    for (const auto& valeur : MeshData.X)
    {
        cout << valeur<< endl;
    }

    std::cout << "Hello, from CMFSN!\n";
}
