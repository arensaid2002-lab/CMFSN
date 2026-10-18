#include <iostream>
#include <fstream>
#include "Errors.h"
#include "MeshHandler.h"

using namespace std;
using namespace Mesh;

namespace IO
{
    void Import_Mesh(string path, MeshData_Euler& Mesh_Euler, MeshHandler::MeshFileFormat Format)
    {

        fstream Meshfile(path);
        if (Meshfile.is_open())
        {
            float position_value;
            Meshfile >> Mesh_Euler.Ni >> Mesh_Euler.Nj;
            Mesh_Euler.N = Mesh_Euler.Ni * Mesh_Euler.Nj;

            for(int i = 0; i < Mesh_Euler.N && Meshfile >> position_value; i++)
            {
                Mesh_Euler.X.push_back(position_value);
            }

            for(int i = 0; i < Mesh_Euler.N && Meshfile >> position_value; i++)
            {
                Mesh_Euler.Y.push_back(position_value);
            }
        }
        else
        {
            throw Errors::FileOpenError("Error: Mesh file could not be opened");
        }
        Meshfile.close();
    }
    
}