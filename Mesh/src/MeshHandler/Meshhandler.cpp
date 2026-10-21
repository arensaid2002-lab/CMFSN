#include <iostream>

#include <MeshHandler.h>

namespace Mesh
{
    MeshData_Euler& MeshHandler::get_Mesh_Euler()
    {
        return Mesh_Euler;
    }
}