#include <iostream>
#include "MeshHandler.h"
using namespace std;

namespace Mesh
{
    MeshData_Euler& MeshHandler::get_Mesh_Euler()
    {
        return Mesh_Euler;
    }
}