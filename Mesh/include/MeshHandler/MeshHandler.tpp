#pragma once

#include <iostream>

#include <MeshHandler.h>

namespace Mesh
{
    template <typename T>
    MeshData_Euler<T>& MeshHandler<T>::get_Mesh_Euler()
    {
        return Mesh_Euler;
    }
}