#pragma once

#include <iostream>
#include <MeshHandler.h>

using namespace std;
using namespace Mesh;
namespace IO
{
    template <typename T>
    void Import_Mesh(string path, MeshData_Structured<T>& Mesh_Euler, typename MeshHandler<T>::MeshFileFormat Format);
}
#include "FileIO.tpp"
