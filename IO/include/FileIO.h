#pragma once

#include <iostream>
#include <MeshHandler.h>

using namespace std;
using namespace Mesh;
namespace IO
{
    template <typename T, typename MeshFileFormat Format = MeshFileFormat::Plot3D>
    void Import_Mesh(string path, MeshData_Structured<T>& MeshData);
}
#include "FileIO.tpp"
