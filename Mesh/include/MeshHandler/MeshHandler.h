#pragma once 

#include <iostream>

#include <MeshFormat.h>

using namespace std;
namespace Mesh
{
    enum class MeshFileFormat
    {
        Plot3D //Ajouter des types de fichier supp
    };
        template <typename T = float>
        class MeshHandler
        {
            public:
        

        
        
            public:
        
        
            MeshHandler() = default;
            ~MeshHandler() = default;
        
            MeshData_Structured<T>& get_Mesh_Euler();

            private:
            MeshData_Structured<T> Mesh_Euler;

        };
}
#include "MeshHandler.tpp"
