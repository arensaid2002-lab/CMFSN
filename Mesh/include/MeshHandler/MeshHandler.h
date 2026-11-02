#pragma once 

#include <iostream>

#include <MeshFormat.h>

using namespace std;
namespace Mesh
{
    template <typename T = float>
    class MeshHandler
    {
        public:
        
        enum class MeshFileFormat
        {
            Plot3D //Ajouter des types de fichier supp
        };
        
        
        public:
        
        
        MeshHandler() = default;
        ~MeshHandler() = default;
        
        MeshData_Euler<T>& get_Mesh_Euler();

        private:
        MeshData_Euler<T> Mesh_Euler;

    };
}
#include "MeshHandler.tpp"
