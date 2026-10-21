#pragma once 

#include <iostream>

#include <MeshFormat.h>

using namespace std;
namespace Mesh
{
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
        
        MeshData_Euler& get_Mesh_Euler();

        private:

        MeshData_Euler Mesh_Euler;

    };
}