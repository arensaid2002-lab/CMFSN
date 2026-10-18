
#include <iostream>
#include <Vector>

using namespace std;
namespace Mesh
{
    struct MeshData_Euler
    {
        int Ni = 0;
        int Nj = 0;
        int N;

        vector<float> X;
        vector<float> Y;
    };
}