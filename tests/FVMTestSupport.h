#pragma once
#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>
#include <string>
#include "SolverFluid/SolverFVM/SolverFVM.h"

namespace Tests
{
    using E=Kokkos::DefaultExecutionSpace;
    using namespace Solver;
    using namespace Solver::CellOps;
    inline void require(bool value,const std::string& message)
    {
        if(!value) throw std::runtime_error(message);
    }
    template<typename T> double tolerance()
    {
        return std::is_same_v<T,float>?3e-4:3e-11;
    }
    inline void near(double actual,double expected,double tol,const std::string& message)
    {
        if(!std::isfinite(actual)||std::abs(actual-expected)>tol*std::max(1.,std::abs(expected)))
            throw std::runtime_error(message+": actual="+std::to_string(actual)+", expected="+std::to_string(expected));
    }
    template<typename F> void invalid(F&& function,const std::string& message)
    {
        bool caught=false;
        try {function();} catch(const std::invalid_argument&) {caught=true;}
        require(caught,message);
    }
    template<typename T> Mesh::MeshData_Structured<T> make_mesh(bool curved=false,int ni=10,int nj=8)
    {
        Mesh::MeshData_Structured<T> mesh;
        mesh.Ni=ni+1;mesh.Nj=nj+1;mesh.N_nodes=(ni+1)*(nj+1);mesh.N_Cells=ni*nj;
        Kokkos::View<T*> x("test_x",mesh.N_nodes),y("test_y",mesh.N_nodes);
        auto hx=Kokkos::create_mirror_view(x),hy=Kokkos::create_mirror_view(y);
        for(int j=0;j<=nj;++j) for(int i=0;i<=ni;++i)
        {
            const int k=i+(ni+1)*j;
            hx(k)=T(.07)*i+T(.01)*j+(curved?T(.0002)*i*j:T{0});
            hy(k)=T(.02)*i+T(.04)*j+(curved?T(.00015)*i*i:T{0});
        }
        Kokkos::deep_copy(x,hx);Kokkos::deep_copy(y,hy);mesh.X=x;mesh.Y=y;
        return mesh;
    }
    template<typename T> struct Fixture
    {
        Mesh::MeshData_Structured<T> mesh;
        CellsData<T> cells;
        initialization_data<T> init;
        Fixture(bool curved=false):mesh(make_mesh<T>(curved)),cells(3,mesh),init(T(.2),T(1.4),T{0},T{250},T(.72))
        {
            eval_GeometricalData<E>(cells,mesh);E().fence();
            uniform(T{1},T(1.2),T(.2),T(.3));
        }
        void uniform(T rho,T theta,T u,T v)
        {
            const T pressure=rho*theta;
            Kokkos::deep_copy(cells.rho_old,rho);
            Kokkos::deep_copy(cells.Momentum_u_old,rho*u);
            Kokkos::deep_copy(cells.Momentum_v_old,rho*v);
            Kokkos::deep_copy(cells.Pressure,pressure);
            Kokkos::deep_copy(cells.SoundSpeed,T(std::sqrt(init.gamma*theta)));
            Kokkos::deep_copy(cells.E_old,pressure/(init.gamma-T{1})+T(.5)*rho*(u*u+v*v));
        }
        int index(int i,int j) const {return Cell_index(i,j,cells.Ni_total_cells);}
    };
}
