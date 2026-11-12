#include "FVMTestSupport.h"
#include "FileIO.h"
#include <array>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>

using namespace Tests;

template<typename T> void affine_fields(const CellsData<T>& cells)
{
    Kokkos::parallel_for("TestAffineFields",Kokkos::RangePolicy<E>(0,cells.N_total_cells),KOKKOS_LAMBDA(int k){
        const T x=cells.Cell_Centers(k,0),y=cells.Cell_Centers(k,1);
        const T u=T(.3)+T(.2)*x-T(.1)*y,v=T(-.2)+T(.15)*x+T(.25)*y;
        const T theta=T(1.5)+T(.05)*x+T(.1)*y;
        cells.rho_old(k)=T{1};cells.Momentum_u_old(k)=u;cells.Momentum_v_old(k)=v;
        cells.Pressure(k)=theta;cells.SoundSpeed(k)=Kokkos::sqrt(T(1.4)*theta);
        cells.E_old(k)=theta/T(.4)+T(.5)*(u*u+v*v);
    });
}
template<typename T> void physical_fluxes(const CellsData<T>& cells)
{
    Kokkos::parallel_for("TestPhysicalFluxes",Kokkos::RangePolicy<E>(0,cells.N_physical_cells),KOKKOS_LAMBDA(int p){
        const int k=Cell_index(p%cells.Ni_physical_cells+cells.N_ghost_layers,p/cells.Ni_physical_cells+cells.N_ghost_layers,cells.Ni_total_cells);
        for(int d=0;d<4;++d) Update_DataField<T,DataFields::Fluxes>(cells,k,d);
    });
}
template<typename T> void indexing()
{
    auto mesh=make_mesh<T>();CellsData<T> cells(mesh);
    require(cells.Ni_nodes==mesh.Ni && cells.Nj_nodes==mesh.Nj,"Node counts changed");
    require(cells.Ni_total_cells==mesh.Ni-1+6 && cells.Nj_total_cells==mesh.Nj-1+6,"Wrong total cell dimensions");
    require(cells.N_total_cells==cells.Ni_total_cells*cells.Nj_total_cells,"Wrong stored cell count");
    require(cells.DJacobien.extent(0)==cells.N_total_cells && cells.DJacobien.extent(1)==4 && cells.DJacobien.extent(2)==4,"Wrong Jacobian shape");
    for(int j=0;j<cells.Nj_total_cells;++j) for(int i=0;i<cells.Ni_total_cells;++i)
    {
        const int previous=i+((mesh.Ni+6)-1)*j;
        require(Cell_index(i,j,cells.Ni_total_cells)==previous,"Index changed during rename");
    }
    const int k=Cell_index(4,4,cells.Ni_total_cells);
    require(Cell_index_sides(k,cells.Ni_total_cells,0)==k-cells.Ni_total_cells,"Bottom neighbour");
    require(Cell_index_sides(k,cells.Ni_total_cells,1)==k+1,"Right neighbour");
    require(Cell_index_sides(k,cells.Ni_total_cells,2)==k+cells.Ni_total_cells,"Top neighbour");
    require(Cell_index_sides(k,cells.Ni_total_cells,3)==k-1,"Left neighbour");
}
template<typename T> void thermodynamics()
{
    const T gamma=T(1.4);
    for(T rho:{T(.8),T(1.3),T(3)}) for(T u:{T(-.3),T(.2)})
    {
        const T v=T(.15),p=T(1.7);
        const T energy=Solver::NumOps::eval_Energy(rho,u,v,p,gamma);
        near(Solver::NumOps::eval_pressure(rho,energy,u,v,gamma),p,tolerance<T>(),"Pressure/energy round trip");
        near(Solver::NumOps::eval_sound_speed(gamma,p,rho),std::sqrt(double(gamma*p/rho)),tolerance<T>(),"Sound speed");
    }
}
template<typename T> void geometry()
{
    for(bool curved:{false,true})
    {
        Fixture<T> f(curved);
        const auto x=Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),f.mesh.X);
        const auto y=Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),f.mesh.Y);
        const auto centers=Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),f.cells.Cell_Centers);
        const auto areas=Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),f.cells.Cell_Areas);
        const auto metrics=Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),f.cells.DJacobien);
        for(int j=0;j<f.cells.Nj_physical_cells;++j) for(int i=0;i<f.cells.Ni_physical_cells;++i)
        {
            const int nodes[4]={i+f.mesh.Ni*j,i+1+f.mesh.Ni*j,i+1+f.mesh.Ni*(j+1),i+f.mesh.Ni*(j+1)};
            long double twice_area=0,cx=0,cy=0;
            for(int a=0;a<4;++a)
            {
                const int b=(a+1)%4;
                const long double cross=(long double)x(nodes[a])*y(nodes[b])-(long double)x(nodes[b])*y(nodes[a]);
                twice_area+=cross;cx+=((long double)x(nodes[a])+x(nodes[b]))*cross;cy+=((long double)y(nodes[a])+y(nodes[b]))*cross;
            }
            const int k=f.index(i+3,j+3);
            near(areas(k),double(twice_area/2),tolerance<T>(),"Cell area");
            near(centers(k,0),double(cx/(3*twice_area)),tolerance<T>(),"Centroid x");
            near(centers(k,1),double(cy/(3*twice_area)),tolerance<T>(),"Centroid y");
            const int offset[4]={-f.cells.Ni_total_cells,1,f.cells.Ni_total_cells,-1};
            for(int d=0;d<4;++d)
            {
                const int other=k+offset[d],p=offset[(d+1)%4],m=offset[(d+3)%4];
                long double a[2],b[2];
                for(int c=0;c<2;++c) {a[c]=(long double)centers(other,c)-centers(k,c);b[c]=((long double)centers(k+p,c)-centers(k+m,c)+centers(other+p,c)-centers(other+m,c))/4;}
                const long double det=a[0]*b[1]-a[1]*b[0],ref[4]={a[0]/det,a[1]/det,b[0]/det,b[1]/det};
                for(int c=0;c<4;++c) near(metrics(k,d,c),double(ref[c]),tolerance<T>(),"Face metric");
            }
        }
    }
}
template<typename T> void gradients()
{
    for(bool curved:{false,true})
    {
        Fixture<T> f(curved);const auto cells=f.cells;affine_fields(cells);
        Kokkos::View<T***> result("test_gradients",cells.N_physical_cells,4,6);
        Kokkos::parallel_for("TestGradients",Kokkos::RangePolicy<E>(0,cells.N_physical_cells),KOKKOS_LAMBDA(int p){
            const int k=Cell_index(p%cells.Ni_physical_cells+3,p/cells.Ni_physical_cells+3,cells.Ni_total_cells);
            for(int d=0;d<4;++d) {
                const auto v=eval_gradient<T,DataFields::Velocity>(cells,k,d);
                const auto t=eval_gradient<T,DataFields::Temperature>(cells,k,d);
                for(int c=0;c<4;++c) result(p,d,c)=v[c];result(p,d,4)=t[0];result(p,d,5)=t[1];
            }
        });
        const auto host=Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),result);
        const double expected[6]={.2,.15,-.1,.25,.05,.1};
        for(int p=0;p<cells.N_physical_cells;++p) for(int d=0;d<4;++d) for(int c=0;c<6;++c)
            near(host(p,d,c),expected[c],tolerance<T>(),"Affine face gradient");
    }
}
template<typename T> void constitutive()
{
    const auto tau=Solver::NumOps::Tau(T{2},Kokkos::Array<T,4>{T(.2),T(.15),T(-.1),T(.25)});
    near(tau[0],.2,tolerance<T>(),"tau_xx");near(tau[1],.4,tolerance<T>(),"tau_yy");near(tau[2],.1,tolerance<T>(),"tau_xy");
}
template<typename T> void periodic()
{
    Fixture<T> f;const auto cells=f.cells;
    Kokkos::parallel_for("TestPeriodicSources",Kokkos::RangePolicy<E>(0,cells.N_total_cells),KOKKOS_LAMBDA(int k){cells.rho_old(k)=T(1)+T(.001)*k;});
    const auto before=Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),cells.rho_old);
    ApplyPeriodicity<E>(cells);
    const auto after=Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),cells.rho_old);
    for(int j=0;j<cells.Nj_total_cells;++j) for(int layer=0;layer<3;++layer)
    {
        const int left=f.index(layer,j),right=f.index(cells.Ni_physical_cells+3+layer,j);
        near(after(left),before(left+cells.Ni_physical_cells),0,"Periodic left including corners");
        near(after(right),before(f.index(3+layer,j)),0,"Periodic right including corners");
    }
    for(int j=3;j<3+cells.Nj_physical_cells;++j) for(int i=3;i<3+cells.Ni_physical_cells;++i)
        near(after(f.index(i,j)),before(f.index(i,j)),0,"Periodicity changed physical state");
}
template<typename T,bool Isothermal> void wall()
{
    Fixture<T> f;const auto cells=f.cells;physical_fluxes(cells);
    constexpr auto type=Isothermal?BoundaryCondition::isothermalWall:BoundaryCondition::AdiabaticWall;
    ApplyWallConditions<E,T,SolverTypes::RANS2D,type>(cells,f.init,T{1});E().fence();
    const auto rho=Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),cells.rho_old);
    const auto pressure=Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),cells.Pressure);
    const auto mx=Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),cells.Momentum_u_old);
    const auto my=Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),cells.Momentum_v_old);
    for(int i=3;i<3+cells.Ni_physical_cells;++i) for(int layer=0;layer<3;++layer)
    {
        const int inside=f.index(i,3+layer),ghost=f.index(i,2-layer);
        near(mx(ghost)/rho(ghost)+mx(inside)/rho(inside),0,tolerance<T>(),"No-slip u");
        near(my(ghost)/rho(ghost)+my(inside)/rho(inside),0,tolerance<T>(),"No-slip v");
        near(pressure(ghost),pressure(inside),tolerance<T>(),"Wall pressure consistency");
        const double ti=pressure(inside)/rho(inside),tg=pressure(ghost)/rho(ghost);
        if constexpr(Isothermal) near((ti+tg)/2,1,tolerance<T>(),"Isothermal face temperature");
        else near(tg,ti,tolerance<T>(),"Adiabatic mirrored temperature");
    }
}
template<typename T> void farfield()
{
    Fixture<T> f;const auto normals=Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),f.cells.Normals);
    const int inside=f.index(3,3+f.cells.Nj_physical_cells-1),ghost=inside+f.cells.Ni_total_cells;
    const double nx=normals(inside,2,0),ny=normals(inside,2,1),ci=std::sqrt(double(f.init.gamma)*1.2);
    for(double mach:{-2.,-.3,.3,2.})
    {
        const double un=mach*ci,u=un*nx-.1*ny,v=un*ny+.1*nx;
        f.uniform(T{1},T(1.2),T(u),T(v));
        const T alpha=T(std::atan2(ny,nx)+(mach<0?std::acos(-1.):0.));
        initialization_data<T> init({T(std::abs(mach)),T(1.4),alpha,T{250},T(.72)});
        ApplyRiemanInvariantConditions<E>(f.cells,init);E().fence();
        const auto rho=Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),f.cells.rho_old);
        const auto p=Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),f.cells.Pressure);
        const auto mx=Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),f.cells.Momentum_u_old);
        const auto my=Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),f.cells.Momentum_v_old);
        const double ug=mx(ghost)/rho(ghost),vg=my(ghost)/rho(ghost),ung=ug*nx+vg*ny,cg=std::sqrt(double(init.gamma*p(ghost)/rho(ghost)));
        require(rho(ghost)>0&&p(ghost)>0,"Nonpositive farfield state");
        if(mach < -1) {
            near(rho(ghost),init.rho_inf,tolerance<T>(),"Supersonic inlet density");near(p(ghost),init.P_inf,tolerance<T>(),"Supersonic inlet pressure");
        } else if(mach > 1) {
            near(ug,u,tolerance<T>(),"Supersonic outlet u");near(vg,v,tolerance<T>(),"Supersonic outlet v");near(p(ghost),1.2,tolerance<T>(),"Supersonic outlet pressure");
        } else {
            near(ung+2*cg/(init.gamma-1),un+2*ci/(init.gamma-1),tolerance<T>(),"Outgoing acoustic invariant");
            near(ung-2*cg/(init.gamma-1),double(init.velocity_u_inf)*nx+double(init.velocity_v_inf)*ny-2*init.c_inf/(init.gamma-1),tolerance<T>(),"Incoming acoustic invariant");
            const double entropy=p(ghost)/std::pow(double(rho(ghost)),double(init.gamma));
            near(entropy,ung>0?1.2:1.,tolerance<T>(),"Farfield entropy donor");
        }
    }
}
template<typename T> void jst()
{
    Fixture<T> f(true);const auto cells=f.cells;affine_fields(cells);physical_fluxes(cells);
    ApplyWallConditions<E,T,SolverTypes::RANS2D,BoundaryCondition::isothermalWall>(cells,f.init,T{1});
    Kokkos::deep_copy(cells.Lambdas,T{2});
    Kokkos::View<T***> result("test_jst",cells.N_physical_cells,4,4);
    Kokkos::parallel_for("TestJST",Kokkos::RangePolicy<E>(0,cells.N_physical_cells),KOKKOS_LAMBDA(int p){
        const int k=Cell_index(p%cells.Ni_physical_cells+3,p/cells.Ni_physical_cells+3,cells.Ni_total_cells);
        for(int d=0;d<4;++d) {
            const auto dissipation=eval_ArtDissipation<T>(cells,k,d);
            for(int c=0;c<4;++c) result(p,d,c)=dissipation[c];
        }
    });
    const auto host=Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),result);
    const auto pressure=Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),cells.Pressure);
    const auto rho=Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),cells.rho_old);
    const auto mx=Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),cells.Momentum_u_old);
    const auto my=Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),cells.Momentum_v_old);
    const auto energy=Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),cells.E_old);
    const int offsets[4]={-cells.Ni_total_cells,1,cells.Ni_total_cells,-1};
    bool nonzero_wall_mass=false;
    for(int p=0;p<cells.N_physical_cells;++p) for(int d=0;d<4;++d)
    {
        const int k=f.index(p%cells.Ni_physical_cells+3,p/cells.Ni_physical_cells+3),step=offsets[d];
        int stencil[6];for(int a=0;a<6;++a) stencil[a]=k+(a-2)*step;
        long double sensor=0;
        for(int a=1;a<5;++a) {
            const long double prev=pressure(stencil[a-1]),center=pressure(stencil[a]),next=pressure(stencil[a+1]);
            sensor=std::max(sensor,std::abs(next-2*center+prev)/(next+2*center+prev));
        }
        const long double second=.6L*sensor,fourth=std::max(0.L,1.L/64-second);
        for(int c=0;c<4;++c) {
            long double w[4];
            for(int a=0;a<4;++a) {const int index=stencil[a+1];w[a]=c==0?rho(index):c==1?mx(index):c==2?my(index):energy(index);}
            const long double expected=2*(second*(w[2]-w[1])-fourth*(w[3]-3*w[2]+3*w[1]-w[0]));
            near(host(p,d,c),double(expected),tolerance<T>(),"JST reference including boundary faces");
        }
        if(p<cells.Ni_physical_cells&&d==0&&std::abs(double(host(p,d,0)))>1e-8) nonzero_wall_mass=true;
    }
    require(nonzero_wall_mass,"JST was unexpectedly suppressed at the isothermal wall");
}

template<typename T> void viscous_flux()
{
    Fixture<T> f(true);const auto cells=f.cells;affine_fields(cells);physical_fluxes(cells);
    const auto init=f.init;
    Kokkos::View<T***> results("test_flux",cells.N_physical_cells,4,12);
    Kokkos::parallel_for("TestFlux",Kokkos::RangePolicy<E>(0,cells.N_physical_cells),KOKKOS_LAMBDA(int p){
        const int k=Cell_index(p%cells.Ni_physical_cells+3,p/cells.Ni_physical_cells+3,cells.Ni_total_cells);
        for(int d=0;d<4;++d) {
            const auto visc=eval_viscuousDissipation(cells,init,k,d);
            const auto flux=eval_Fluxes<T,SolverTypes::RANS2D>(cells,init,k,d);
            const auto inviscid=eval_Fluxes<T,SolverTypes::Euler2D>(cells,init,k,d);
            for(int c=0;c<4;++c) {results(p,d,c)=visc[c];results(p,d,4+c)=flux[c];results(p,d,8+c)=inviscid[c]-flux[c];}
        }
    });
    const auto result=Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),results);
    const auto normals=Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),cells.Normals);
    const auto positions=Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),cells.Cell_Centers);
    for(int p=0;p<cells.N_physical_cells;++p) for(int d=0;d<4;++d)
    {
        const int i=p%cells.Ni_physical_cells,j=p/cells.Ni_physical_cells,k=f.index(i+3,j+3);
        const double nx=normals(k,d,0),ny=normals(k,d,1),mu=init.mu_inf;
        near(result(p,d,1),mu*(.1*nx+.05*ny),tolerance<T>(),"Viscous x traction");
        near(result(p,d,2),mu*(.05*nx+.2*ny),tolerance<T>(),"Viscous y traction");
        const bool interior=(d==0&&j>0)||(d==1&&i+1<cells.Ni_physical_cells)||(d==2&&j+1<cells.Nj_physical_cells)||(d==3&&i>0);
        // Both states, metrics, normals and fluxes are prepared on these interior faces.
        if(interior) for(int c=0;c<4;++c) near(result(p,d,8+c),result(p,d,c),tolerance<T>(),"Viscous flux sign");
        if((d==1&&i+1<cells.Ni_physical_cells)||(d==2&&j+1<cells.Nj_physical_cells))
        {
            const int neighbour=p+(d==1?1:cells.Ni_physical_cells);
            for(int c=0;c<4;++c) near(result(p,d,4+c)+result(neighbour,(d+2)%4,4+c),0,tolerance<T>(),"Shared-face conservation");
        }
    }
}
template<typename T> void time_step()
{
    Fixture<T> f;const auto cells=f.cells;Kokkos::deep_copy(cells.Lambdas,T{1});
    const int k=f.index(5,5);Kokkos::View<T*> out("test_dt",3);const auto init=f.init;
    Kokkos::parallel_for("TestStep",1,KOKKOS_LAMBDA(int){out(0)=eval_TimeStep<T,SolverTypes::Euler2D>(cells,init,k,T(.5));out(1)=eval_TimeStep<T,SolverTypes::RANS2D>(cells,init,k,T(.5));out(2)=eval_TimeStep<T,SolverTypes::RANS2D>(cells,init,k,T(.25));});
    const auto result=Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),out);
    const auto area=Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),cells.Cell_Areas);
    near(result(0),.5*area(k)/2,tolerance<T>(),"Convective time step");
    require(result(1)>0&&result(1)<result(0),"Viscous restriction missing");near(result(2),result(1)/2,tolerance<T>(),"CFL scaling");
}
template<typename T> void constructors()
{
    const auto mesh=make_mesh<T>();
    using Viscous=SolverFVM<T,SolverTypes::RANS2D>;
    SolverFVM<T> euler({T(.2),T(1.4),T{0}},mesh);
    Viscous iso({T(.2),T(1.4),T{0},T{250},T(.72),BoundaryCondition::isothermalWall,T(.9)},mesh);
    require(iso.get_InitializationData().wall_condition==BoundaryCondition::isothermalWall,"Wall selector lost");
    near(iso.get_InitializationData().wall_temperature,.9,tolerance<T>(),"Wall temperature lost");
    near(iso.get_InitializationData().mu_inf,double(T(.2)*std::sqrt(T(1.4))/T{250}),tolerance<T>(),"Reynolds not propagated");
    invalid([&]{Viscous bad({T(.2),T(1.4),T{0},T{0},T(.72)},mesh);},"Zero Reynolds accepted");
    invalid([&]{Viscous bad({T(.2),T(1.4),T{0},T{250},T{0}},mesh);},"Zero Prandtl accepted");
    invalid([&]{Viscous bad({T(.2),T{1},T{0},T{250},T(.72)},mesh);},"Invalid gamma accepted");
    invalid([&]{Viscous bad({T(.2),T(1.4),T{0},T{250},T(.72),BoundaryCondition::isothermalWall,T{0}},mesh);},"Nonpositive wall temperature accepted");
    invalid([&]{euler.run(-1);},"Negative iteration limit accepted");
}
template<typename T,bool Carpenter> void rk()
{
    const auto mesh=make_mesh<T>();
    constexpr auto scheme=Carpenter?RungeKuttaScheme::CarpenterKennedy:RungeKuttaScheme::Jameson;
    SolverFVM<T,SolverTypes::RANS2D,TimeEvolutionScheme::UnSteady,scheme> solver({T(.2),T(1.4),T{0},T{250},T(.72)},mesh);
    std::ostringstream history;
    solver.run_until(1e-6,T(.05),1000,&history,{},1,100);
    near(solver.get_StateInfo().physical_time,1e-6,1e-12,"Final time not reached");
    require(solver.get_StateInfo().iteration>0,"No physical step executed");
    const auto diagnostics=solver.eval_diagnostics(solver.get_CellsData());
    require(diagnostics.invalid_cells==0&&std::isfinite(diagnostics.residual_rms),"Invalid laminar smoke state");
}
template<typename T> void export_solution()
{
    Fixture<T> f;const auto cells=f.cells;
    Kokkos::parallel_for("ExportPattern",Kokkos::RangePolicy<E>(0,cells.N_total_cells),KOKKOS_LAMBDA(int k){cells.rho_old(k)=T(1)+T(.001)*k;});
    const auto path=std::filesystem::path(std::is_same_v<T,float>?"test-export-float.vts":"test-export-double.vts");
    IO::Export_Solution<T>(path.string(),f.mesh,cells,f.init,7,.25,.01);
    std::ifstream input(path,std::ios::binary);const std::string bytes((std::istreambuf_iterator<char>(input)),{});
    const std::string marker="<AppendedData encoding=\"raw\">_";const auto pos=bytes.find(marker);
    require(pos!=std::string::npos,"Missing VTS binary payload");
    std::uint64_t length=0;std::memcpy(&length,bytes.data()+pos+marker.size(),sizeof(length));
    require(length==std::uint64_t(cells.N_physical_cells)*sizeof(T),"Export included ghost cells");
    for(int j=0;j<cells.Nj_physical_cells;++j) for(int i=0;i<cells.Ni_physical_cells;++i)
    {
        T actual;const int p=i+cells.Ni_physical_cells*j;
        std::memcpy(&actual,bytes.data()+pos+marker.size()+sizeof(length)+p*sizeof(T),sizeof(T));
        near(actual,T(1)+T(.001)*f.index(i+3,j+3),tolerance<T>(),"Export indexing changed");
    }
    input.close();std::filesystem::remove(path);
}
template<typename T> void dispatch(const std::string& name)
{
    if(name=="indexing") indexing<T>();else if(name=="thermodynamics") thermodynamics<T>();
    else if(name=="geometry") geometry<T>();else if(name=="gradients") gradients<T>();
    else if(name=="constitutive") constitutive<T>();else if(name=="periodic") periodic<T>();
    else if(name=="wall_adiabatic") wall<T,false>();else if(name=="wall_isothermal") wall<T,true>();
    else if(name=="farfield") farfield<T>();else if(name=="viscous_flux") viscous_flux<T>();
    else if(name=="jst") jst<T>();
    else if(name=="time_step") time_step<T>();else if(name=="constructors") constructors<T>();
    else if(name=="rk_jameson") rk<T,false>();else if(name=="rk_ck") rk<T,true>();
    else if(name=="export") export_solution<T>();else throw std::invalid_argument("Unknown test case: "+name);
}
int main(int argc,char** argv)
{
    Kokkos::ScopeGuard guard(argc,argv);
    try {
        require(argc>=3,"Expected a test case and scalar precision");
        if(std::string(argv[2])=="float") dispatch<float>(argv[1]);
        else if(std::string(argv[2])=="double") dispatch<double>(argv[1]);
        else throw std::invalid_argument("Unknown precision");
        std::cout<<argv[1]<<" ("<<argv[2]<<") passed\n";return 0;
    } catch(const std::exception& error) {std::cerr<<error.what()<<'\n';return 1;}
}
