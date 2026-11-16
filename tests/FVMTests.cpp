#include "FVMTestSupport.h"
#include "FileIO.h"
#include "RunSimulation.h"
#include <array>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <numbers>
#include <sstream>
#include <utility>
#include <vector>

using namespace Tests;

template<typename T> void affine_fields(const CellsData<T>& cells,const initialization_data<T>& init)
{
    Kokkos::parallel_for("TestAffineFields",Kokkos::RangePolicy<E>(0,cells.N_total_cells),KOKKOS_LAMBDA(int k){
        const T x=cells.Cell_Centers(k,0),y=cells.Cell_Centers(k,1);
        const T u=T(.3)+T(.2)*x-T(.1)*y,v=T(-.2)+T(.15)*x+T(.25)*y;
        const T theta=T(1.5)+T(.05)*x+T(.1)*y;
        cells.rho_old(k)=T{1};cells.Momentum_u_old(k)=u;cells.Momentum_v_old(k)=v;
        cells.Pressure(k)=theta;cells.SoundSpeed(k)=Kokkos::sqrt(T(1.4)*theta);
        cells.Temperature(k)=theta;cells.Viscosity(k)=Solver::NumOps::eval_viscosity(theta,init.T_inf,init.mu_inf);
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
// Inviscid fluxes of the ghost rows touching the wall and the farfield, through the face they share with the domain.
template<typename T> void boundary_ghost_fluxes(const CellsData<T>& cells)
{
    Kokkos::parallel_for("TestGhostFluxes",Kokkos::RangePolicy<E>(0,cells.Ni_physical_cells),KOKKOS_LAMBDA(int i){
        Update_DataField<T,DataFields::Fluxes>(cells,Cell_index(i+3,2,cells.Ni_total_cells),2);
        Update_DataField<T,DataFields::Fluxes>(cells,Cell_index(i+3,3+cells.Nj_physical_cells,cells.Ni_total_cells),0);
    });
}
// Metrics stored in the ghost row j for its face d (the one shared with the domain),
// against a reference rebuilt from the cell centres as seen from the ghost cell.
template<typename T> void ghost_metrics(const Fixture<T>& f,int j,int d,const std::string& message)
{
    const auto centers=Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),f.cells.Cell_Centers);
    const auto metrics=Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),f.cells.DJacobien);
    const int offset[4]={-f.cells.Ni_total_cells,1,f.cells.Ni_total_cells,-1};
    for(int i=3;i<3+f.cells.Ni_physical_cells;++i)
    {
        const int k=f.index(i,j),other=k+offset[d],p=offset[(d+1)%4],m=offset[(d+3)%4];
        long double a[2],b[2];
        for(int c=0;c<2;++c) {a[c]=(long double)centers(other,c)-centers(k,c);b[c]=((long double)centers(k+p,c)-centers(k+m,c)+centers(other+p,c)-centers(other+m,c))/4;}
        const long double det=a[0]*b[1]-a[1]*b[0],ref[4]={a[0]/det,a[1]/det,b[0]/det,b[1]/det};
        for(int c=0;c<4;++c)
        {
            near(metrics(k,d,c),double(ref[c]),tolerance<T>(),message);
            near(metrics(k,d,c)+metrics(other,(d+2)%4,c),0,tolerance<T>(),message+" (opposite of the physical cell)");
        }
    }
}
// Closed O-mesh around a circular wall of radius R (node ni coincides with node 0, periodic in i).
// i runs clockwise so that the cells have a positive area, like make_mesh.
template<typename T> Mesh::MeshData_Structured<T> make_circle_mesh(int ni,int nj,double R,double dr)
{
    Mesh::MeshData_Structured<T> mesh;
    mesh.Ni=ni+1;mesh.Nj=nj+1;mesh.N_nodes=(ni+1)*(nj+1);mesh.N_Cells=ni*nj;
    Kokkos::View<T*> x("circle_x",mesh.N_nodes),y("circle_y",mesh.N_nodes);
    auto hx=Kokkos::create_mirror_view(x),hy=Kokkos::create_mirror_view(y);
    const double step=-2*std::numbers::pi/ni;
    for(int j=0;j<=nj;++j) for(int i=0;i<=ni;++i)
    {
        const int k=i+(ni+1)*j;const double theta=step*(i%ni),r=R+dr*j;
        hx(k)=T(r*std::cos(theta));hy(k)=T(r*std::sin(theta));
    }
    Kokkos::deep_copy(x,hx);Kokkos::deep_copy(y,hy);mesh.X=x;mesh.Y=y;
    return mesh;
}
// Wall metrics on a circular wall: x_eps is tangent, x_eps_eps points to the centre (periodic wrap included),
// x_eta joins the adjacent ghost centre to the wall cell centre and is normal to the wall, pointing outwards.
template<typename T> void wall_metrics()
{
    const int ni=16,nj=4;const double R=1,dr=.1,step=-2*std::numbers::pi/ni;
    const auto mesh=make_circle_mesh<T>(ni,nj,R,dr);CellsData<T> cells(mesh);
    eval_GeometricalData<E>(cells,mesh);E().fence();
    require(cells.Wall_Metrics.extent(0)==size_t(cells.Ni_total_cells)&&cells.Wall_Metrics.extent(1)==6,"Wrong wall metrics shape");
    Kokkos::View<T**> read("read_wall_metrics",cells.Ni_physical_cells,6);
    Kokkos::parallel_for("ReadWallMetrics",Kokkos::RangePolicy<E>(0,cells.Ni_physical_cells),KOKKOS_LAMBDA(int i){
        // Any cell of the column reads the metrics of its wall face: use the second physical row.
        const int k=Cell_index(i+3,4,cells.Ni_total_cells);
        const auto a=extract_DataField<T,DataFields::Wall_xy_epsilon>(cells,k);
        const auto b=extract_DataField<T,DataFields::Wall_xy_epsilon_epsilon>(cells,k);
        const auto c=extract_DataField<T,DataFields::Wall_xy_eta>(cells,k);
        read(i,0)=a[0];read(i,1)=a[1];read(i,2)=b[0];read(i,3)=b[1];read(i,4)=c[0];read(i,5)=c[1];
    });
    const auto h=Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),read);
    const auto centers=Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),cells.Cell_Centers);
    const double rho_m=R*std::cos(step/2); // Radius of the wall face midpoints.
    for(int i=0;i<ni;++i)
    {
        const double phi=step*(i+.5),c=std::cos(phi),s=std::sin(phi);
        near(h(i,0),-rho_m*std::sin(step)*s,tolerance<T>(),"Wall x_eps, x");
        near(h(i,1),rho_m*std::sin(step)*c,tolerance<T>(),"Wall x_eps, y");
        near(h(i,2),2*rho_m*(std::cos(step)-1)*c,tolerance<T>(),"Wall x_eps_eps, x");
        near(h(i,3),2*rho_m*(std::cos(step)-1)*s,tolerance<T>(),"Wall x_eps_eps, y");
        const int k=Cell_index(i+3,3,cells.Ni_total_cells);
        near(h(i,4),2*(centers(k,0)-rho_m*c),tolerance<T>(),"Wall x_eta, x");
        near(h(i,5),2*(centers(k,1)-rho_m*s),tolerance<T>(),"Wall x_eta, y");
        near(-h(i,4)*s+h(i,5)*c,0,tolerance<T>(),"Wall x_eta is not normal to the wall");
        require(h(i,4)*c+h(i,5)*s>0,"Wall x_eta points into the wall");
    }
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
        // The geometry alone must provide the ghost-side metrics of the wall and farfield faces.
        ghost_metrics(f,2,2,"Wall ghost face metric");
        ghost_metrics(f,3+f.cells.Nj_physical_cells,0,"Farfield ghost face metric");
    }
}
template<typename T> void gradients()
{
    for(bool curved:{false,true})
    {
        Fixture<T> f(curved);const auto cells=f.cells;affine_fields(cells,f.init);
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
    // Nondimensional Sutherland law: mu(theta) = mu_inf*theta^1.5*(1+S)/(theta+S), S = 110.4/T_inf. Compared in units of mu_inf.
    const T mu_inf=T(.003),T_inf=T(288.15);
    near(Solver::NumOps::eval_viscosity(T{1},T_inf,mu_inf)/mu_inf,1.,tolerance<T>(),"Sutherland: mu(theta_inf) is not mu_inf");
    near(Solver::NumOps::eval_viscosity(T{2},T_inf,mu_inf)/mu_inf,1.6415751136851524,tolerance<T>(),"Sutherland reference value");
    // Independent reference: ratio of the dimensional law mu = C1*T^1.5/(T+110.4) at T_inf*theta and at T_inf.
    const auto dimensional=[](double t){return 1.458e-6*t*std::sqrt(t)/(t+110.4);};
    for(double theta:{.5,1.3,2.,4.})
        near(Solver::NumOps::eval_viscosity(T(theta),T_inf,mu_inf)/mu_inf,dimensional(theta*288.15)/dimensional(288.15),tolerance<T>(),"Sutherland ratio");
    near(Solver::NumOps::eval_thermalconductivity(mu_inf,T(1.4),T(.72))/mu_inf,4.861111111111111,tolerance<T>(),"Thermal conductivity");
}
template<typename T> void periodic()
{
    Fixture<T> f;const auto cells=f.cells;
    Kokkos::parallel_for("TestPeriodicSources",Kokkos::RangePolicy<E>(0,cells.N_total_cells),KOKKOS_LAMBDA(int k){
        cells.rho_old(k)=T(1)+T(.001)*k;cells.Temperature(k)=T(2)+T(.001)*k;cells.Viscosity(k)=T(.01)+T(1e-5)*k;});
    // Density, and the derived temperature and viscosity, which the periodic copy must carry like the sound speed.
    const std::pair<Kokkos::View<T*>,std::string> fields[3]={{cells.rho_old,"density"},{cells.Temperature,"temperature"},{cells.Viscosity,"viscosity"}};
    std::vector<decltype(Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),cells.rho_old))> before;
    for(const auto& field:fields) before.push_back(Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),field.first));
    ApplyPeriodicity<E>(cells);
    for(int n=0;n<3;++n)
    {
        const auto after=Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),fields[n].first);
        const auto& name=fields[n].second;
        for(int j=0;j<cells.Nj_total_cells;++j) for(int layer=0;layer<3;++layer)
        {
            const int left=f.index(layer,j),right=f.index(cells.Ni_physical_cells+3+layer,j);
            near(after(left),before[n](left+cells.Ni_physical_cells),0,"Periodic left including corners ("+name+")");
            near(after(right),before[n](f.index(3+layer,j)),0,"Periodic right including corners ("+name+")");
        }
        for(int j=3;j<3+cells.Nj_physical_cells;++j) for(int i=3;i<3+cells.Ni_physical_cells;++i)
            near(after(f.index(i,j)),before[n](f.index(i,j)),0,"Periodicity changed physical state ("+name+")");
    }
}
// The derived fields refreshed from a conservative state: theta = p/rho and the Sutherland viscosity mu(theta).
template<typename T> void derived_fields()
{
    Fixture<T> f;const auto cells=f.cells;const auto init=f.init;
    const T rho=T(.8),u=T(.3),v=T(-.1),p=T(1.4);
    const int k=f.index(5,4);
    const Kokkos::Array<T,4> state={rho,rho*u,rho*v,Solver::NumOps::eval_Energy<T>(rho,u,v,p,init.gamma)};
    Kokkos::parallel_for("TestDerivedFields",1,KOKKOS_LAMBDA(int){Update_DataField<T,DataFields::Conservative_Old>(state,cells,init,k);});
    const auto temperature=Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),cells.Temperature);
    const auto viscosity=Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),cells.Viscosity);
    const double theta=double(p)/double(rho);
    near(temperature(k),theta,tolerance<T>(),"Temperature is not p/rho");
    const double mu_ref=Solver::NumOps::eval_viscosity<double>(theta,double(init.T_inf),double(init.mu_inf));
    near(viscosity(k)/init.mu_inf,mu_ref/init.mu_inf,tolerance<T>(),"Viscosity is not mu(theta)");
    // Initial conditions: the free stream, theta_inf = 1 and mu = mu_inf.
    ApplyInitialConditions<E>(cells,init);E().fence();
    const auto t0=Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),cells.Temperature);
    const auto mu0=Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),cells.Viscosity);
    for(int n=0;n<cells.N_total_cells;++n)
    {
        near(t0(n),1.,tolerance<T>(),"Initial temperature");
        near(mu0(n)/init.mu_inf,1.,tolerance<T>(),"Initial viscosity");
    }
    // Inviscid conditions built without the constructor's validation: mu_inf = 0 gives exactly zero, even with a
    // T_inf that would make Sutherland NaN (S = 110.4/0).
    FlowConditions<T> euler_flow{T(.2),T(1.4)};euler_flow.T_inf=T{0};
    const initialization_data<T> euler_init(euler_flow);
    Kokkos::parallel_for("TestInviscidViscosity",1,KOKKOS_LAMBDA(int){Update_DataField<T,DataFields::Conservative_Old>(state,cells,euler_init,k);});
    const auto mu_euler=Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),cells.Viscosity);
    require(mu_euler(k)==T{0},"Inviscid viscosity is not exactly zero");
}
template<typename T> void wall()
{
    Fixture<T> f;const auto cells=f.cells;physical_fluxes(cells);
    ApplyWallConditions<E,T,SolverTypes::RANS2D,BoundaryCondition::AdiabaticWall>(cells,f.init);E().fence();
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
        near(tg,ti,tolerance<T>(),"Adiabatic mirrored temperature");
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
    Fixture<T> f(true);const auto cells=f.cells;affine_fields(cells,f.init);physical_fluxes(cells);
    ApplyWallConditions<E,T,SolverTypes::RANS2D,BoundaryCondition::AdiabaticWall>(cells,f.init);
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
    for(int p=0;p<cells.N_physical_cells;++p) for(int d=0;d<4;++d)
    {
        const int k=f.index(p%cells.Ni_physical_cells+3,p/cells.Ni_physical_cells+3),step=offsets[d];
        int stencil[6];for(int a=0;a<6;++a) stencil[a]=k+(a-2)*step;
        long double sensor=0;
        for(int a=1;a<5;++a) {
            const long double prev=pressure(stencil[a-1]),center=pressure(stencil[a]),next=pressure(stencil[a+1]);
            sensor=std::max(sensor,std::abs(next-2*center+prev)/(next+2*center+prev));
        }
        const long double second=.5L*sensor,fourth=std::max(0.L,1.L/64-second);
        for(int c=0;c<4;++c) {
            long double w[4];
            for(int a=0;a<4;++a) {const int index=stencil[a+1];w[a]=c==0?rho(index):c==1?mx(index):c==2?my(index):energy(index);}
            const long double expected=2*(second*(w[2]-w[1])-fourth*(w[3]-3*w[2]+3*w[1]-w[0]));
            near(host(p,d,c),double(expected),tolerance<T>(),"JST reference including boundary faces");
        }
    }
}

template<typename T> void viscous_flux()
{
    Fixture<T> f(true);const auto cells=f.cells;const auto init=f.init;
    // The boundary conditions give the ghost rows their normals (their face metrics come from the geometry);
    // the affine field then replaces every state, ghosts included, so the analytic tractions hold up to the boundaries.
    ApplyWallConditions<E,T,SolverTypes::RANS2D,BoundaryCondition::AdiabaticWall>(cells,init);
    ApplyRiemanInvariantConditions<E>(cells,init);E().fence();
    affine_fields(cells,f.init);physical_fluxes(cells);boundary_ghost_fluxes(cells);
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
    const auto pressure=Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),cells.Pressure);
    const auto rho=Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),cells.rho_old);
    // Sutherland viscosity of a cell, from its temperature theta = p/rho (the affine field makes theta vary in space).
    const auto mu_at=[&](int k){return Solver::NumOps::eval_viscosity<double>(double(pressure(k))/double(rho(k)),double(init.T_inf),double(init.mu_inf));};
    const int offset[4]={-cells.Ni_total_cells,1,cells.Ni_total_cells,-1};
    for(int p=0;p<cells.N_physical_cells;++p) for(int d=0;d<4;++d)
    {
        const int i=p%cells.Ni_physical_cells,j=p/cells.Ni_physical_cells,k=f.index(i+3,j+3);
        const double nx=normals(k,d,0),ny=normals(k,d,1),mu=mu_at(k),mu_side=mu_at(k+offset[d]);
        require(std::abs(mu/init.mu_inf-1)>1e-2,"Test field too close to theta_inf: Sutherland would go unnoticed");
        near(result(p,d,1)/mu,.1*nx+.05*ny,tolerance<T>(),"Viscous x traction");
        near(result(p,d,2)/mu,.05*nx+.2*ny,tolerance<T>(),"Viscous y traction");
        // Wall (d==0, j==0) and farfield (d==2, last j) faces are included: the full viscous flux must cross them, not half of it.
        // Only the i-boundaries are left out: this test grid is not periodic.
        const bool assembled=d==0||d==2||(d==1&&i+1<cells.Ni_physical_cells)||(d==3&&i>0);
        const bool boundary=(d==0&&j==0)||(d==2&&j+1==cells.Nj_physical_cells);
        // Each side evaluates the face flux with its own viscosity and the face averages them: with exact (affine) gradients
        // the assembled viscous flux is the cell's own one scaled by (mu + mu_side) / (2 mu).
        const double face_scale=(mu+mu_side)/(2*mu);
        // Compared in units of mu: the viscous fluxes are of order mu, well below the tolerance of an unscaled comparison in float.
        if(assembled) for(int c=0;c<4;++c) near(result(p,d,8+c)/mu,face_scale*result(p,d,c)/mu,20*tolerance<T>(),boundary?"Boundary viscous flux is not the full flux":"Viscous flux sign");
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
    // Viscous part of the spectral radius, recovered from the step: mu(theta)/rho*(4/3+gamma/Pr)*sum(L^2)/A, with the local
    // Sutherland viscosity of the uniform fixture (theta = 1.2, rho = 1), about 15 % above mu_inf.
    const auto normals=Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),cells.Normals);
    double sum_l2=0;for(int d=0;d<4;++d) sum_l2+=double(normals(k,d,2))*double(normals(k,d,2));
    const double mu=Solver::NumOps::eval_viscosity<double>(1.2,double(init.T_inf),double(init.mu_inf));
    const double viscous_expected=mu*(4./3+double(init.gamma)/double(init.Pr_inf))*sum_l2/double(area(k));
    const double viscous_measured=.5*double(area(k))/double(result(1))-2;
    near(viscous_measured/viscous_expected,1.,tolerance<T>(),"Viscous time step does not use the local viscosity");
}
template<typename T> void constructors()
{
    const auto mesh=make_mesh<T>();
    using Viscous=SolverFVM<T,SolverTypes::RANS2D>;
    SolverFVM<T> euler({T(.2),T(1.4),T{0}},mesh);
    Viscous viscous({T(.2),T(1.4),T{0},T{250},T(.72),BoundaryCondition::AdiabaticWall},mesh);
    require(viscous.get_InitializationData().wall_condition==BoundaryCondition::AdiabaticWall,"Wall selector lost");
    near(viscous.get_InitializationData().mu_inf,double(T(.2)*std::sqrt(T(1.4))/T{250}),tolerance<T>(),"Reynolds not propagated");
    // T_inf comes last: the positional initialisations above keep the default free-stream temperature.
    near(viscous.get_InitializationData().T_inf/288.15,1.,tolerance<T>(),"Default free-stream temperature lost");
    invalid([&]{Viscous bad({T(.2),T(1.4),T{0},T{0},T(.72)},mesh);},"Zero Reynolds accepted");
    invalid([&]{Viscous bad({T(.2),T(1.4),T{0},T{250},T{0}},mesh);},"Zero Prandtl accepted");
    invalid([&]{Viscous bad({T(.2),T{1},T{0},T{250},T(.72)},mesh);},"Invalid gamma accepted");
    invalid([&]{Viscous bad({T(.2)},mesh);},"Viscous solver accepted the default infinite Reynolds");
    // Sutherland divides by T_inf (S = 110.4/T_inf): a nonpositive, infinite or NaN free-stream temperature is refused,
    // by Euler2D too, which does not use it.
    FlowConditions<T> cold{T(.2),T(1.4),T{0},T{250},T(.72)};
    for(T t_inf:{T{0},T(-100),std::numeric_limits<T>::infinity(),std::numeric_limits<T>::quiet_NaN()})
    {
        cold.T_inf=t_inf;
        invalid([&]{Viscous bad(cold,mesh);},"Invalid free-stream temperature accepted by RANS2D");
        invalid([&]{SolverFVM<T> bad(cold,mesh);},"Invalid free-stream temperature accepted by Euler2D");
    }
    // Euler2D: zero viscosity everywhere, ghosts included, after the initial state and after one iteration.
    cold.T_inf=T(288.15);SolverFVM<T> inviscid_run(cold,mesh);
    const auto zero_viscosity=[&](const char* when){
        const auto mu=Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),inviscid_run.get_CellsData().Viscosity);
        for(std::size_t k=0;k<mu.extent(0);++k) require(mu(k)==T{0},std::string("Euler2D viscosity is not zero ")+when);
    };
    zero_viscosity("after initialisation");
    inviscid_run.run(1,T{0});zero_viscosity("after one iteration");
    FlowConditions<T> named;named.Mach_inf=T(.3);named.Re_inf=T{100};
    Viscous by_name(named,mesh);
    near(by_name.get_InitializationData().gamma,1.4,tolerance<T>(),"Default gamma lost");
    near(by_name.get_InitializationData().Pr_inf,1,tolerance<T>(),"Default Prandtl lost");
    require(by_name.get_InitializationData().wall_condition==BoundaryCondition::AdiabaticWall,"Default wall is not adiabatic");
    near(by_name.get_InitializationData().velocity_u_inf,double(T(.3))*std::sqrt(1.4),tolerance<T>(),"Free-stream velocity not derived");
    invalid([&]{euler.run(-1,T{0});},"Negative iteration limit accepted");
    // Residual tolerance: NaN, negative or infinite values would silently disable the stop, they are refused.
    for(T tol:{std::numeric_limits<T>::quiet_NaN(),T(-1e-6),std::numeric_limits<T>::infinity()})
        invalid([&]{euler.run(1,tol);},"Invalid residual tolerance accepted");
    // tol = 0 never stops on the residual; a tolerance above the initial residual stops before the first step.
    SolverFVM<T> no_stop({T(.2),T(1.4),T{0}},mesh);no_stop.run(3,T{0});
    require(no_stop.get_StateInfo().iteration==3,"tol = 0 stopped on the residual");
    SolverFVM<T> early_stop({T(.2),T(1.4),T{0}},mesh);early_stop.run(3,T(1e30));
    require(early_stop.get_StateInfo().iteration==0,"A residual below the tolerance did not stop the steady run");
}
// Density RMS residual (dRHO/dT), its history rows, and the steady stop on both residuals.
template<typename T> void residuals()
{
    const auto mesh=make_mesh<T>();
    SolverFVM<T> solver({T(.2),T(1.4),T{0}},mesh);
    solver.run(2,T{0});
    // Reference: per-cell residual R/A of every component, reduced on the host.
    const auto cells=solver.get_CellsData();const auto init=solver.get_InitializationData();
    Kokkos::View<double**> per_cell("test_residuals",cells.N_physical_cells,4);
    Kokkos::parallel_for("TestResiduals",Kokkos::RangePolicy<E>(0,cells.N_physical_cells),KOKKOS_LAMBDA(int p){
        const int k=Cell_index(p%cells.Ni_physical_cells+3,p/cells.Ni_physical_cells+3,cells.Ni_total_cells);
        const auto fd=eval_All_FluxesDissipation<T,SolverTypes::Euler2D>(cells,init,k);
        for(int c=0;c<4;++c) per_cell(p,c)=double(fd.first[c]-fd.second[c])/double(cells.Cell_Areas(k));
    });
    const auto host=Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),per_cell);
    double all=0,rho=0;
    for(int p=0;p<cells.N_physical_cells;++p) {for(int c=0;c<4;++c) all+=host(p,c)*host(p,c);rho+=host(p,0)*host(p,0);}
    const double rms_ref=std::sqrt(all/(4.*cells.N_physical_cells)),rho_ref=std::sqrt(rho/cells.N_physical_cells);
    const auto d=solver.eval_diagnostics(cells);
    require(rho_ref>0,"Test state has no density residual");
    near(d.residual_rms_rho/rho_ref,1.,1e-12,"Density RMS residual");
    near(d.residual_rms/rms_ref,1.,1e-12,"Global RMS residual");

    // History: a row every history_interval iterations, plus the initial and the last state; rms_rho in the third column.
    SolverFVM<T> logged({T(.2),T(1.4),T{0}},mesh);std::ostringstream history;
    logged.run(5,T{0},T{1},1000,2,&history);
    std::istringstream lines(history.str());std::string line;std::vector<std::string> rows;
    std::getline(lines,line);require(line.rfind("iteration,rms,rms_rho,",0)==0,"History header lacks rms_rho");
    while(std::getline(lines,line)) rows.push_back(line);
    require(rows.size()==4,"Expected history rows at iterations 0, 2, 4 and 5, got "+std::to_string(rows.size()));
    const int expected_iterations[4]={0,2,4,5};
    for(int r=0;r<4;++r) require(std::stoi(rows[r])==expected_iterations[r],"History row at the wrong iteration");
    const auto last=logged.eval_diagnostics(logged.get_CellsData());
    std::istringstream fields(rows.back());std::string iteration,rms,rms_rho;
    std::getline(fields,iteration,',');std::getline(fields,rms,',');std::getline(fields,rms_rho,',');
    near(std::stod(rms_rho)/last.residual_rms_rho,1.,1e-12,"History rms_rho is not the last density residual");
    invalid([&]{logged.run(1,T{0},T{1},1000,0);},"Zero history interval accepted");

    // Steady stop: only the density residual is compared with tol. Just below it does not stop; a tol between
    // the density and the global residuals stops even though the global one is still above it.
    SolverFVM<T> probe({T(.2),T(1.4),T{0}},mesh);
    const auto initial=probe.eval_diagnostics(probe.get_CellsData());
    require(initial.residual_rms>1.01*initial.residual_rms_rho,"Test state needs a global residual above the density one");
    probe.run(1,T(.999*initial.residual_rms_rho));
    require(probe.get_StateInfo().iteration==1,"Stopped while the density residual was still above tol");
    SolverFVM<T> stopped({T(.2),T(1.4),T{0}},mesh);
    stopped.run(1,T(std::sqrt(initial.residual_rms_rho*initial.residual_rms)));
    require(stopped.get_StateInfo().iteration==0,"Did not stop with the density residual below tol");
}
// Grid study: one output directory per mesh inside results/study_NNNN, a failing mesh skipped and listed in study.csv.
template<typename T> void study()
{
    const auto root=std::filesystem::temp_directory_path()/(std::string("cmfsn-study-")+(std::is_same_v<T,float>?"float":"double"));
    std::filesystem::remove_all(root);std::filesystem::create_directories(root);
    // Mesh files in the import format: "Ni Nj", then every X, then every Y.
    const auto write_mesh=[&](const std::string& name,int ni,int nj){
        const auto mesh=make_mesh<T>(false,ni,nj);
        const auto x=Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),mesh.X),y=Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),mesh.Y);
        std::ofstream file(root/name);file<<std::setprecision(17)<<mesh.Ni<<' '<<mesh.Nj<<'\n';
        for(int n=0;n<mesh.N_nodes;++n) file<<x(n)<<'\n';
        for(int n=0;n<mesh.N_nodes;++n) file<<y(n)<<'\n';
        return root/name;
    };
    const std::vector<std::filesystem::path> meshes={write_mesh("a.x",10,8),root/"missing.x",write_mesh("b.x",12,8)};
    CMFSN::SimulationParameters<T> parameters;
    parameters.flow={T(.2),T(1.4),T{0}};parameters.iterations=2;parameters.tol=T{0};parameters.output_root=root/"results";
    const int failures=CMFSN::RunStudy<T,SolverTypes::Euler2D,TimeEvolutionScheme::Steady,RungeKuttaScheme::Jameson>(parameters,meshes);
    require(failures==1,"Expected exactly one failed mesh, got "+std::to_string(failures));
    const auto study_dir=root/"results"/"study_0001";
    for(const char* name:{"a.x","b.x"})
    {
        require(std::filesystem::exists(study_dir/name/"solution.vts"),std::string("Missing solution for ")+name);
        std::ifstream history(study_dir/name/"residuals.csv");std::string line;int rows=-1;
        while(std::getline(history,line)) ++rows;
        require(rows==3,std::string("Expected history rows at iterations 0, 1 and 2 for ")+name);
    }
    require(!std::filesystem::exists(study_dir/"missing.x"/"solution.vts"),"A solution was written for the missing mesh");
    // The failure does not stop the study: the mesh after it was still computed (checked above), and every mesh is listed.
    std::ifstream summary(study_dir/"study.csv");std::string line;std::vector<std::string> rows;
    std::getline(summary,line);require(line=="mesh,directory,status,message","Unexpected study.csv header");
    while(std::getline(summary,line)) rows.push_back(line);
    require(rows.size()==3,"study.csv must list every mesh");
    require(rows[0].find(",ok,")!=std::string::npos&&rows[1].find(",failed,")!=std::string::npos&&rows[2].find(",ok,")!=std::string::npos,"Wrong mesh status in study.csv");
    // A second study goes to the next directory; duplicated mesh names are refused before any computation.
    CMFSN::RunStudy<T,SolverTypes::Euler2D,TimeEvolutionScheme::Steady,RungeKuttaScheme::Jameson>(parameters,{meshes[0]});
    require(std::filesystem::exists(root/"results"/"study_0002"/"a.x"/"solution.vts"),"Second study did not get its own directory");
    invalid([&]{CMFSN::RunStudy<T,SolverTypes::Euler2D,TimeEvolutionScheme::Steady,RungeKuttaScheme::Jameson>(parameters,{meshes[0],meshes[0]});},"Duplicated mesh names accepted");
    summary.close();std::filesystem::remove_all(root);
}
template<typename T,bool Carpenter> void rk()
{
    const auto mesh=make_mesh<T>();
    constexpr auto scheme=Carpenter?RungeKuttaScheme::CarpenterKennedy:RungeKuttaScheme::Jameson;
    SolverFVM<T,SolverTypes::RANS2D,TimeEvolutionScheme::UnSteady,scheme> solver({T(.2),T(1.4),T{0},T{250},T(.72)},mesh);
    std::ostringstream history;
    solver.run_until(1e-6,T(.05),1000,1,&history,{},1,100);
    near(solver.get_StateInfo().physical_time,1e-6,1e-12,"Final time not reached");
    require(solver.get_StateInfo().iteration>0,"No physical step executed");
    const auto diagnostics=solver.eval_diagnostics(solver.get_CellsData());
    require(diagnostics.invalid_cells==0&&std::isfinite(diagnostics.residual_rms),"Invalid laminar smoke state");

    // Stage storage left by a previous step must not leak into the next one.
    Kokkos::deep_copy(solver.get_CellsData().D_history,std::numeric_limits<T>::quiet_NaN());
    solver.run_until(2e-6,T(.05),1000,1,nullptr,{},1,100);
    const auto after=solver.eval_diagnostics(solver.get_CellsData());
    require(after.invalid_cells==0&&std::isfinite(after.residual_rms),"Stale stage storage leaked into the new step");

    if constexpr(Carpenter)
    {
        // One short step from the same state: both schemes are consistent, so their increments must agree to leading order.
        SolverFVM<T,SolverTypes::RANS2D,TimeEvolutionScheme::UnSteady,RungeKuttaScheme::CarpenterKennedy> carpenter({T(.2),T(1.4),T{0},T{250},T(.72)},mesh);
        SolverFVM<T,SolverTypes::RANS2D,TimeEvolutionScheme::UnSteady,RungeKuttaScheme::Jameson> jameson({T(.2),T(1.4),T{0},T{250},T(.72)},mesh);
        const auto fields=[](const CellsData<T>& cells){return std::array{Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),cells.rho_old),Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),cells.Momentum_u_old),Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),cells.Momentum_v_old),Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),cells.E_old)};};
        const auto start=fields(carpenter.get_CellsData());
        carpenter.run(1,T{0},T(.05),1000,1,nullptr,{},1,1e-4);jameson.run(1,T{0},T(.05),1000,1,nullptr,{},1,1e-4);
        near(carpenter.get_StateInfo().time_step,jameson.get_StateInfo().time_step,0,"Schemes did not take the same step");
        const auto ck=fields(carpenter.get_CellsData()),reference=fields(jameson.get_CellsData());
        const auto& cells=carpenter.get_CellsData();
        double increment=0,difference=0;
        for(int c=0;c<4;++c) for(int j=3;j<3+cells.Nj_physical_cells;++j) for(int i=3;i<3+cells.Ni_physical_cells;++i)
        {
            const int k=Cell_index(i,j,cells.Ni_total_cells);
            increment=std::max(increment,std::abs(double(reference[c](k))-double(start[c](k))));
            difference=std::max(difference,std::abs(double(ck[c](k))-double(reference[c](k))));
        }
        require(increment>0,"The reference step did not change the state");
        require(difference<=.05*increment,"Carpenter-Kennedy increment disagrees with Jameson: difference="+std::to_string(difference/increment)+" of the increment");
    }
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
    // Temperature and viscosity, located by name, are rebuilt from the conservative state like the pressure.
    const auto array_value=[&](const std::string& name,int p){
        const auto tag=bytes.find("Name=\""+name+"\"");
        require(tag!=std::string::npos,"Missing exported field "+name);
        const auto at=bytes.find("offset=\"",tag)+8;
        const std::uint64_t offset=std::stoull(bytes.substr(at,bytes.find('"',at)-at));
        T value;std::memcpy(&value,bytes.data()+pos+marker.size()+offset+sizeof(std::uint64_t)+p*sizeof(T),sizeof(T));
        return double(value);
    };
    require(bytes.find("Name=\"t_inf_kelvin\"")!=std::string::npos,"Free-stream temperature missing from the metadata");
    const auto rho=Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),cells.rho_old);
    const auto mx=Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),cells.Momentum_u_old);
    const auto my=Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),cells.Momentum_v_old);
    const auto energy=Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),cells.E_old);
    for(int p=0;p<cells.N_physical_cells;++p)
    {
        const int k=f.index(p%cells.Ni_physical_cells+3,p/cells.Ni_physical_cells+3);
        const double r=rho(k),u=mx(k)/r,v=my(k)/r;
        const double theta=(double(f.init.gamma)-1)*(energy(k)-.5*r*(u*u+v*v))/r;
        near(array_value("temperature",p),theta,tolerance<T>(),"Exported temperature");
        const double mu=Solver::NumOps::eval_viscosity<double>(theta,double(f.init.T_inf),double(f.init.mu_inf));
        near(array_value("viscosity",p)/f.init.mu_inf,mu/f.init.mu_inf,tolerance<T>(),"Exported viscosity");
    }
    input.close();std::filesystem::remove(path);
    // Euler2D has no viscosity: neither the field nor T_inf is exported, the temperature still is.
    const auto euler_path=std::filesystem::path(std::is_same_v<T,float>?"test-export-euler-float.vts":"test-export-euler-double.vts");
    IO::Export_Solution<T>(euler_path.string(),f.mesh,cells,initialization_data<T>(inviscid(FlowConditions<T>(f.init))),7);
    std::ifstream euler_input(euler_path,std::ios::binary);const std::string euler_bytes((std::istreambuf_iterator<char>(euler_input)),{});
    require(euler_bytes.find("Name=\"temperature\"")!=std::string::npos,"Euler export lost the temperature");
    require(euler_bytes.find("Name=\"viscosity\"")==std::string::npos&&euler_bytes.find("t_inf_kelvin")==std::string::npos,"Euler export wrote a viscosity");
    euler_input.close();std::filesystem::remove(euler_path);
}
template<typename T> void dispatch(const std::string& name)
{
    if(name=="indexing") indexing<T>();else if(name=="thermodynamics") thermodynamics<T>();
    else if(name=="geometry") geometry<T>();else if(name=="wall_metrics") wall_metrics<T>();
    else if(name=="gradients") gradients<T>();
    else if(name=="constitutive") constitutive<T>();else if(name=="periodic") periodic<T>();
    else if(name=="derived_fields") derived_fields<T>();
    else if(name=="wall_adiabatic") wall<T>();
    else if(name=="farfield") farfield<T>();else if(name=="viscous_flux") viscous_flux<T>();
    else if(name=="jst") jst<T>();
    else if(name=="time_step") time_step<T>();else if(name=="constructors") constructors<T>();
    else if(name=="residuals") residuals<T>();else if(name=="study") study<T>();
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
