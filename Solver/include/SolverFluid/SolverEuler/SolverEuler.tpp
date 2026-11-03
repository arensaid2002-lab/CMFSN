#pragma once

#include "SolverEuler.h"
#include "Kokkos_Core.hpp"

using namespace Kokkos;
namespace Solver
{
	template <typename T>
	SolverFVM<T>::SolverFVM(T Mach_inf, T gamma, T alpha,  int N_gLayer, const MeshData_Structured<T>& Mesh_Euler) : cells_euler(N_gLayer, Mesh_Euler)
	{
		gamma_ = gamma;
		Mach_inf_ = Mach_inf;
		alpha_ = alpha; // Radians
		velocity_u_inf = Mach_inf * sqrt(gamma_) * cos(alpha_);
		velocity_v_inf = Mach_inf * sqrt(gamma_) * sin(alpha_);
		
	}

	void init_cells()
	{

	}

	//template<typename ExecSpace>
	//void SolverEuler<T>::inite_ghostcells(CellsData<T> cellsdata, const MeshData_Euler<T>& meshdata)
	//{
	//	CellOps::ApplyPeriodicity<ExecSpace>(cellsdata, meshdata);
	//}


}