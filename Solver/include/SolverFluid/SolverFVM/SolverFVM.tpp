#pragma once

#include "SolverFVM.h"
#include "Kokkos_Core.hpp"


namespace Solver
{
	template <typename T, SolverTypes SType>
	SolverFVM<T, SType>::SolverFVM(T Mach_inf, T gamma, T alpha, int N_gLayer, const MeshData_Structured<T>& Mesh_Euler) : cells_handler(Mach_inf, gamma, alpha, N_gLayer, Mesh_Euler)
	{

	}

	//void init_cells()
	//{

	//}

	//template<typename ExecSpace>
	//void SolverFVM<T>::inite_ghostcells(CellsData<T> cellsdata, const MeshData_Euler<T>& meshdata)
	//{
	//	CellOps::ApplyPeriodicity<ExecSpace>(cellsdata, meshdata);
	//}


}