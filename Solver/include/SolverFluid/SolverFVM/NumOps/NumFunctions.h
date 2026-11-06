#pragma once
#include <Kokkos_Core.hpp>
#include <iostream>
#include <NumOps/ArrayOps.h>

namespace Solver::NumOps
{
	template <typename T>
	KOKKOS_INLINE_FUNCTION
	T eval_pressure(T rho, T rhoE, T u, T v, T gamma)
	{
		T kinetic_energy = 0.5 * rho * (u * u + v * v);
		T pressure = (gamma - 1.0) *(rhoE - kinetic_energy);
		return pressure;
	}
	template <typename T>
	KOKKOS_INLINE_FUNCTION
	T eval_Energy(T rho, T u, T v, T P, T gamma)
	{
		T kinetic_energy = 0.5 * rho * (u * u + v * v);
		T energy = P / (gamma - 1.0) + kinetic_energy;
		return energy;
	}
	template <typename T>
	KOKKOS_INLINE_FUNCTION
	T eval_sound_speed(T gamma, T pressure, T rho)
	{
		T sound_speed = sqrt(gamma * pressure / rho);
		return sound_speed;
	}
	template <typename T>
	KOKKOS_INLINE_FUNCTION
	T eval_velocity_magnitude(T u, T v)
	{
		T velocity = sqrt(u * u + v * v);
		return velocity;
	}
	template <typename T>
	KOKKOS_INLINE_FUNCTION
	T eval_MachNumber(T u, T v, T c)
	{
		T velocity = eval_velocity_magnitude(u, v);
		T MachNumber = velocity / c;
		return MachNumber;
	}
	template <typename T>
	KOKKOS_INLINE_FUNCTION
	T eval_MachNumber(T Vn, T c)
	{
		T MachNumber = Vn / c;
		return MachNumber;
	}
	template <typename T>
	KOKKOS_INLINE_FUNCTION
	T mean(const T& A, const T& B)
	{
		return T(0.5) * (A + B);
	}
	template <typename T, size_t N>
	KOKKOS_INLINE_FUNCTION
	T max(const Kokkos::Array<T, N>& A)
	{
		static_assert(N > 0, "max requires a non-empty array");

		T result = A[0];
		for (size_t i = 1; i < N; ++i)
		{
			if (A[i] > result)
				result = A[i];
		}
		return result;
	}



	template <typename T>
	KOKKOS_INLINE_FUNCTION
	Kokkos::pair< Kokkos::Array<T, 2>, Kokkos::pair<T, T> > eval_RiemannBoundaryState(Kokkos::Array<T, 2> V_in, Kokkos::Array<T, 2> V_inf, Kokkos::Array<T, 2> normal, T c_in, T c_inf, T gamma, T rho_in, T rho_inf)
	{
		T Vn_in = dot(V_in, normal);
		T Vn_inf = dot(V_inf, normal);

		T Mach = eval_MachNumber(Vn_in, c_in);
		T J_plus;
		T J_minus;

		if (Mach < - 1.0)
		{
			J_plus = Vn_inf + 2 * c_inf / (gamma - 1.0);
			J_minus = Vn_inf - 2 * c_inf / (gamma - 1.0);
		}
		else if (Mach < 1)
		{
			J_plus = Vn_in + 2 * c_in / (gamma - 1.0);
			J_minus = Vn_inf - 2 * c_inf / (gamma - 1.0);
		}
		else
		{
			J_plus = Vn_in + 2 * c_in / (gamma - 1.0);
			J_minus = Vn_in - 2 * c_in / (gamma - 1.0);
		}
		T Vbn = (J_plus + J_minus) / T{2};
		T cb = (gamma - T{1}) / T{4} * (J_plus - J_minus);
		T Entropie_k;
		Kokkos::Array<T, 2> V_boundary;

		if (Vbn > T{0})
		{
			Entropie_k = c_in * c_in /(gamma * Kokkos::pow(rho_in, gamma - T{1}));
			V_boundary = V_in + (Vbn - Vn_in) * normal;
		}
		else
		{
			Entropie_k = c_inf * c_inf / (gamma * Kokkos::pow(rho_inf, gamma - T{1}));
			V_boundary = V_inf + (Vbn - Vn_inf) * normal;
		}

		T rho_boundary = Kokkos::pow(cb * cb / (gamma * Entropie_k), T{1} / (gamma - T{1}));
		T pressure_boundary = rho_boundary * cb * cb / gamma;
		return Kokkos::make_pair(V_boundary, Kokkos::make_pair(rho_boundary, pressure_boundary));
	}
	
}
