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
    T eval_Lambda(Kokkos::Array<T, 2> Velocity, Kokkos::Array<T, 2> normal, T SoundSpeed, T EdgeLength)
    {
        return (Kokkos::abs(dot(Velocity, normal)) + SoundSpeed) * EdgeLength;
    }

	template <typename T>
	KOKKOS_INLINE_FUNCTION
	T mean(const T& A, const T& B)
	{
		return T(0.5) * (A + B);
	}
	template <typename T, size_t N>
	KOKKOS_INLINE_FUNCTION
	Kokkos::Array<T, N> mean(const Kokkos::Array<T, N>& A, const Kokkos::Array<T, N>& B)
	{
		Kokkos::Array<T, N> mean;
		mean = 0.5 * (A + B);
		return mean;
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
	Kokkos::Array<T, 3> Tau(T mu, const Kokkos::Array<T, 4>& DV)
	{
		T Tau_xx = mu * ((T{ 4 } / T{3})*DV[0] - (T{ 2 } / T{ 3 })*DV[3]);
		T Tau_yy = mu * ((T{ 4 } / T{ 3 }) * DV[3] - (T{ 2 } / T{ 3 }) * DV[0]);
		T Tau_xy = mu * (DV[1] + DV[2]);
		return Kokkos::Array<T, 3> {Tau_xx, Tau_yy, Tau_xy};
	}



}
