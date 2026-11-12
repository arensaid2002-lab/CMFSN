#pragma once

#include <Kokkos_Core.hpp>
#include <iostream>
#include <type_traits>

namespace Solver::NumOps
{

    #pragma region Array_Operations
    template <typename T, size_t N>
    KOKKOS_INLINE_FUNCTION
    Kokkos::Array<T, N>& operator+=(Kokkos::Array<T, N>& A, const Kokkos::Array<T, N>& B)
    {
        for (size_t i = 0; i < N; i++)
        {
            A[i] = A[i] + B[i];
        }
        return A;
    }
    template <typename T, size_t N>
    KOKKOS_INLINE_FUNCTION
    Kokkos::Array<T, N> operator+(const Kokkos::Array<T, N>& A,const Kokkos::Array<T, N>& B)
    {
        Kokkos::Array<T, N> C{};

        for (size_t i = 0; i < N; i++)
        {
            C[i] = A[i] + B[i];
        }
        return C;
    }

    template <typename T, size_t N>
    KOKKOS_INLINE_FUNCTION
    Kokkos::Array<T, N> operator-(const Kokkos::Array<T, N>& A, const Kokkos::Array<T, N>& B)
    {
        Kokkos::Array<T, N> C{};

        for (size_t i = 0; i < N; i++)
        {
            C[i] = A[i] - B[i];
        }
        return C;
    }

    template <typename T, size_t N>
    KOKKOS_INLINE_FUNCTION
    Kokkos::Array<T, N> operator*(const Kokkos::Array<T, N>& A, const Kokkos::Array<T, N>& B)
    {
        Kokkos::Array<T, N> C{};

        for (size_t i = 0; i < N; i++)
        {
            C[i] = A[i] * B[i];
        }
        return C;
    }

    template <typename T, size_t N>
    KOKKOS_INLINE_FUNCTION
    Kokkos::Array<T, N> operator/(const Kokkos::Array<T, N>& A, const Kokkos::Array<T, N>& B)
    {
        Kokkos::Array<T, N> C{};

        for (size_t i = 0; i < N; i++)
        {
            C[i] = A[i]/B[i];
        }
        return C;
    }
    #pragma endregion

    #pragma region Scalar_Array_Operations
    template <typename T, typename S, size_t N>
    requires std::is_arithmetic_v<S>
    KOKKOS_INLINE_FUNCTION
    Kokkos::Array<T, N> operator+(const Kokkos::Array<T, N>& A, S b)
    {
        Kokkos::Array<T, N> C;

        for (size_t i = 0; i < N; i++)
        {
            C[i] = A[i] + b;
        }
        return C;
    }
    template <typename T, typename S, size_t N>
    requires std::is_arithmetic_v<S>
    KOKKOS_INLINE_FUNCTION
    Kokkos::Array<T, N> operator+(S b, const Kokkos::Array<T, N>& A)
    {
        Kokkos::Array<T, N> C;

        for (size_t i = 0; i < N; i++)
        {
            C[i] = A[i] + b;
        }
        return C;
    }

    template <typename T, typename S, size_t N>
    requires std::is_arithmetic_v<S>
    KOKKOS_INLINE_FUNCTION
    Kokkos::Array<T, N> operator-(const Kokkos::Array<T, N>& A, S b)
    {
        Kokkos::Array<T, N> C;

        for (size_t i = 0; i < N; i++)
        {
            C[i] = A[i] - b;
        }
        return C;
    }
    template <typename T, typename S, size_t N>
    requires std::is_arithmetic_v<S>
    KOKKOS_INLINE_FUNCTION
    Kokkos::Array<T, N> operator-(S b, const Kokkos::Array<T, N>& A)
    {
        Kokkos::Array<T, N> C;

        for (size_t i = 0; i < N; i++)
        {
            C[i] = b - A[i];
        }
        return C;
    }

    template <typename T, typename S, size_t N>
    requires std::is_arithmetic_v<S>
    KOKKOS_INLINE_FUNCTION
    Kokkos::Array<T, N> operator *(const Kokkos::Array<T, N>& A, S b)
    {
        Kokkos::Array<T, N> C;

        for (size_t i = 0; i < N; i++)
        {
            C[i] = A[i] * b;
        }
        return C;
    }
    template <typename T, typename S, size_t N>
    requires std::is_arithmetic_v<S>
    KOKKOS_INLINE_FUNCTION
    Kokkos::Array<T, N> operator *(S b, const Kokkos::Array<T, N>& A)
    {
        Kokkos::Array<T, N> C;

        for (size_t i = 0; i < N; i++)
        {
            C[i] = A[i] * b;
        }
        return C;
    }

    template <typename T, typename S, size_t N>
    requires std::is_arithmetic_v<S>
    KOKKOS_INLINE_FUNCTION
    Kokkos::Array<T, N> operator/(const Kokkos::Array<T, N>& A, S b)
    {
        Kokkos::Array<T, N> C;

        for (size_t i = 0; i < N; i++)
        {
            C[i] = A[i] / b;
        }
        return C;
    }

    template <typename T, typename S, size_t N>
    requires std::is_arithmetic_v<S>
    KOKKOS_INLINE_FUNCTION
    Kokkos::Array<T, N> operator/(S b, const Kokkos::Array<T, N>& A)
    {
        Kokkos::Array<T, N> C;

        for (size_t i = 0; i < N; i++)
        {
            C[i] = b / A[i];
        }
        return C;
    }
    #pragma endregion

    #pragma region ArrayOps
    template <typename T, size_t N>
    KOKKOS_INLINE_FUNCTION
    T dot(const Kokkos::Array<T, N>& A, const Kokkos::Array<T, N>& B)
    {
        T dot_product = T{ 0 };

        for (size_t i = 0; i < N; i++)
        {
           dot_product += A[i] * B[i];
        }
        return dot_product;
    }
    template <typename T>
    KOKKOS_INLINE_FUNCTION
    T det2D(const Kokkos::Array<T, 2>& A, const Kokkos::Array<T, 2>& B)
    {
        T det2D = A[0] * B[1] - A[1] * B[0];

        
        return det2D;
    }
    #pragma endregion

}