#pragma once

#include <Kokkos_Core.hpp>


namespace NumOps
{

    #pragma region Array_Operations

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
    template <typename T, size_t N>   
    KOKKOS_INLINE_FUNCTION
    Kokkos::Array<T, N> operator+(const Kokkos::Array<T, N>& A, const T& b)
    {
        Kokkos::Array<T, N> C{};

        for (size_t i = 0; i < N; i++)
        {
            C[i] = A[i] + b;
        }
        return C;
    }
    template <typename T, size_t N>   
    KOKKOS_INLINE_FUNCTION
    Kokkos::Array<T, N> operator+(const T& b, const Kokkos::Array<T, N>& A)
    {
        Kokkos::Array<T, N> C{};

        for (size_t i = 0; i < N; i++)
        {
            C[i] = A[i] + b;
        }
        return C;
    }

    template <typename T, size_t N>   
    KOKKOS_INLINE_FUNCTION
    Kokkos::Array<T, N> operator-(const Kokkos::Array<T, N>& A, const T& b)
    {
        Kokkos::Array<T, N> C{};

        for (size_t i = 0; i < N; i++)
        {
            C[i] = A[i] - b;
        }
        return C;
    }
    template <typename T, size_t N>   
    KOKKOS_INLINE_FUNCTION
    Kokkos::Array<T, N> operator-(const T& b, const Kokkos::Array<T, N>& A)
    {
        Kokkos::Array<T, N> C{};

        for (size_t i = 0; i < N; i++)
        {
            C[i] = b - A[i];
        }
        return C;
    }
    
    template <typename T, size_t N>   
    KOKKOS_INLINE_FUNCTION
    Kokkos::Array<T, N> operator *(const Kokkos::Array<T, N>& A, const T& b)
    {
        Kokkos::Array<T, N> C{};

        for (size_t i = 0; i < N; i++)
        {
            C[i] = A[i] * b;
        }
        return C;
    }
    template <typename T, size_t N>   
    KOKKOS_INLINE_FUNCTION
    Kokkos::Array<T, N> operator *(const T& b, const Kokkos::Array<T, N>& A)
    {
        Kokkos::Array<T, N> C{};

        for (size_t i = 0; i < N; i++)
        {
            C[i] = A[i] * b;
        }
        return C;
    }

    template <typename T, size_t N>   
    KOKKOS_INLINE_FUNCTION
    Kokkos::Array<T, N> operator/(const T& b, const Kokkos::Array<T, N>& A)
    {
        Kokkos::Array<T, N> C{};

        for (size_t i = 0; i < N; i++)
        {
            C[i] = A[i] / b;
        }
        return C;
    }
    #pragma endregion
    
    #pragma region Array_functions
    template <typename T, size_t N>   
    KOKKOS_INLINE_FUNCTION
    Kokkos::Array<T, N> dot(const Kokkos::Array<T, N> A, const Kokkos::Array<T, N> B)
    {
        T dot_product{};

        for (size_t i = 0; i < N; i++)
        {
           dot_product += A[i] + B[i];
        }
        return C;
    }
    #pragma endregion

}