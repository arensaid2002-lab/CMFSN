#pragma once

// CMake's CUDA sources need the same preincluded runtime and compiler
// version information in Visual Studio's editor as they receive from NVCC.
// This file has no effect during compilation.
#ifdef __INTELLISENSE__
#include <cuda.h>
#ifndef __CUDACC__
#define __CUDACC__ 1
#endif
#ifndef __NVCC__
#define __NVCC__ 1
#endif
#ifndef __CUDACC_VER_MAJOR__
#define __CUDACC_VER_MAJOR__ (CUDA_VERSION / 1000)
#endif
#ifndef __CUDACC_VER_MINOR__
#define __CUDACC_VER_MINOR__ ((CUDA_VERSION % 1000) / 10)
#endif
#include <cuda_runtime.h>
// The editor does not implement this CUDA parameter attribute. The real
// CUDA compiler still sees the original definition from its runtime.
#ifdef __grid_constant__
#undef __grid_constant__
#endif
#define __grid_constant__
#endif
