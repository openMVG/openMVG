#pragma once

#include <cstdint>

#include "cuda_runtime.h"

#ifdef __INTELLISENSE__
#define asm(x)
#include "device_launch_parameters.h"
#define __CUDACC__
#include "device_functions.h"
#undef __CUDACC__
#endif

void cudaBruteForceMatcher(const void* const __restrict d_t, const int num_t, const cudaTextureObject_t tex_q, const int num_q, int* d_m, const int threshold, const cudaStream_t stream);
