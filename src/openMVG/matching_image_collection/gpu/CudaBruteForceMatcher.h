/*******************************************************************
*   CudaBruteForceMatcher.h
*   CudaBruteForceMatcher
*
*	Author: Kareem Omar
*	kareem.omar@uah.edu
*	https://github.com/komrad36
*
*	Last updated Oct 8, 2016
*******************************************************************/
//
// Fastest GPU implementation of a brute-force
// Hamming-weight matrix for 512-bit binary descriptors.
//
// Yes, that means the DIFFERENCE in popcounts is used
// for thresholding, NOT the ratio. This is the CORRECT
// approach for binary descriptors.
//
// This laboriously crafted kernel is EXTREMELY fast.
// 43 BILLION comparisons per second on a stock GTX1080,
// enough to match nearly 38,000 descriptors per frame at 30 fps (!)
//
// A key insight responsible for much of the performance of
// this insanely fast CUDA kernel is due to
// Christopher Parker (https://github.com/csp256), to whom
// I am extremely grateful.
//
// CUDA CC 3.0 or higher is required.
//
// All functionality is contained in the files CUDAK2NN.h
// and CUDAK2NN.cu. 'main.cpp' is simply a sample test harness
// with example usage and performance testing.
//

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

void cudaBruteForceMatcher(const void* const __restrict d_t, const int num_t, const cudaTextureObject_t tex_q, const int num_q, uint32_t* const __restrict, const cudaStream_t stream);
