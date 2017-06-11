// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2016 Kareem Omar (kareem.omar@uah.edu - https://github.com/komrad36)

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

//
// Fastest GPU implementation of a brute-force
// Hamming-weight matrix for 512-bit binary descriptors.
//
// Yes, that means the DIFFERENCE in popcounts is used
// for thresholding, NOT the ratio. This is the CORRECT
// approach for binary descriptors.
//
// This laboriously crafted kernel is EXTREMELY fast.
// 63 BILLION comparisons per second on a stock GTX1080,
// enough to match nearly 48,000 descriptors per frame at 30 fps (!)
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
#define min(x) 0
#include "device_launch_parameters.h"
#define __CUDACC__
#include "device_functions.h"
#undef __CUDACC__
#endif

void cudaBruteForceMatcher(const void* const __restrict d_t, const int num_t, const cudaTextureObject_t tex_q, const int num_q, int* const __restrict d_m, const int threshold, const cudaStream_t stream);
