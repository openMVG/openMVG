// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 Pierre MOULON, Matthew Daiter

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "LatchBitMatcher.hpp"

#include <cstring>
#include <iostream>

#include "openMVG/matching_image_collection/gpu/CudaBruteForceMatcher.h"

/* Helper functions. */

#define cudaCalloc(A, B, STREAM) \
    do { \
        cudaError_t __cudaCalloc_err = cudaMalloc(A, B); \
        if (__cudaCalloc_err == cudaSuccess) cudaMemsetAsync(*A, 0, B, STREAM); \
    } while (0)

#define checkError(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true) {
   if (code != cudaSuccess) {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

namespace openMVG {
namespace matching_image_collection {
namespace gpu {

LatchBitMatcher::LatchBitMatcher(unsigned int keypoints) :
    m_maxKP(keypoints) {
  if (cudaStreamCreate(&m_stream1) == cudaErrorInvalidValue
      || cudaStreamCreate(&m_stream2) == cudaErrorInvalidValue)
    std::cerr << "Unable to create stream" << std::endl;

  const size_t sizeD = m_maxKP * 64; // D for Descriptor
  const size_t sizeMatches = m_maxKP * sizeof(uint32_t); // M for Matches

  cudaCalloc((void**) &m_dQuery, sizeD, m_stream1);
  cudaCalloc((void**) &m_dTraining, sizeD, m_stream2);

  cudaCalloc((void**) &m_dMatches1, sizeMatches, m_stream1);
  cudaCalloc((void**) &m_dMatches2, sizeMatches, m_stream2);

  cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);
  cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeEightByte);

  struct cudaResourceDesc resDescQuery;
  memset(&resDescQuery, 0, sizeof(resDescQuery));
  resDescQuery.resType = cudaResourceTypeLinear;
  resDescQuery.res.linear.devPtr = m_dQuery;
  resDescQuery.res.linear.desc.f = cudaChannelFormatKindUnsigned;
  resDescQuery.res.linear.desc.x = 32;
  resDescQuery.res.linear.desc.y = 32;
  resDescQuery.res.linear.sizeInBytes = sizeD;
  
  struct cudaTextureDesc texDescQuery;
  memset(&texDescQuery, 0, sizeof(texDescQuery));
  texDescQuery.addressMode[0] = cudaAddressModeBorder;
  texDescQuery.addressMode[1] = cudaAddressModeBorder;
  texDescQuery.filterMode = cudaFilterModePoint;
  texDescQuery.readMode = cudaReadModeElementType;
  texDescQuery.normalizedCoords = 0;
  
  struct cudaResourceDesc resDescTraining;
  memset(&resDescTraining, 0, sizeof(resDescTraining));
  resDescTraining.resType = cudaResourceTypeLinear;
  resDescTraining.res.linear.devPtr = m_dTraining;
  resDescTraining.res.linear.desc.f = cudaChannelFormatKindUnsigned;
  resDescTraining.res.linear.desc.x = 32;
  resDescTraining.res.linear.desc.y = 32;
  resDescTraining.res.linear.sizeInBytes = sizeD;
  
  struct cudaTextureDesc texDescTraining;
  memset(&texDescTraining, 0, sizeof(texDescTraining));
  texDescTraining.addressMode[0] = cudaAddressModeBorder;
  texDescTraining.addressMode[1] = cudaAddressModeBorder;
  texDescTraining.filterMode = cudaFilterModePoint;
  texDescTraining.readMode = cudaReadModeElementType;
  texDescTraining.normalizedCoords = 0;

  m_texQuery = 0;
  m_texTraining = 0;

  checkError(cudaCreateTextureObject(&m_texQuery, &resDescQuery, &texDescQuery, nullptr));
  checkError(cudaCreateTextureObject(&m_texTraining, &resDescTraining, &texDescTraining, nullptr));
}

void LatchBitMatcher::match(void* h_descriptorsQuery, void* h_descriptorsTraining, int numKPQuery, int numKPTraining) {
  m_numKPQuery = numKPQuery;
  m_numKPTraining = numKPTraining;

  const size_t sizeDQuery = numKPQuery * 64; // D1 for descriptor1
  const size_t sizeDTraining = numKPTraining * 64; // D2 for descriptor2
  const size_t sizeMatches = m_maxKP * sizeof(int); // M for Matches

  cudaMemsetAsync(m_dQuery, 0, sizeDQuery, m_stream1);
  cudaMemsetAsync(m_dTraining, 0, sizeDTraining, m_stream2);

  cudaMemsetAsync(m_dMatches1, 0, sizeMatches, m_stream1);
  cudaMemsetAsync(m_dMatches2, 0, sizeMatches, m_stream2);

  cudaMemcpyAsync(m_dQuery, h_descriptorsQuery, sizeDQuery, cudaMemcpyHostToDevice, m_stream1);
  cudaMemcpyAsync(m_dTraining, h_descriptorsTraining, sizeDTraining, cudaMemcpyHostToDevice, m_stream2);
 
  cudaStreamSynchronize(m_stream1);
  cudaStreamSynchronize(m_stream2);

  cudaBruteForceMatcher(m_dTraining, m_numKPTraining, m_texQuery, m_numKPQuery, m_dMatches1, m_stream1);
  cudaBruteForceMatcher(m_dQuery, m_numKPQuery, m_texTraining, m_numKPTraining, m_dMatches2, m_stream2);
}

openMVG::matching::IndMatches LatchBitMatcher::retrieveMatches(float ratio) {
  int h_Matches1[m_maxKP];
  int h_Matches2[m_maxKP];
  cudaMemcpyAsync(h_Matches1, m_dMatches1, m_maxKP * sizeof(int), cudaMemcpyDeviceToHost, m_stream1);
  cudaMemcpyAsync(h_Matches2, m_dMatches2, m_maxKP * sizeof(int), cudaMemcpyDeviceToHost, m_stream2);

  cudaStreamSynchronize(m_stream1);
  cudaStreamSynchronize(m_stream2);

  openMVG::matching::IndMatches matches;
  const uint32_t threshold = std::ceil(ratio * 16);

#ifdef OPENMVG_USE_OMP
  #pragma omp parallel for schedule(dynamic)
#endif
  for (size_t i = 0; i < m_numKPQuery; i++) {
   if (h_Matches1[i] >= ratio && h_Matches1[i] < m_numKPTraining && h_Matches2[h_Matches1[i]] == i)
      matches.push_back(openMVG::matching::IndMatch(i, h_Matches1[i]));
  }
  return matches;
}

LatchBitMatcher::~LatchBitMatcher() {
  cudaFree(m_dQuery);
  cudaFree(m_dTraining);
  cudaFree(m_dMatches1);
  cudaFree(m_dMatches2);
  cudaStreamDestroy(m_stream1);
  cudaStreamDestroy(m_stream2);
}

} // namespace gpu
} // namespace matching_image_collection
} // namespace openMVG
