#pragma once

#include <cstdint>
#include <vector>

#include <cuda.h>
#include <cuda_runtime.h>

#include "openMVG/features/latch/LatchBitMatcherMatch.hpp"

class LatchBitMatcher {
  public:
    LatchBitMatcher();
    void match(void*, void*, int, int);
	  std::vector<LatchBitMatcherMatch> retrieveMatches();
    ~LatchBitMatcher();
	private:
    const int m_maxKP;
    const int m_matchThreshold;

    void* m_dQuery;
    void* m_dTraining;

		unsigned int m_numKPQuery;
		unsigned int m_numKPTraining;

    int* m_dMatches1;
    int* m_dMatches2;

    cudaTextureObject_t m_texQuery;
    cudaTextureObject_t m_texTraining;
    // CUDA stuff
    cudaStream_t m_stream1;
    cudaStream_t m_stream2;
};
