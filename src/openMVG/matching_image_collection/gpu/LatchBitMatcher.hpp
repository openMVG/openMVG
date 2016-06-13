#pragma once

#include <vector>

#include <cuda.h>
#include <cuda_runtime.h>

#include "openMVG/features/latch/LatchBitMatcherMatch.hpp"

class LatchBitMatcher {
    public:
        LatchBitMatcher();
       	void match(unsigned int*, unsigned int*, int, int);
		std::vector<LatchBitMatcherMatch> retrieveMatches();
        ~LatchBitMatcher();
    private:
        const int m_maxKP;
        const int m_matchThreshold;

        unsigned int* m_dD1;
        unsigned int* m_dD2;

		unsigned int m_numKP0;
		unsigned int m_numKP1;

        int* m_dM1;
        int* m_dM2;

        // CUDA stuff
        cudaStream_t m_stream1;
        cudaStream_t m_stream2;
};
