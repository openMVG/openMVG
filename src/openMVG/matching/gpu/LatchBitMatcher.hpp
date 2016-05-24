#pragma once

#include <vector>

#include <cuda.h>
#include <cuda_runtime.h>

struct LatchBitMatcherMatch {
    LatchBitMatcherMatch() :
        queryIdx(0),
        trainIdx(0),
        distance(0.0) {
    }

    LatchBitMatcherMatch(int _queryIdx, int _trainIdx, int _distance) :
        queryIdx(_queryIdx),
        trainIdx(_trainIdx),
        distance(_distance) {

    }

    int queryIdx;
    int trainIdx;
    float distance;
};

class LatchBitMatcher {
    public:
        LatchBitMatcher();
        std::vector<LatchBitMatcherMatch> match(unsigned int*, unsigned int*, int, int);
        ~LatchBitMatcher();
    private:
        const int m_maxKP;
        const int m_matchThreshold;

        unsigned int* m_dD1;
        unsigned int* m_dD2;

        int* m_dM1;
        int* m_dM2;

        // CUDA stuff
        cudaStream_t m_stream1;
        cudaStream_t m_stream2;
        cudaEvent_t m_finished;
};
