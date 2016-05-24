#ifndef LATCH_BIT_MATCHER_H
#define LATCH_BIT_MATCHER_H

#include <vector>

#include <cuda.h>
#include <cuda_runtime.h>

struct LatchBitMatcherMatches {
    LatchBitMatcherMatches() :
        queryIdx(0),
        trainIdx(0),
        accuracy(0) {
    }

    LatchBitMatcherMatches(int _queryIdx, int _trainIdx, int _accuracy) :
        queryIdx(_queryIdx),
        trainIdx(_trainIdx),
        accuracy(_accuracy) {

    }

    int queryIdx;
    int trainIdx;
    int accuracy;
};

class LatchBitMatcher {
    public:
        LatchBitMatcher();
        std::vector<LatchBitMatcherMatches> match(unsigned int*, unsigned int*, int, int);
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

#endif
