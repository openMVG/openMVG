#include "LatchBitMatcher.hpp"

#include <iostream>

#include "bitMatcher.h"
#include "params.hpp"

/* Helper functions. */

#define cudaCalloc(A, B) \
    do { \
        cudaError_t __cudaCalloc_err = cudaMalloc(A, B); \
        if (__cudaCalloc_err == cudaSuccess) cudaMemset(*A, 0, B); \
    } while (0)

#define checkError(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true) {
   if (code != cudaSuccess) {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

LatchBitMatcher::LatchBitMatcher() :
    m_maxKP(512 * NUM_SM),
    m_matchThreshold(MATCH_THRESHOLD) {
    if (cudaStreamCreate(&m_stream1) == cudaErrorInvalidValue
        || cudaStreamCreate(&m_stream2) == cudaErrorInvalidValue) {
        std::cerr << "Unable to create streams" << std::endl;
    }
    if (cudaEventCreate(&m_finished) != cudaSuccess) {
        std::cerr << "Unable to create event" << std::endl;
    }
    size_t sizeD = m_maxKP * (2048 / 32) * sizeof(unsigned int); // D for Descriptor
	size_t sizeM = m_maxKP * sizeof(int); // M for Matches

    cudaCalloc((void**) &m_dD1, sizeD);
    cudaCalloc((void**) &m_dD2, sizeD);
	cudaCalloc((void**) &m_dM1, sizeM);
	cudaCalloc((void**) &m_dM2, sizeM);
}

std::vector<LatchBitMatcherMatches> LatchBitMatcher::match(unsigned int* h_descriptors1, unsigned int* h_descriptors2, int numKP0, int numKP1) {
    size_t sizeD = m_maxKP * (2048 / 32) * sizeof(unsigned int); // D for descriptor

    cudaMemcpyAsync(m_dD1, h_descriptors1, sizeD, cudaMemcpyHostToDevice, m_stream1);
    cudaMemcpyAsync(m_dD2, h_descriptors2, sizeD, cudaMemcpyHostToDevice, m_stream2);

    bitMatcher(m_dD1, m_dD2, numKP0, numKP1, m_maxKP, m_dM1, m_matchThreshold, m_stream1, m_finished);
    bitMatcher(m_dD2, m_dD1, numKP1, numKP0, m_maxKP, m_dM2, m_matchThreshold, m_stream2, m_finished);

    cudaStreamSynchronize(m_stream1);
    cudaStreamSynchronize(m_stream2);

    int h_M1[m_maxKP];
    int h_M2[m_maxKP];
    getMatches(m_maxKP, h_M1, m_dM1);
    getMatches(m_maxKP, h_M2, m_dM2);

    std::vector<LatchBitMatcherMatches> matches;

    for (size_t i = 0; i < numKP0; i++) {
        if (h_M1[i] >= 0 && h_M1[i] < numKP1 && h_M2[h_M1[i]] == i) {
            // Add matches
            matches.push_back(LatchBitMatcherMatches(i, h_M1[i], 0));
        }
    }

    return matches;
}

LatchBitMatcher::~LatchBitMatcher() {
    cudaEventDestroy(m_finished);
    cudaStreamDestroy(m_stream1);
    cudaStreamDestroy(m_stream2);
}
