#include "LatchBitMatcher.hpp"

#include <iostream>

#include "bitMatcher.h"
#include "params.hpp"

#define DESCRIPTOR_SIZE 2048

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
    size_t sizeD = m_maxKP * (DESCRIPTOR_SIZE / 32) * sizeof(unsigned int); // D for Descriptor
	size_t sizeM = m_maxKP * sizeof(int); // M for Matches

    cudaCalloc((void**) &m_dD1, sizeD);
    cudaCalloc((void**) &m_dD2, sizeD);

	cudaCalloc((void**) &m_dM1, sizeM);
	cudaCalloc((void**) &m_dM2, sizeM);
}

std::vector<LatchBitMatcherMatch> LatchBitMatcher::match(unsigned int* h_descriptors1, unsigned int* h_descriptors2, int numKP0, int numKP1) {
    size_t sizeD1 = numKP0 * (DESCRIPTOR_SIZE / 32) * sizeof(unsigned int); // D1 for descriptor1
    size_t sizeD2 = numKP1 * (DESCRIPTOR_SIZE / 32) * sizeof(unsigned int); // D2 for descriptor2
	size_t sizeM = m_maxKP * sizeof(int); // M for Matches

	cudaMemsetAsync(m_dD1, 0, sizeD1, m_stream1);
	cudaMemsetAsync(m_dD2, 0, sizeD2, m_stream2);

	cudaMemsetAsync(m_dM1, 0, sizeM, m_stream1);
	cudaMemsetAsync(m_dM2, 0, sizeM, m_stream2);
/*
	cudaStreamSynchronize(m_stream1);
	cudaStreamSynchronize(m_stream2);
	cudaDeviceSynchronize();
*/
    cudaMemcpyAsync(m_dD1, h_descriptors1, sizeD1, cudaMemcpyHostToDevice, m_stream1);
    cudaMemcpyAsync(m_dD2, h_descriptors2, sizeD2, cudaMemcpyHostToDevice, m_stream2);
/*    
	cudaStreamSynchronize(m_stream1);
	cudaStreamSynchronize(m_stream2);
	cudaDeviceSynchronize();
*/
    bitMatch(m_dD1, m_dD2, numKP0, numKP1, m_maxKP, m_dM1, m_matchThreshold, m_stream1);
    bitMatch(m_dD2, m_dD1, numKP1, numKP0, m_maxKP, m_dM2, m_matchThreshold, m_stream2);

//    int *h_M1 = new int[m_maxKP];
//    int *h_M2 = new int[m_maxKP];
	int h_M1[m_maxKP];
	int h_M2[m_maxKP];

	cudaDeviceSynchronize();
    getMatches(m_maxKP, h_M1, m_dM1, m_stream1);
    getMatches(m_maxKP, h_M2, m_dM2, m_stream2);

    std::vector<LatchBitMatcherMatch> matches;

	cudaStreamSynchronize(m_stream1);
	cudaStreamSynchronize(m_stream2);

	cudaDeviceSynchronize();

    for (size_t i = 0; i < numKP0; i++) {
        if (h_M1[i] >= 0 && h_M1[i] < numKP1 && h_M2[h_M1[i]] == i) {
            // Add matches
            matches.push_back(LatchBitMatcherMatch(i, h_M1[i], 0));
        }
    }

//	delete [] h_M1;
//	delete [] h_M2;

    return matches;
}

LatchBitMatcher::~LatchBitMatcher() {
	cudaFree(m_dD1);
	cudaFree(m_dD2);
	cudaFree(m_dM1);
	cudaFree(m_dM2);
    cudaStreamDestroy(m_stream1);
    cudaStreamDestroy(m_stream2);
}
