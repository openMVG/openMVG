#include <stdio.h>
#include <iostream>
#include <errno.h>

#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <device_launch_parameters.h>

#include "bitMatcher.h"

using namespace std;

// Number of values each thread in a warp gets per vector.
#define chunksPerVector (2)
#define vectorsPerWarp (16)
// Vectors per group is used to increase ILP. it must divide vectorsPerWarp. This implementation is specialized for vectorsPerGroup==8.
#define vectorsPerGroup (8)
#define warpsPerBlock (32)
// The total number of int32's needed to store a vector. We should drop this down to 16 for an optimized implementation for canonical LATCH.
#define vectorDimension (16)
#define _warpSize (32)
#define cacheSize (128)
#define halfCacheSize (64)

#define checkError(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true) {
   if (code != cudaSuccess) {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

// Launch as 32x32
__global__ void __launch_bounds__(1024, 1)
                bitMatch(   const unsigned int *g_query,
                            const unsigned int *g_training,
                            int *g_match,
                            const int trainingSize,
                            const int threshold) {
    // Load query vectors
    register unsigned int query[vectorsPerWarp][chunksPerVector];

    volatile __shared__ unsigned int s_training[cacheSize][chunksPerVector][_warpSize]; // We have enough room to load extra query vectors from shared memory...
    {
        int offset = threadIdx.x;
        offset += blockIdx.x * vectorDimension * warpsPerBlock * vectorsPerWarp;
        offset +=              vectorDimension *  threadIdx.y  * vectorsPerWarp;

        #pragma unroll
        for (int i=0; i<vectorsPerWarp; i++) {
            #pragma unroll
            for (int j=0; j<chunksPerVector; j++, offset += warpSize) {
                query[i][j] = g_query[offset];
            }
        }
    }

    // Load the first training vectors.
    int trainingOffset = threadIdx.y * vectorDimension;
    if (threadIdx.y < halfCacheSize) {
        for (int i=0; i<chunksPerVector; i++, trainingOffset += warpSize) {
            s_training[threadIdx.y][i][threadIdx.x] = g_training[trainingOffset + threadIdx.x];
        }
    }
    __threadfence_block();

    register int bestIndex = -1;
    register int best       =  9999999;
    register int secondBest = 99999999;
    #pragma unroll 4
    for (int t=0; t < trainingSize; t+= cacheSize) {
        // Synchronize halfway through using shared memory...
        // So you can freely write to the other half.
        #pragma unroll
        for (int half=0; half < 2; half++) { // Half will be 0 when you should be working with top half, and loading into bottom half.
            register unsigned int prefetch = 0.0f;
            #pragma unroll
            for (int st=0; st < halfCacheSize; st++) { // Every iteration of this loop must load a single training vector into shared memory.
                {
                    // Stream a new pair of training vectors to registers at start of every even loop (and write them to shared memory at end of every odd loop)
                    if (st % 2 == 0) {
                        if (threadIdx.y < 2*chunksPerVector) {
                            const int index = (t + (half+1)*halfCacheSize + st)*vectorDimension + threadIdx.y*_warpSize + threadIdx.x;
                            if (index < trainingSize*vectorDimension) {
                                prefetch = g_training[index];
                            }
                        }
                    }
                }
                {
                    // This is the offset into our shared memory cache of training vectors.
                    const register int trainingOffset = half*halfCacheSize + st;
                    register unsigned int train[chunksPerVector];

                    // Load training vector into registers.
                    #pragma unroll
                    for (int chunk = 0; chunk < chunksPerVector; chunk++) {
                        train[chunk] = s_training[trainingOffset][chunk][threadIdx.x];
                    }

                    // The compiler throws a hissy fit if you try to make dist an array, and tosses everything into local memory.
                    register int dist0, dist1, dist2, dist3, dist4, dist5, dist6, dist7;
                    // Also, the compiler does not like this being in a (fully unrolled) loop... drama queen.
                    dist0 = __popc(query[0][0] ^ train[0]);// + __popc(query[0][1] ^ train[1]);
                    dist1 = __popc(query[1][0] ^ train[0]);// + __popc(query[1][1] ^ train[1]);
                    dist2 = __popc(query[2][0] ^ train[0]);// + __popc(query[2][1] ^ train[1]);
                    dist3 = __popc(query[3][0] ^ train[0]);// + __popc(query[3][1] ^ train[1]);
                    dist4 = __popc(query[4][0] ^ train[0]);// + __popc(query[4][1] ^ train[1]);
                    dist5 = __popc(query[5][0] ^ train[0]);// + __popc(query[5][1] ^ train[1]);
                    dist6 = __popc(query[6][0] ^ train[0]);// + __popc(query[6][1] ^ train[1]);
                    dist7 = __popc(query[7][0] ^ train[0]);// + __popc(query[7][1] ^ train[1]);
                    dist0 |= (__popc(query[ 8][0] ^ train[0]) /*+ __popc(query[ 8][1] ^ train[1])*/)<<16;
                    dist1 |= (__popc(query[ 9][0] ^ train[0]) /*+ __popc(query[ 9][1] ^ train[1])*/)<<16;
                    dist2 |= (__popc(query[10][0] ^ train[0]) /*+ __popc(query[10][1] ^ train[1])*/)<<16;
                    dist3 |= (__popc(query[11][0] ^ train[0]) /*+ __popc(query[11][1] ^ train[1])*/)<<16;
                    dist4 |= (__popc(query[12][0] ^ train[0]) /*+ __popc(query[12][1] ^ train[1])*/)<<16;
                    dist5 |= (__popc(query[13][0] ^ train[0]) /*+ __popc(query[13][1] ^ train[1])*/)<<16;
                    dist6 |= (__popc(query[14][0] ^ train[0]) /*+ __popc(query[14][1] ^ train[1])*/)<<16;
                    dist7 |= (__popc(query[15][0] ^ train[0]) /*+ __popc(query[15][1] ^ train[1])*/)<<16;

                    dist0 += __shfl_xor(dist0,   1);
                    dist1 += __shfl_xor(dist1,   1);
                    if (threadIdx.x & 1) dist0 = dist1;
                    dist2 += __shfl_xor(dist2,   1);
                    dist3 += __shfl_xor(dist3,   1);
                    if (threadIdx.x & 1) dist2 = dist3;
                    dist4 += __shfl_xor(dist4,   1);
                    dist5 += __shfl_xor(dist5,   1);
                    if (threadIdx.x & 1) dist4 = dist5;
                    dist6 += __shfl_xor(dist6,   1);
                    dist7 += __shfl_xor(dist7,   1);
                    if (threadIdx.x & 1) dist6 = dist7;
                    dist0 += __shfl_xor(dist0,   2);
                    dist2 += __shfl_xor(dist2,   2);
                    if (threadIdx.x & 2) dist0 = dist2;
                    dist4 += __shfl_xor(dist4,   2);
                    dist6 += __shfl_xor(dist6,   2);
                    if (threadIdx.x & 2) dist4 = dist6;
                    dist0 += __shfl_xor(dist0,   4);
                    dist4 += __shfl_xor(dist4,   4);
                    if (threadIdx.x & 4) dist0 = dist4;
                    dist0 += __shfl_xor(dist0,   8);
                    dist0 += __shfl_xor(dist0,  16);
                    if (threadIdx.x < 8) dist0 &= 2047;
                    else dist0 >>= 16;

                    if (dist0 < secondBest) {
                        if (dist0 < best) {
                            secondBest = best;
                            best = dist0;
                            bestIndex = t + trainingOffset;
                        } else {
                            secondBest = dist0;
                        }
                    }
                }
                { // Write new training vectors prefetched into registers to shared memory cache at end of every even loop.
                    if (st % 2 == 1) {
                        if (threadIdx.y < chunksPerVector) { // We can load identically for each chunk, but not so for write to shared memory differently.
                            s_training[(half^1)*halfCacheSize + (st-1)    ][threadIdx.y                  ][threadIdx.x] = prefetch;
                        } else if (threadIdx.y < 2*chunksPerVector) {
                            s_training[(half^1)*halfCacheSize + (st-1) + 1][threadIdx.y - chunksPerVector][threadIdx.x] = prefetch;
                        }
                    }
                }
            }
            __syncthreads();
        }
    }
    if (threadIdx.x < vectorsPerWarp) {
        if (secondBest - best < threshold) {
            bestIndex = -1; // Failed hard threshold test.
        }
        // We can trash what is in shared memory now... it is called s_training, but here it is just scratch space.
        // I guess I should use a union for this?
        const register int packing = _warpSize / vectorsPerWarp; // NOTE: This assumes vectorsPerWarp divides _warpSize. If it doesnt, you'll have to handle this differently.
        s_training[0][threadIdx.y / packing][(threadIdx.y%packing)*vectorsPerWarp + threadIdx.x] = bestIndex;
    }
    __threadfence_block();
    if (threadIdx.y < vectorsPerWarp) {
        g_match[blockIdx.x*vectorsPerWarp*warpsPerBlock + threadIdx.y*warpsPerBlock + threadIdx.x] = s_training[0][threadIdx.y][threadIdx.x];
    }
}


void bitMatcher(unsigned int* d_Q, unsigned int* d_T, int keypointsQ, int keypointsT, int maxKP, int* d_M, const int threshold, cudaStream_t stream, cudaEvent_t event) {
    dim3 threadsPerBlock(_warpSize, warpsPerBlock);
    const int neededBlocks = (keypointsQ + (vectorsPerWarp * warpsPerBlock) - 1) / (vectorsPerWarp * warpsPerBlock); // This is the "round up integer division" pattern
    dim3 blocksPerGrid(neededBlocks, 1, 1);

    cudaStreamWaitEvent(stream, event, 0);
    bitMatch<<<blocksPerGrid, threadsPerBlock, 0, stream>>>(d_Q, d_T, d_M, keypointsT, threshold);
}

void getMatches(int maxKP, int* h_M, int* d_M) {
    size_t sizeM = maxKP * sizeof(int);
    cudaMemcpyAsync(h_M, d_M, sizeM, cudaMemcpyDeviceToHost);
};
