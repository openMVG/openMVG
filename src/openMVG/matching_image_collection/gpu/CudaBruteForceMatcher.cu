#include "CudaBruteForceMatcher.h"

__global__ void CUDABRUTE_kernel(const cudaTextureObject_t tex_q, const int num_q, const uint64_t* const __restrict__ g_training, const int num_t, int* const __restrict__ g_match, const int threshold) {
	uint64_t q[8];
	for (int i = 0, offset = ((threadIdx.x & 24) << 3) + (threadIdx.x & 7) + (blockIdx.x << 11) + (threadIdx.y << 8); i < 8; ++i, offset += 8) {
		const uint2 buf = tex1Dfetch<uint2>(tex_q, offset);
		asm("mov.b64 %0, {%1,%2};" : "=l"(q[i]) : "r"(buf.x), "r"(buf.y));
	}
	int best_i, best_v = 100000, second_v = 200000;
#pragma unroll 7
	for (int t = 0; t < num_t; ++t) {
		const uint64_t train = g_training[(t << 3) + (threadIdx.x & 7)];
		uint32_t dist[4];
		for (int i = 0; i < 4; ++i) dist[i] = __byte_perm(__popcll(q[i] ^ train), __popcll(q[i + 4] ^ train), 0x5410);
#pragma unroll
		for (int i = 0; i < 4; ++i) dist[i] += __shfl_xor(dist[i], 1);
		if (threadIdx.x & 1) dist[0] = dist[1];
		if (threadIdx.x & 1) dist[2] = dist[3];
		dist[0] += __shfl_xor(dist[0], 2);
		dist[2] += __shfl_xor(dist[2], 2);
		if (threadIdx.x & 2) dist[0] = dist[2];
		dist[0] = __byte_perm(dist[0] + __shfl_xor(dist[0], 4), 0, threadIdx.x & 4 ? 0x5432 : 0x5410);
		if (dist[0] < second_v) second_v = dist[0];
		if (dist[0] < best_v) {
			second_v = best_v;
			best_i = t;
			best_v = dist[0];
		}
	}
	const int idx = (blockIdx.x << 8) + (threadIdx.y << 5) + threadIdx.x;
	if (idx < num_q) g_match[idx] = second_v - best_v > threshold ? best_i : -1;
}

void cudaBruteForceMatcher(const void* const __restrict d_t, const int num_t, const cudaTextureObject_t tex_q, const int num_q, int* d_m, const int threshold, const cudaStream_t stream) {
	CUDABRUTE_kernel<<<((num_q - 1)>>8) + 1, { 32, 8 }, 0, stream>>>(tex_q, num_q, reinterpret_cast<const uint64_t*>(d_t), num_t, d_m, threshold);
}
