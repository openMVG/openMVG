/*******************************************************************
*   CLATCH.cu
*   KORAL
*
*	Author: Kareem Omar
*	kareem.omar@uah.edu
*	https://github.com/komrad36
*
*	Last updated Dec 27, 2016
*******************************************************************/
//
// ## Summary ##
// KORAL is a novel, extremely fast, highly accurate, scale- and
// rotation-invariant, CPU-GPU cooperative detector-descriptor.
//
// Detection is based on the author's custom multi-scale KFAST corner
// detector, with rapid bilinear interpolation performed by the GPU
// asynchronously while the CPU works on KFAST.
//
// ## Usage ##
// Basic use of KORAL is extremely easy, although, of course, for a
// larger high-performance pipeline, users will benefit from
// calling KORAL functions directly and modifying it to suit their needs.
//
// To detect and describe, simply #include "KORAL.h" and
// then do:
//
// 	    KORAL koral(scale_factor, scale_levels);
//      koral.go(image, width, height, KFAST_threshold);
//
// where scale_factor is the factor by which each scale leve
// is reduced from the previous, scale_levels is the total
// number of such scale levels used, image is a pointer to
// uint8_t (grayscale) image data, and KFAST_threshold
// is the threshold supplied to the KFAST feature detector.
//
// After this call, keypoints are avaiable in a vector at 
// koral.kps, while descriptors are available at
// koral.desc.
//
// Portions of KORAL require SSE, AVX, AVX2, and CUDA.
// The author is working on reduced-performance versions
// with lesser requirements, but as the intent of this work
// is primarily novel performance capability, modern
// hardware and this full version are highly recommended.
//
// Description is performed by the GPU using the novel CLATCH
// (CUDA LATCH) binary descriptor kernel.
//
// Rotation invariance is provided by a novel vectorized
// SSE angle weight detector.
//
// All components have been written and carefully tuned by the author
// for maximum performance and have no external dependencies. Some have
// been modified for integration into KORAL,
// but the original standalone projects are all availble on
// the author's GitHub (https://github.com/komrad36).
//
// These individual components are:
// -KFAST        (https://github.com/komrad36/KFAST)
// -CUDALERP     (https://github.com/komrad36/CUDALERP)
// -FeatureAngle (https://github.com/komrad36/FeatureAngle)
// -CLATCH       (https://github.com/komrad36/CLATCH)
//
// In addition, the natural next step of matching descriptors
// is available in the author's currently separate
// project, CUDAK2NN (https://github.com/komrad36/CUDAK2NN).
//
// A key insight responsible for much of the performance of
// this insanely fast system is due to Christopher Parker
// (https://github.com/csp256), to whom I am extremely grateful.
// 
// The file 'main.cpp' is a simple test driver
// illustrating example usage.It requires OpenCV
// for image read and keypoint display.KORAL itself,
// however, does not require OpenCV or any other
// external dependencies.
//
// Note that KORAL is a work in progress.
// Suggestions and improvements are welcomed.
//
// ## License ##
// The FAST detector was created by Edward Rosten and Tom Drummond
// as described in the 2006 paper by Rosten and Drummond:
// "Machine learning for high-speed corner detection"
//         Edward Rosten and Tom Drummond
// https://www.edwardrosten.com/work/rosten_2006_machine.pdf
//
// The FAST detector is BSD licensed:
// 
// Copyright(c) 2006, 2008, 2009, 2010 Edward Rosten
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met :
// 
// 
// *Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
// 
// *Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and / or other materials provided with the distribution.
// 
// *Neither the name of the University of Cambridge nor the names of
// its contributors may be used to endorse or promote products derived
// from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED.IN NO EVENT SHALL THE COPYRIGHT OWNER OR
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES(INCLUDING, BUT NOT LIMITED TO,
// 	PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// 	PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// 	LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT(INCLUDING
// 		NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// 	SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
//
//
//
// KORAL is licensed under the MIT License : https://opensource.org/licenses/mit-license.php
// 
// Copyright(c) 2016 Kareem Omar, Christopher Parker
// 
// Permission is hereby granted, free of charge,
// to any person obtaining a copy of this software and associated documentation
// files(the "Software"), to deal in the Software without restriction, including
// without limitation the rights to use, copy, modify, merge, publish, distribute,
// sublicense, and / or sell copies of the Software, and to permit persons to whom
// the Software is furnished to do so, subject to the following conditions :
// 
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
// INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
// PARTICULAR PURPOSE AND NONINFRINGEMENT.IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
// HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
// SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//
// 
// Note again that KORAL is a work in progress.
// Suggestions and improvements are welcomed.
//

#include "CLATCH.h"

__global__ void
#ifndef __INTELLISENSE__
__launch_bounds__(512, 4)
#endif
CLATCH_kernel(const cudaTextureObject_t d_all_tex[8], const cudaTextureObject_t d_triplets, const Keypoint* const __restrict d_kps, uint32_t* const __restrict__ d_desc) {
	volatile __shared__ uint8_t s_ROI[4608];
	const Keypoint pt = d_kps[blockIdx.x];
	const cudaTextureObject_t d_img_tex = d_all_tex[pt.scale];
	const float s = sin(pt.angle), c = cos(pt.angle);
	for (int32_t i = 0; i <= 48; i += 16) {
		for (int32_t k = 0; k <= 32; k += 32) {
			const float x_offset = static_cast<float>(static_cast<int>(threadIdx.x) + k - 32);
			const float y_offset = static_cast<float>(static_cast<int>(threadIdx.y) + i - 32);
			s_ROI[(threadIdx.y + i) * 72 + threadIdx.x + k] = tex2D<uint8_t>(d_img_tex, static_cast<int>((pt.x + (x_offset*c - y_offset*s)) + 0.5f), static_cast<int>((pt.y + (x_offset*s + y_offset*c)) + 0.5f));
		}
	}
	uint32_t ROI_base = 144 * (threadIdx.x & 3) + (threadIdx.x >> 2), triplet_base = threadIdx.y << 5, desc = 0;
	__syncthreads();
	for (int32_t i = 0; i < 4; ++i, triplet_base += 8) {
		int32_t accum[8];
		for (uint32_t j = 0; j < 8; ++j) {
			const ushort4 t = tex1D<ushort4>(d_triplets, triplet_base + j);
			const int32_t b1 = s_ROI[ROI_base + t.y],      b2 = s_ROI[ROI_base + t.y + 72]     ;
			const int32_t a1 = s_ROI[ROI_base + t.x] - b1, a2 = s_ROI[ROI_base + t.x + 72] - b2;
			const int32_t c1 = s_ROI[ROI_base + t.z] - b1, c2 = s_ROI[ROI_base + t.z + 72] - b2;
			accum[j] = a1 * a1 - c1 * c1 + a2 * a2 - c2 * c2;
		}
		for (int32_t k = 1; k <= 4; k <<= 1) {
			for (int32_t s = 0; s < 8; s += k) accum[s] += __shfl_xor(accum[s], k);
			if (threadIdx.x & k) for (int32_t s = 0; s < 8; s += k << 1) accum[s] = accum[s + k];
		}
		accum[0] += __shfl_xor(accum[0], 8);
		desc |= (accum[0] + __shfl_xor(accum[0], 16) < 0) << ((i << 3) + (threadIdx.x & 7));
	}
	for (int32_t s = 1; s <= 4; s <<= 1) desc |= __shfl_xor(desc, s);
	if (threadIdx.x == 0) d_desc[(blockIdx.x << 4) + threadIdx.y] = desc;
}

void CLATCH(cudaTextureObject_t d_all_tex[8], const cudaTextureObject_t d_triplets, const Keypoint* const __restrict d_kps, const int num_kps, uint64_t* const __restrict d_desc) {
	CLATCH_kernel<<<num_kps, { 32, 16 } >>>(d_all_tex, d_triplets, d_kps, reinterpret_cast<uint32_t*>(d_desc));
}
