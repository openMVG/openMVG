/*******************************************************************
*   KFAST.h
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

#pragma once

#include <algorithm>
#include <cstdint>
#include <future>
#include <immintrin.h>
#include <vector>

#include "Keypoint.h"

// Yes, this function MUST be inlined.
// Even if your compiler thinks otherwise.
// 2000 -> 2600 microseconds without forced inlining.
template<const bool full, const bool nonmax_suppression>
#ifdef _MSC_VER
__forceinline
#else
inline __attribute__((always_inline))
#endif
void processCols(int32_t& num_corners, const uint8_t* __restrict & ptr, int32_t& j,
	const int32_t* const __restrict offsets, const __m256i& ushft, const __m256i& t, const int32_t cols,
	const __m256i& consec, int32_t* const __restrict corners, uint8_t* const __restrict cur,
	std::vector<Keypoint>& keypoints, const int32_t i, const int32_t start_row) {
	// ppt is an integer vector that now holds 32 of point p
	__m256i ppt = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(ptr));

	// we subtract (and clamp) the threshold value from all 32 pixels
	// pmt represents p - t
	__m256i pmt = _mm256_xor_si256(_mm256_subs_epu8(ppt, t), ushft);

	// we add (and clamp) the threshold value to all 32 pixels
	// ppt represents p + t
	ppt = _mm256_xor_si256(_mm256_adds_epu8(ppt, t), ushft);

	// Rosten's point 9 for all 32 pixels in consideration
	__m256i p9 = _mm256_xor_si256(_mm256_loadu_si256(reinterpret_cast<const __m256i*>(ptr + *offsets)), ushft);

	// Rosten's point 5
	__m256i p5 = _mm256_xor_si256(_mm256_loadu_si256(reinterpret_cast<const __m256i*>(ptr + offsets[4])), ushft);

	// Rosten's point 1
	__m256i p1 = _mm256_xor_si256(_mm256_loadu_si256(reinterpret_cast<const __m256i*>(ptr + offsets[8])), ushft);

	// Rosten's point 13
	__m256i p13 = _mm256_xor_si256(_mm256_loadu_si256(reinterpret_cast<const __m256i*>(ptr + offsets[12])), ushft);

	// compare p's against p9's
	// compare p's against p5's
	// AND the result into ppt_accum
	__m256i ppt_accum = _mm256_and_si256(_mm256_cmpgt_epi8(p9, ppt), _mm256_cmpgt_epi8(p5, ppt));

	// compare p's against p9's
	// compare p's against p5's
	// AND the result into pmt_accum
	__m256i pmt_accum = _mm256_and_si256(_mm256_cmpgt_epi8(pmt, p9), _mm256_cmpgt_epi8(pmt, p5));

	// compare p's against p5's
	// compare p's against p1's
	// AND the result, and then OR that with ppt_accum into ppt_accum
	ppt_accum = _mm256_or_si256(ppt_accum, _mm256_and_si256(_mm256_cmpgt_epi8(p5, ppt), _mm256_cmpgt_epi8(p1, ppt)));

	// compare p's against p5's
	// compare p's against p1's
	// AND the result, and then OR that with pmt_accum into pmt_accum
	pmt_accum = _mm256_or_si256(pmt_accum, _mm256_and_si256(_mm256_cmpgt_epi8(pmt, p5), _mm256_cmpgt_epi8(pmt, p1)));

	// compare p's against p1's
	// compare p's against p13's
	// AND the result, and then OR that with ppt_accum into ppt_accum
	ppt_accum = _mm256_or_si256(ppt_accum, _mm256_and_si256(_mm256_cmpgt_epi8(p1, ppt), _mm256_cmpgt_epi8(p13, ppt)));

	// compare p's against p1's
	// compare p's against p13's
	// AND the result, and then OR that with pmt_accum into pmt_accum
	pmt_accum = _mm256_or_si256(pmt_accum, _mm256_and_si256(_mm256_cmpgt_epi8(pmt, p1), _mm256_cmpgt_epi8(pmt, p13)));

	// compare p's against p13's
	// compare p's against p9's
	// AND the result, and then OR that with ppt_accum into ppt_accum
	ppt_accum = _mm256_or_si256(ppt_accum, _mm256_and_si256(_mm256_cmpgt_epi8(p13, ppt), _mm256_cmpgt_epi8(p9, ppt)));

	// compare p's against p13's
	// compare p's against p9's
	// AND the result, and then OR that with pmt_accum into pmt_accum
	pmt_accum = _mm256_or_si256(pmt_accum, _mm256_and_si256(_mm256_cmpgt_epi8(pmt, p13), _mm256_cmpgt_epi8(pmt, p9)));

	// 'full' is known by the template.
	// this and all following ternaries and ifs involving full
	// are optimized away, allowing efficient code generation for
	// the normal full 32 columns, or special handling for the last few columns if they
	// don't divide up evenly into 32
	uint32_t last_cols_mask;
	if (!full) last_cols_mask = (1 << (cols - j - 3)) - 1;

	// 'mask' now contains one bit for each element
	// which is SET if that element COULD be a corner based on the 2 consective cardinal point test
	// CAREFUL - returns signed int32_t but you NEED all 32 bits so cast immediately to unsigned
	uint32_t m = full ?
		static_cast<uint32_t>(_mm256_movemask_epi8(_mm256_or_si256(ppt_accum, pmt_accum))) :
		static_cast<uint32_t>(_mm256_movemask_epi8(_mm256_or_si256(ppt_accum, pmt_accum))) & last_cols_mask;

	// if none of the elements can be corners, bail
	if (m == 0) return;

	// if none of the left half can be corners, retreat 16 pixels to the left and bail,
	// so that after the 'continue' the total change will be forward by 16 pixels
	if (full) {
		if ((m & 0xFFFF) == 0) {
			j -= 16;
			ptr -= 16;
			return;
		}
	}

	// profiling suggests it's not worth further bailout checks for 8, 4, 2, 1
	__m256i ppt_cnt = _mm256_setzero_si256();
	__m256i pmt_cnt = _mm256_setzero_si256();
	__m256i ppt_max = _mm256_setzero_si256();
	__m256i pmt_max = _mm256_setzero_si256();

	// for each of the 24 pixels in the circle (wrapping around extra 8 at the end)
	for (int32_t k = 0; k < 24; ++k) {
		// x is a vector of the kth member of the circle
		__m256i p = _mm256_xor_si256(_mm256_loadu_si256(reinterpret_cast<const __m256i*>(ptr + offsets[k])), ushft);

		ppt_accum = _mm256_cmpgt_epi8(p, ppt);
		pmt_accum = _mm256_cmpgt_epi8(pmt, p);

		// add 1 to the chain count for all salient pixels,
		// then AND with pmt_accum to destroy any exiting chain that was broken
		// at this offset
		ppt_max = _mm256_max_epu8(ppt_max, ppt_cnt = _mm256_and_si256(_mm256_sub_epi8(ppt_cnt, ppt_accum), ppt_accum));
		pmt_max = _mm256_max_epu8(pmt_max, pmt_cnt = _mm256_and_si256(_mm256_sub_epi8(pmt_cnt, pmt_accum), pmt_accum));
	}

	// 'm' now contains one bit for whether each element
	// is a corner!
	m = full ?
		static_cast<uint32_t>(_mm256_movemask_epi8(_mm256_cmpgt_epi8(_mm256_max_epu8(ppt_max, pmt_max), consec))) :
		static_cast<uint32_t>(_mm256_movemask_epi8(_mm256_cmpgt_epi8(_mm256_max_epu8(ppt_max, pmt_max), consec))) & last_cols_mask;

	// visit each corner in the mask
	while (m) {
		const uint32_t x = _tzcnt_u32(m);
		m = _blsr_u32(m);
		if (nonmax_suppression) {
			// add it!
			corners[num_corners++] = j + x;

			// --- BEGIN COMPUTE CORNER SCORE ---

			// inlining gives measurably better performance

			const uint8_t* ptrpk = ptr + x;

			// the actual offsets value of point p
			const int16_t p = static_cast<int16_t>(*ptrpk);

			int16_t ring[24];
			for (int n = 0; n < 24; ++n) ring[n] = p - static_cast<int16_t>(ptrpk[offsets[n]]);

			int16_t* ringp = ring;

			// points 0-15
			__m256i ringv = _mm256_loadu_si256(reinterpret_cast<__m256i*>(ringp++));

			// points 1-16
			__m256i ringv2 = _mm256_loadu_si256(reinterpret_cast<__m256i*>(ringp++));
			__m256i minv = _mm256_min_epi16(ringv, ringv2);
			__m256i maxv = _mm256_max_epi16(ringv, ringv2);

			// points 2-17
			ringv = _mm256_loadu_si256(reinterpret_cast<__m256i*>(ringp++));
			minv = _mm256_min_epi16(minv, ringv);
			maxv = _mm256_max_epi16(maxv, ringv);

			// points 3-18
			ringv = _mm256_loadu_si256(reinterpret_cast<__m256i*>(ringp++));
			minv = _mm256_min_epi16(minv, ringv);
			maxv = _mm256_max_epi16(maxv, ringv);

			// points 4-19
			ringv = _mm256_loadu_si256(reinterpret_cast<__m256i*>(ringp++));
			minv = _mm256_min_epi16(minv, ringv);
			maxv = _mm256_max_epi16(maxv, ringv);

			// points 5-20
			ringv = _mm256_loadu_si256(reinterpret_cast<__m256i*>(ringp++));
			minv = _mm256_min_epi16(minv, ringv);
			maxv = _mm256_max_epi16(maxv, ringv);

			// points 6-21
			ringv = _mm256_loadu_si256(reinterpret_cast<__m256i*>(ringp++));
			minv = _mm256_min_epi16(minv, ringv);
			maxv = _mm256_max_epi16(maxv, ringv);

			// points 7-22
			ringv = _mm256_loadu_si256(reinterpret_cast<__m256i*>(ringp++));
			minv = _mm256_min_epi16(minv, ringv);
			maxv = _mm256_max_epi16(maxv, ringv);

			// points 8-23
			ringv = _mm256_loadu_si256(reinterpret_cast<__m256i*>(ringp));
			minv = _mm256_min_epi16(minv, ringv);
			maxv = _mm256_max_epi16(maxv, ringv);

			// minv now has the smallest of [0-8], [1-9], ..., [14-6], [15-7] (all 16 possible regions of 9 pixels)
			// maxv now has the largest of  [0-8], [1-9], ..., [14-6], [15-7] (all 16 possible regions of 9 pixels)

			// inside expression is just the negation of maxv
			// to get expression of absolute deviation from center offsets, resulting in
			// the greatest deviation among:
			// [0-8], [1-9], ..., [14-6], [15-7] (all 16 possible regions of 9 pixels)

			maxv = _mm256_max_epi16(minv, _mm256_sub_epi16(_mm256_setzero_si256(), maxv));
			
			// The overall single max is now found through a horizontal reduction of 'maxv'.
			// This score represents the deviation of the most deviant region of 9 pixels.
			// _mm_minpos_epu16() emits the phminposuw instruction from SSE4. Have to
			// correct for signed->unsigned, and also for max, not min. Can shift into
			// the correct space with just a single subtract operation.
			cur[j + x] = static_cast<uint8_t>(_mm_cvtsi128_si32(_mm_sub_epi16(_mm_set1_epi16(32767),
				_mm_minpos_epu16(_mm_sub_epi16(_mm_set1_epi16(32767),
					_mm_max_epi16(_mm256_extracti128_si256(maxv, 1), _mm256_castsi256_si128(maxv)))))));

			// --- END COMPUTE CORNER SCORE ---

		}
		else {
			keypoints.emplace_back(j + x, start_row + i, 0);
		}
	}
}

template <const bool nonmax_suppression, const bool first_thread, const bool last_thread>
void _KFAST(const uint8_t* __restrict const data, const int32_t cols, const int32_t start_row, const int32_t rows, const int32_t stride,
	std::vector<Keypoint>& keypoints, const uint8_t threshold) {
	keypoints.reserve(8500);

	// Rosten's circle pixels in the order 9, 8, 7, 6, 5, 4, 3, 2, 1, 16, 15, 14, 13, 12, 11, 10, then repeat 9, 8, 7, 6, 5, 4, 3, 2
	const int32_t offsets[24] = { 3 * stride, 3 * stride + 1, 2 * stride + 2, stride + 3, 3, -stride + 3, -2 * stride + 2,
		-3 * stride + 1, -3 * stride, -3 * stride - 1, -2 * stride - 2, -stride - 3, -3, stride - 3, 2 * stride - 2,
		3 * stride - 1, 3 * stride, 3 * stride + 1, 2 * stride + 2, stride + 3, 3, -stride + 3, -2 * stride + 2, -3 * stride + 1 };

	// no epu8 comparisons so in order to do them we must use epi8 comparisons and shift everything down by 128
	// add, xor, and sub 128 all do the same thing. Do it to both comparands before epi8 comparison to get
	// the equivalent unshifted epu8 comparison.
	const __m256i ushft = _mm256_set1_epi8(-128);

	// the threshold value repeated 32 times
	const __m256i t = _mm256_set1_epi8(threshold);

	// the value 8 repeated 32 times
	// will be used for comparing number of consecutive salient pixels - greater than 8 means corner!
	const __m256i consec = _mm256_set1_epi8(8);

	uint8_t* rawbuf;
	uint8_t* rowbuf[3];
	int32_t* cornerbuf[3];
	if (nonmax_suppression) {
		// allocate enough buffer for 3 rows of uint8_t and then 3 rows of int32_t
		rawbuf = reinterpret_cast<uint8_t*>(_mm_malloc(cols * 3 * (sizeof(int32_t) + sizeof(uint8_t)) + 4 * sizeof(int32_t), 4096));

		// each rowbuf entry is a pointer to a uint8_t row buffer
		rowbuf[0] = rawbuf;
		rowbuf[1] = rowbuf[0] + cols;
		rowbuf[2] = rowbuf[1] + cols;

		// each cornerbuf entry is a pointer to an int32_t row buffer
		// make sure it's aligned to a 4-byte boundary
		// and that there's space for an extra element to store num_corners
		cornerbuf[0] = reinterpret_cast<int32_t*>((reinterpret_cast<uintptr_t>(rowbuf[2] + cols) + 3) & ~3) + 1;
		cornerbuf[1] = cornerbuf[0] + cols + 1;
		cornerbuf[2] = cornerbuf[1] + cols + 1;

		// zero out the uint8_t row buffers
		memset(rowbuf[0], 0, cols * 3);
	}

	int32_t j;
	for (int32_t i = 3; i < rows - 2; ++i) {

		// ptr points to the first valid offsets in the row but hasn't retrieved it yet
		const uint8_t* ptr = data + i*stride + 3;

		uint8_t* cur = nullptr;
		int32_t* corners = nullptr;
		int32_t num_corners;
		if (nonmax_suppression) {
			// cur points to which row buffer we're in - pattern goes 0, 1, 2, 0, 1, 2, 0, 1, 2, ...
			cur = rowbuf[i % 3];

			// corners does the same thing for the int32_t buffer
			corners = cornerbuf[i % 3];

			// zero out the current rowbuffer
			memset(cur, 0, cols);
			num_corners = 0;
		}

		if (i < rows - 3) {
			// for col (3) to (cols - 35)
			// jumping forward 32 cols at a time and also moving ptr forward 32 cols each time with it
			// these calls to processCols MUST be inlined for best performance, even if your compiler thinks otherwise
			for (j = 3; j < cols - 35; j += 32, ptr += 32) {
				processCols<true, nonmax_suppression>(num_corners, ptr, j, offsets, ushft, t,
					cols, consec, corners, cur, keypoints, i, start_row);
			}
			// handle last few columns
			processCols<false, nonmax_suppression>(num_corners, ptr, j, offsets, ushft, t,
				cols, consec, corners, cur, keypoints, i, start_row);
		}

		if (nonmax_suppression) {
			corners[-1] = num_corners;

			// for first thread: skip first and last row
			// for inner threads: skip first, second, last row
			// for last thread: skip first and second row
			if ((!last_thread && i == rows - 3) || (!first_thread && i == 4) || i == 3) continue;

			// last buffered row
			const uint8_t* last = rowbuf[(i - 1) % 3];

			// last last buffered row
			const uint8_t* last2 = rowbuf[(i - 2) % 3];

			// set corners to previous buffered row
			corners = cornerbuf[(i - 1) % 3];

			// retrieve previous num_corners
			num_corners = corners[-1];
			// for each corner from the previous row
			for (int32_t k = 0; k < num_corners; ++k) {
				// corner was at col j
				j = corners[k];

				const uint8_t score = last[j];

				// if score is higher than score of all 8 surrounding pixels, add keypoint!
				// NOTE: too many branches for short-circuit evaluation to be worth it here.
				if ((score > last[j - 1]) & (score > last[j + 1]) & (score > cur[j - 1]) & (score > cur[j]) & (score > cur[j + 1]) & (score > last2[j - 1]) & (score > last2[j]) & (score > last2[j + 1])) {
					keypoints.emplace_back(j, start_row + i - 1, score);
				}
			}
		}
	}

	if (nonmax_suppression) _mm_free(rawbuf);
}

template <const bool multithreading, const bool nonmax_suppression>
void KFAST(const uint8_t* __restrict const data, const int32_t cols, const int32_t rows, const int32_t stride,
        std::vector<Keypoint>& keypoints, const uint8_t threshold) {
        if (multithreading) {
                const int32_t hw_concur = std::min(rows >> 4, static_cast<int32_t>(std::thread::hardware_concurrency()));
                std::vector<std::vector<Keypoint>> thread_kps(hw_concur);
                std::vector<std::future<void>> fut(hw_concur);

                if (hw_concur <= 1) {
                        keypoints.clear();
                        keypoints.reserve(8500);
                        _KFAST<nonmax_suppression, true, true>(data, cols, 0, rows, stride, keypoints, threshold);
                }
                else {
                        int row = (rows - 1) / hw_concur + 1;
                        fut[0] = std::async(std::launch::async, _KFAST<nonmax_suppression, true, false>, data, cols, 0, row + (3 + nonmax_suppression), stride, std::ref(thread_kps[0]), threshold);
                        int i = 1;
                        for (; i < hw_concur - 1; ++i) {
                                const int start_row = row - (3 + nonmax_suppression);
                                const int delta = (rows - row - 1) / (hw_concur - i) + 1;
                                fut[i] = std::async(std::launch::async, _KFAST<nonmax_suppression, false, false>, data + start_row*stride, cols, start_row, delta + ((3 + nonmax_suppression) << 1), stride, std::ref(thread_kps[i]), threshold);
                                row += delta;
                        }
                        int start_row = row - (3 + nonmax_suppression);
                        fut[i] = std::async(std::launch::async, _KFAST<nonmax_suppression, false, true>, data + start_row*stride, cols, start_row, rows - start_row, stride, std::ref(thread_kps[i]), threshold);
                        keypoints.clear();
                        keypoints.reserve(8500);
                        for (int j = 0; j <= i; ++j) {
                                fut[j].wait();
                                keypoints.insert(keypoints.end(), thread_kps[j].begin(), thread_kps[j].end());
                        }
                }
        }
        else {
                keypoints.clear();
                keypoints.reserve(8500);
                _KFAST<nonmax_suppression, true, true>(data, cols, 0, rows, stride, keypoints, threshold);
        }
}

