/*******************************************************************
*   FeatureAngle.h
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

#include <cfloat>
#include <cmath>
#include <cstdint>
#include <immintrin.h>

constexpr float PI = 3.1415927f;

float fastAtan2(float y, float x) {
	const float ax = fabs(x);
	const float ay = fabs(y);
	float a;
	if (ax >= ay) {
		const float c = ay / (ax + FLT_MIN);
		const float cc = c*c;
		a = (((-0.0443265555479f*cc + 0.1555786518f)*cc - 0.325808397f)*cc + 0.9997878412f)*c;
	}
	else {
		const float c = ax / (ay + FLT_MIN);
		const float cc = c*c;
		a = PI * 0.5f - (((-0.0443265555479f*cc + 0.1555786518f)*cc - 0.325808397f)*cc + 0.9997878412f)*c;
	}
	if (x < 0.0f) a = PI - a;
	if (y < 0.0f) a = -a;
	return a;
}

//     0 1 2 3 4 5 6
//   +--------------
// 0 | - - x x x - -
// 1 | - x x x x x -
// 2 | x x x x x x x
// 3 | x x x o x x x
// 4 | x x x x x x x
// 5 | - x x x x x -
// 6 | - - x x x - -

static const __m128i xwt0 = _mm_setr_epi16(0, 0, -1, 0, 1, 0, 0, 0);
static const __m128i xwt1 = _mm_setr_epi16(0, -2, -1, 0, 1, 2, 0, 0);
static const __m128i xwt2 = _mm_setr_epi16(-3, -2, -1, 0, 1, 2, 3, 0);

static const __m128i ywt0 = _mm_setr_epi16(0, 0, 3, 3, 3, 0, 0, 0);
static const __m128i ywt1 = _mm_setr_epi16(0, 2, 2, 2, 2, 2, 0, 0);
static const __m128i ywt2 = _mm_setr_epi16(1, 1, 1, 1, 1, 1, 1, 0);

float featureAngle(const uint8_t* const __restrict image, const int px, const int py, const int step) {
	const uint8_t* __restrict p = image + (py - 3)*step + (px - 3);
	__m128i x = _mm_setzero_si128();
	__m128i y = _mm_setzero_si128();

	__m128i r;
	r = _mm_cvtepu8_epi16(_mm_loadl_epi64(reinterpret_cast<const __m128i*>(p)));
	x = _mm_add_epi16(x, _mm_mullo_epi16(r, xwt0));
	y = _mm_sub_epi16(y, _mm_mullo_epi16(r, ywt0));
	p += step;

	r = _mm_cvtepu8_epi16(_mm_loadl_epi64(reinterpret_cast<const __m128i*>(p)));
	x = _mm_add_epi16(x, _mm_mullo_epi16(r, xwt1));
	y = _mm_sub_epi16(y, _mm_mullo_epi16(r, ywt1));
	p += step;

	r = _mm_cvtepu8_epi16(_mm_loadl_epi64(reinterpret_cast<const __m128i*>(p)));
	x = _mm_add_epi16(x, _mm_mullo_epi16(r, xwt2));
	y = _mm_sub_epi16(y, _mm_mullo_epi16(r, ywt2));
	p += step;

	r = _mm_cvtepu8_epi16(_mm_loadl_epi64(reinterpret_cast<const __m128i*>(p)));
	x = _mm_add_epi16(x, _mm_mullo_epi16(r, xwt2));
	p += step;

	r = _mm_cvtepu8_epi16(_mm_loadl_epi64(reinterpret_cast<const __m128i*>(p)));
	x = _mm_add_epi16(x, _mm_mullo_epi16(r, xwt2));
	y = _mm_add_epi16(y, _mm_mullo_epi16(r, ywt2));
	p += step;

	r = _mm_cvtepu8_epi16(_mm_loadl_epi64(reinterpret_cast<const __m128i*>(p)));
	x = _mm_add_epi16(x, _mm_mullo_epi16(r, xwt1));
	y = _mm_add_epi16(y, _mm_mullo_epi16(r, ywt1));
	p += step;

	r = _mm_cvtepu8_epi16(_mm_loadl_epi64(reinterpret_cast<const __m128i*>(p)));
	x = _mm_add_epi16(x, _mm_mullo_epi16(r, xwt0));
	y = _mm_add_epi16(y, _mm_mullo_epi16(r, ywt0));

	x = _mm_add_epi16(x, _mm_shuffle_epi32(x, 78));
	x = _mm_hadd_epi16(x, x);
	x = _mm_add_epi16(x, _mm_shufflelo_epi16(x, 225));
	const float x_sum = static_cast<float>(static_cast<int16_t>(_mm_cvtsi128_si32(x)));

	y = _mm_add_epi16(y, _mm_shuffle_epi32(y, 78));
	y = _mm_hadd_epi16(y, y);
	y = _mm_add_epi16(y, _mm_shufflelo_epi16(y, 225));
	const float y_sum = static_cast<float>(static_cast<int16_t>(_mm_cvtsi128_si32(y)));

	return fastAtan2(y_sum, x_sum);
}