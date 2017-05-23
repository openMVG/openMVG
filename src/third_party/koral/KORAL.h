/*******************************************************************
*   KORAL.h
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

#include "CLATCH.h"
#include "CUDALERP.h"
#include "FeatureAngle.h"
#include "KFAST.h"

#include <cstdint>
#include <cstring>
#include <cuda_runtime.h>
#include <vector>

class KORAL {
	// public member variables
public:
	std::vector<Keypoint> kps;
	std::vector<uint64_t> desc;

	// private member variables
private:

	struct Level {
		uint8_t* d_img;
		size_t pitch;
		const uint8_t* h_img;
		uint32_t w;
		uint32_t h;
		size_t total;

		Level() : d_img(nullptr), h_img(nullptr) {}
	};

	Level* levels;
	cudaTextureObject_t *all_tex;
	cudaChannelFormatDesc chandesc_img;
	struct cudaTextureDesc texdesc_img;
	cudaTextureObject_t d_trip_tex;
	const float scale_factor;
	const uint8_t scale_levels;
	uint64_t* d_desc;
	Keypoint* d_kps;

	// public methods
public:
	KORAL(const float _scale_factor, const uint8_t _scale_levels) : scale_factor(_scale_factor), scale_levels(_scale_levels) {
		// setting cache and shared modes
		cudaDeviceSetCacheConfig(cudaFuncCachePreferEqual);
		cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeFourByte);

		// allocating and transferring triplets and binding to texture object
		// for CLATCH
		uint32_t* d_triplets;
		cudaMalloc(&d_triplets, 2048 * sizeof(uint16_t));
		cudaMemcpy(d_triplets, triplets, 2048 * sizeof(uint16_t), cudaMemcpyHostToDevice);
		cudaChannelFormatDesc chandesc_trip = cudaCreateChannelDesc(16, 16, 16, 16, cudaChannelFormatKindUnsigned);
		cudaArray* d_trip_arr;
		cudaMallocArray(&d_trip_arr, &chandesc_trip, 512);
		cudaMemcpyToArray(d_trip_arr, 0, 0, d_triplets, 2048 * sizeof(uint16_t), cudaMemcpyHostToDevice);
		struct cudaResourceDesc resdesc_trip;
		memset(&resdesc_trip, 0, sizeof(resdesc_trip));
		resdesc_trip.resType = cudaResourceTypeArray;
		resdesc_trip.res.array.array = d_trip_arr;
		struct cudaTextureDesc texdesc_trip;
		memset(&texdesc_trip, 0, sizeof(texdesc_trip));
		texdesc_trip.addressMode[0] = cudaAddressModeClamp;
		texdesc_trip.filterMode = cudaFilterModePoint;
		texdesc_trip.readMode = cudaReadModeElementType;
		texdesc_trip.normalizedCoords = 0;
		cudaCreateTextureObject(&d_trip_tex, &resdesc_trip, &texdesc_trip, nullptr);

		memset(&texdesc_img, 0, sizeof(texdesc_img));
		texdesc_img.addressMode[0] = cudaAddressModeClamp;
		texdesc_img.addressMode[1] = cudaAddressModeClamp;
		texdesc_img.filterMode = cudaFilterModePoint;
		texdesc_img.normalizedCoords = 0;

		chandesc_img = cudaCreateChannelDesc(8, 0, 0, 0, cudaChannelFormatKindUnsigned);

		levels = new Level[scale_levels];
		all_tex = new cudaTextureObject_t[scale_levels];
	}

	~KORAL() {
		delete[] levels;
		delete[] all_tex;
	}

	void go(const uint8_t* image, const uint32_t width, const uint32_t height, const uint8_t KFAST_thresh) {
		kps.clear();
		levels[0].h_img = image;
		levels[0].w = width;
		levels[0].h = height;
		levels[0].total = static_cast<size_t>(width) * static_cast<size_t>(height);

		// allocating and transferring original image as cudaArray
		// and binding to texture object, one as normalized float (for LERP),
		// one as ElementType (for CLATCH)
		cudaTextureObject_t d_img_tex_nf;
		{
			cudaArray* d_img_array;
			cudaMallocArray(&d_img_array, &chandesc_img, width, height, cudaArrayTextureGather);
			cudaMemcpyToArray(d_img_array, 0, 0, image, levels[0].total, cudaMemcpyHostToDevice);
			struct cudaResourceDesc resdesc_img;
			memset(&resdesc_img, 0, sizeof(resdesc_img));
			resdesc_img.resType = cudaResourceTypeArray;
			resdesc_img.res.array.array = d_img_array;

			// first as normalized float
			texdesc_img.readMode = cudaReadModeNormalizedFloat;
			cudaCreateTextureObject(&d_img_tex_nf, &resdesc_img, &texdesc_img, nullptr);

			// then as ElementType
			texdesc_img.readMode = cudaReadModeElementType;
			cudaCreateTextureObject(&all_tex[0], &resdesc_img, &texdesc_img, nullptr);
		}

		cudaStream_t* stream = new cudaStream_t[scale_levels - 1];
		for (int i = 0; i < scale_levels - 1; ++i) {
			cudaStreamCreate(stream + i);
		}

		// prepare 7 more scales as 2D pitched linear
		// and bind to ElementType textures
		float f = 1.0f;
		for (int i = 1; i < scale_levels; ++i) {
			f *= scale_factor;
			levels[i].w = static_cast<uint32_t>(static_cast<float>(width) / f + 0.5f);
			levels[i].h = static_cast<uint32_t>(static_cast<float>(height) / f + 0.5f);
			levels[i].total = static_cast<size_t>(levels[i].w)*static_cast<size_t>(levels[i].h);

			levels[i].h_img = reinterpret_cast<uint8_t*>(malloc(levels[i].total + 1));

			cudaMallocPitch(&levels[i].d_img, &levels[i].pitch, levels[i].w, levels[i].h);

			struct cudaResourceDesc resdesc_img;
			memset(&resdesc_img, 0, sizeof(resdesc_img));
			resdesc_img.resType = cudaResourceTypePitch2D;
			resdesc_img.res.pitch2D.desc = chandesc_img;
			resdesc_img.res.pitch2D.devPtr = levels[i].d_img;
			resdesc_img.res.pitch2D.height = levels[i].h;
			resdesc_img.res.pitch2D.pitchInBytes = levels[i].pitch;
			resdesc_img.res.pitch2D.width = levels[i].w;
			cudaCreateTextureObject(&all_tex[i], &resdesc_img, &texdesc_img, nullptr);

			// GPU: non-blocking launch of resize kernels
			CUDALERP(d_img_tex_nf, f, f, levels[i].d_img, levels[i].pitch, levels[i].w, levels[i].h, stream[i - 1]);
		}

		// meanwhile, CPU, get started on KFAST

		// bring in downscale results from GPU (except for first level) and operate on them
		// as they arrive
		for (uint8_t i = 0; i < scale_levels; ++i) {
			std::vector<Keypoint> local_kps;
			if (i) {
				cudaMemcpy2DAsync(const_cast<uint8_t*>(levels[i].h_img), levels[i].w, levels[i].d_img, levels[i].pitch, levels[i].w, levels[i].h, cudaMemcpyDeviceToHost, stream[i - 1]);
				cudaStreamSynchronize(stream[i - 1]);
			}
			KFAST<true, true>(levels[i].h_img, levels[i].w, levels[i].h, levels[i].w, local_kps, KFAST_thresh);

			// set scale and compute angles
			for (auto& kp : local_kps) {
				kp.scale = i;
				kp.angle = featureAngle(levels[i].h_img, kp.x, kp.y, static_cast<int>(levels[i].w));
			}
			//std::cout << "Got " << local_kps.size() << " keypoints from level " << +i << '.' << std::endl;
			kps.insert(kps.end(), local_kps.begin(), local_kps.end());
		}

		// Describe

		// allocating space for descriptors
		cudaMalloc(&d_desc, 64 * kps.size());

		// allocating and transferring keypoints and binding to texture object
		cudaMalloc(&d_kps, kps.size() * sizeof(Keypoint));
		cudaMemcpy(d_kps, kps.data(), kps.size() * sizeof(Keypoint), cudaMemcpyHostToDevice);

		// transfer tex
		cudaTextureObject_t* d_all_tex;
		cudaMalloc(&d_all_tex, scale_levels * sizeof(cudaTextureObject_t));
		cudaMemcpy(d_all_tex, all_tex, scale_levels * sizeof(cudaTextureObject_t), cudaMemcpyHostToDevice);

		CLATCH(d_all_tex, d_trip_tex, d_kps, static_cast<int>(kps.size()), d_desc);

		// transfer descriptors

		desc.clear();
		desc.resize(8 * kps.size());
		cudaMemcpy(&desc[0], d_desc, 64 * kps.size(), cudaMemcpyDeviceToHost);

		//for (int i = 0; i < scale_levels; ++i) {
		//	std::vector<cv::KeyPoint> converted_kps;
		//	for (const auto& kp : kps) if (kp.scale == i) converted_kps.emplace_back(static_cast<float>(kp.x), static_cast<float>(kp.y), 14.0f, 180.0f/3.1415926535f*kp.angle, static_cast<float>(kp.score));
		//	cv::Mat image_with_kps;
		//	cv::Mat raw_img(static_cast<int>(levels[i].h), static_cast<int>(levels[i].w), CV_8U, const_cast<uint8_t*>(levels[i].h_img), levels[i].w);
		//	cv::drawKeypoints(raw_img, converted_kps, image_with_kps, { 255.0, 0.0, 0.0 }, 4);
		//	cv::namedWindow(std::to_string(i), CV_WINDOW_AUTOSIZE);
		//	cv::imshow(std::to_string(i), image_with_kps);
		//}
		//cv::waitKey(0);

		cudaDeviceSynchronize();
		//cudaDeviceReset();
	}

	// private methods
private:



};