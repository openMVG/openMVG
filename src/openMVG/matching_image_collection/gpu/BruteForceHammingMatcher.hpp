#pragma once

#include <vector>

#include <cuda.h>
#include <cuda_runtime.h>

#include <opencv2/opencv.hpp>
#include <opencv2/core/cuda.hpp>
#include <opencv2/cudaimgproc.hpp>
#include <opencv2/cudafeatures2d.hpp>

#include "openMVG/features/latch/LatchBitMatcherMatch.hpp"

template <typename featureType, unsigned int featureLength>
class GPUBruteForceHammingMatcher {
	public:
		GPUBruteForceHammingMatcher(const float);
		void match(featureType*, featureType*, int, int);
		const std::vector<LatchBitMatcherMatch> retrieveMatches();
		~GPUBruteForceHammingMatcher();
	private:
		const float m_matchThreshold;
		std::vector<LatchBitMatcherMatch> m_goodMatches;

		cv::Ptr<cv::cuda::DescriptorMatcher> m_matcher;
		
		cudaStream_t m_stream;
};
