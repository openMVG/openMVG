#pragma once

/*
 * Author: Matthew Daiter
 * This class serves as a client class to access DeepDescriptors.
 * This class primarily serves to abstract away the OpenMVG interface
 * from the rest of the library.
 * */

#include <vector>

#include "DeepClassifierKeypoint.hpp"
#include "openMVG/features/feature.hpp"
#include "openMVG/image/image_container.hpp"

#include "opencv2/core.hpp"

using namespace openMVG;
using namespace image;

class DeepClassifier {
	public:
		virtual const std::vector<DeepClassifierKeypoint>& extractKeypoints(const Image<unsigned char>&) = 0;
		virtual const std::vector<openMVG::features::AffinePointFeature>& extractKeypointsOpenMVG(const openMVG::image::Image<unsigned char>&) = 0;
		virtual float* describe(const openMVG::image::Image<unsigned char>&, const std::vector<DeepClassifierKeypoint>&) = 0;
		virtual float* describeOpenMVG(const openMVG::image::Image<unsigned char>&, const std::vector<openMVG::features::AffinePointFeature>&) = 0;
		// Just for testing!
		virtual float* describeOpenMVG(const openMVG::image::Image<unsigned char>&, std::vector<cv::KeyPoint>&) = 0;
		virtual ~DeepClassifier() { }
	protected:
		const bool m_useGPU = true;
};
