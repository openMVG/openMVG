#pragma once

#include "DeepClassifier.hpp"

#include "openMVG/features/feature.hpp"
#include "openMVG/image/image.hpp"

#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/features2d.hpp>


#include <string>
#include <vector>
#include <cereal/cereal.hpp>

extern "C" {
	#include "../thnets/thnetwork.h"
}

enum NetworkType {
	CUDA,
	OpenCL,
	CPU
};

enum FloatType {
	FULL,
	HALF
};

class DeepClassifierTHNets : public DeepClassifier {
	public:
		DeepClassifierTHNets(const std::string&);
		const std::vector<DeepClassifierKeypoint>& extractKeypoints(const openMVG::image::Image<unsigned char>&) override;
		const std::vector<openMVG::features::AffinePointFeature>& extractKeypointsOpenMVG(const openMVG::image::Image<unsigned char>&) override;
		float* describe(const openMVG::image::Image<unsigned char>&, const std::vector<DeepClassifierKeypoint>&) override;
		float* describeOpenMVG(const openMVG::image::Image<unsigned char>&, const std::vector<openMVG::features::AffinePointFeature>&) override;
		// JUST FOR TESTING!
		float* describeOpenMVG(const openMVG::image::Image<unsigned char>&, std::vector<cv::KeyPoint>&) override;
		virtual ~DeepClassifierTHNets();
	protected:
		void extractDescriptors(std::vector<cv::Mat>&);
		void extractDescriptorsOpenMVG(std::vector<openMVG::image::Image<float>>&);

		THNETWORK* m_network;
		
		int m_count;

		void extractPatches(const cv::Mat&, const std::vector<cv::KeyPoint>&, std::vector<cv::Mat>&);
		void extractPatches(const cv::Mat&, const std::vector<openMVG::features::AffinePointFeature>&, std::vector<cv::Mat>&);

		Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> m_descriptorsOpenMVG;
		cv::Mat m_descriptors;
};
