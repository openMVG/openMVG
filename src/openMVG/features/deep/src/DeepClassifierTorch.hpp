#pragma once

#include "DeepClassifier.hpp"

#include "loader.h"
#include "openMVG/features/feature.hpp"
#include "openMVG/image/image.hpp"
#include "wrapper.h"

#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/features2d.hpp>


#include <string>
#include <vector>
#include <cereal/cereal.hpp>

class DeepClassifierTorch : public DeepClassifier {
	public:
		DeepClassifierTorch(const std::string&);
		const std::vector<DeepClassifierKeypoint>& extractKeypoints(const openMVG::image::Image<unsigned char>&) override;
		const std::vector<openMVG::features::AffinePointFeature>& extractKeypointsOpenMVG(const openMVG::image::Image<unsigned char>&) override;
		float* describe(const openMVG::image::Image<unsigned char>&, const std::vector<DeepClassifierKeypoint>&) override;
		float* describeOpenMVG(const openMVG::image::Image<unsigned char>&, const std::vector<openMVG::features::AffinePointFeature>&) override;
		// JUST FOR TESTING!
		float* describeOpenMVG(const openMVG::image::Image<unsigned char>&, std::vector<cv::KeyPoint>&) override;
		virtual ~DeepClassifierTorch();
	protected:
		void extractDescriptors(std::vector<cv::Mat>&);
		void extractDescriptorsOpenMVG(std::vector<openMVG::image::Image<float>>&);
		
		int m_count;

		Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> m_descriptorsOpenMVG;
		cv::Mat m_descriptors;
		cunn::Sequential::Ptr m_net;
		THCState* m_state;
};
