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
		
		Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> m_descriptorsOpenMVG;
		cv::Mat m_descriptors;
		cunn::Sequential::Ptr m_net;
		THCState* m_state;

	// ONLY FOR TESTING PURPOSES FOR CVPR16 PAPER.
				int m_testingArr[48] = { 
					409,
				 	254,
					407,
					616,
					639,
					898,
					713,
					312,
					227,
					150,
					126,
					101,
					2332,
					2343,
					2018,
					1363,
					1153,
					1334,
					505,
					618,
					695,
					664,
					644,
					648,
					730,
					581,
					447,
					399,
					323,
					219,
					4780,
					5034,
					4530,
					2625,
					1268,
					696,
					1611,
					1517,
					1470,
					1439,
					1258,
					798,
					3678,
					2658,
					2443,
					2560,
					2555,
					2770
				};
			int m_count;

};
