#ifndef ORB_CLASSIFIER_OPENMVG_H
#define ORB_CLASSIFIER_OPENMVG_H

#include <tuple>
#include <vector>

#include <opencv2/core/types.hpp>
#include <opencv2/cudafeatures2d.hpp>
#include <Eigen/Core>

#include "OrbClassifierKeypoint.hpp"

class OrbClassifierOpenMVG  {
    public:
        OrbClassifierOpenMVG();
        
        std::vector<OrbClassifierKeypoint> identifyFeaturePointsOpenMVG(Eigen::Matrix<unsigned char, Eigen::Dynamic,
        Eigen::Dynamic, Eigen::RowMajor>);

				Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> returnDescriptors() { return m_descriptors; }
        
        ~OrbClassifierOpenMVG();

		protected:
				std::vector<OrbClassifierKeypoint> convertCVKeypointsToCustom(std::vector<cv::KeyPoint>&);

		private:
				cv::cuda::Stream m_stream;
				cv::Ptr<cv::cuda::ORB> m_orbClassifier;
				Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic>  m_descriptors;
				const size_t m_maxKP;
};
#endif
