#include "OrbClassifierOpenMVG.hpp"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits>
#include <iostream>
#include "opencv2/core/mat.hpp"
#include "opencv2/core/eigen.hpp"
#include "opencv2/core/cuda.hpp"
#include "opencv2/core/cuda_stream_accessor.hpp"
#include "opencv2/cudaimgproc.hpp"
#include "opencv2/xfeatures2d.hpp"

/* Helper functions. */

#define checkError(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true) {
   if (code != cudaSuccess) {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

#define checkLaunchError()                                            \
do {                                                                  \
    /* Check synchronous errors, i.e. pre-launch */                   \
    cudaError_t err = cudaGetLastError();                             \
    if (cudaSuccess != err) {                                         \
        fprintf (stderr, "Cuda error in file '%s' in line %i : %s.\n",\
                 __FILE__, __LINE__, cudaGetErrorString(err) );       \
        exit(EXIT_FAILURE);                                           \
    }                                                                 \
    /* Check asynchronous errors, i.e. kernel failed (ULF) */         \
    err = cudaThreadSynchronize();                                    \
    if (cudaSuccess != err) {                                         \
        fprintf (stderr, "Cuda error in file '%s' in line %i : %s.\n",\
                 __FILE__, __LINE__, cudaGetErrorString( err) );      \
        exit(EXIT_FAILURE);                                           \
    }                                                                 \
} while (0)


/* Main class definition */

OrbClassifierOpenMVG::OrbClassifierOpenMVG() :
	m_maxKP(30720),
	m_stream(cv::cuda::Stream())
{
    m_orbClassifier = cv::cuda::ORB::create(m_maxKP);
    m_orbClassifier->setBlurForDescriptor(true);
}

std::vector<OrbClassifierKeypoint> OrbClassifierOpenMVG::identifyFeaturePointsOpenMVG(Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> img) {
    cv::Mat imgConverted;
    cv::eigen2cv(img, imgConverted);
		
		// Convert image to grayscale
    cv::cuda::GpuMat img1g;
		{
      cv::cuda::GpuMat imgGpu;
      imgGpu.upload(imgConverted, m_stream);

      imgConverted.channels() == 3 ? cv::cuda::cvtColor(imgGpu, img1g, CV_BGR2GRAY, 0, m_stream) : img1g.upload(imgConverted, m_stream);
    }
		std::vector<cv::KeyPoint> returnedKeypoints;
		cudaStream_t copiedStream = cv::cuda::StreamAccessor::getStream(m_stream);
		{
      cv::cuda::GpuMat d_keypoints;
			cv::cuda::GpuMat d_descriptors;
	
      m_orbClassifier->detectAndComputeAsync(img1g, cv::cuda::GpuMat(), d_keypoints, d_descriptors);
		  m_orbClassifier->convert(d_keypoints, returnedKeypoints);
			
			std::cout << "Hit this point. Downloading." << std::endl;

			cv::Mat h_descriptors;
			d_descriptors.download(h_descriptors);

			cv::cv2eigen(h_descriptors, m_descriptors);
		}

    return convertCVKeypointsToCustom(returnedKeypoints);
}

std::vector<OrbClassifierKeypoint> OrbClassifierOpenMVG::convertCVKeypointsToCustom(std::vector<cv::KeyPoint>& keypointsCV) {
    std::vector<OrbClassifierKeypoint> keypoints;
    for (size_t i = 0; i < keypointsCV.size(); i++) {
         OrbClassifierKeypoint kp(
            keypointsCV[i].pt.x,
            keypointsCV[i].pt.y,
            keypointsCV[i].angle * M_PI / 180.0,
            keypointsCV[i].size
        );
        keypoints.push_back(kp);
    }
    return keypoints;
}

OrbClassifierOpenMVG::~OrbClassifierOpenMVG() {

}
