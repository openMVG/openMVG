#include "DeepClassifierTHNets.hpp"

#include "openMVG/features/tbmr/tbmr.hpp"
#include "openMVG/image/image_converter.hpp"

#include <unsupported/Eigen/MatrixFunctions>

#include <opencv2/core/eigen.hpp>
#include <opencv2/core/mat.hpp>
#include <opencv2/xfeatures2d.hpp>

#include <cuda.h>
#include <cuda_runtime.h>
#include "opencv2/core/cuda.hpp"
#include "opencv2/cudafeatures2d.hpp"
#include <opencv2/features2d.hpp>
#include <opencv2/xfeatures2d.hpp>

#include <cmath>
#include <iostream>

#include "../thnets/thnets.h"

using namespace openMVG::features;
using namespace openMVG::image;

void DeepClassifierTHNets::extractPatches(const cv::Mat& image,
    const std::vector<cv::KeyPoint>& kp,
    std::vector<cv::Mat>& patches) 
{
  const size_t M = 32;
  for (auto &it : kp) {
    cv::Mat patch(M, M, CV_8UC1);
    cv::Mat buf;
    const float size = it.size;
    cv::getRectSubPix(image, cv::Size(size*1.3, size*1.3), cv::Point2f(it.pt.x, it.pt.y), buf);
    cv::resize(buf, patch, cv::Size(M, M));
    patches.push_back(patch);
  }
}

void DeepClassifierTHNets::extractPatches(const cv::Mat& image,
    const std::vector<openMVG::features::AffinePointFeature>& kp,
    std::vector<cv::Mat>& patches) 
{
  const size_t M = 32;
  for (auto &it : kp) {
    cv::Mat patch(M, M, CV_8UC1);
    cv::Mat buf;
    const float size = 2 * std::sqrt(it.l1() * it.l2() );
    cv::getRectSubPix(image, cv::Size(size*1.3, size*1.3), cv::Point2f(it.x(), it.y()), buf);
    cv::resize(buf, patch, cv::Size(M, M));
    patches.push_back(patch);
  }
}

void DeepClassifierTHNets::extractDescriptors(std::vector<cv::Mat>& patches)
{
  const size_t batch_size = 1024;
  const size_t N = patches.size();
  const unsigned int color_dim = m_network->grayscale ? 1 : 3;

  std::cout << "Extracting descriptors" << std::endl;

  for (size_t j = 0; j < ceil(N / batch_size); j++) {
    unsigned char* images[batch_size];
    size_t k = 0;
    for (size_t i = j*batch_size; i < std::min((j+1)*batch_size, N); i++, k++) {
      unsigned char* image_ptr = (unsigned char*) malloc(sizeof(unsigned char) * 32 * 32);
      memcpy(image_ptr, patches[i].data, sizeof(unsigned char) * 32 * 32);
      images[k] = image_ptr;
    }
    
    float* result;
    int outWidth, outHeight;
    THProcessImages(m_network, images, batch_size, 32, 32, 32, &result, &outWidth, &outHeight, 0);
    
    const size_t feature_dim = 128;
    if (m_descriptors.cols != feature_dim || m_descriptors.rows != N)
      m_descriptors.create(N, feature_dim, CV_32F);
    memcpy(m_descriptors.data + j * feature_dim * batch_size * sizeof(float), result, sizeof(float) * feature_dim * k);
  }
}

void DeepClassifierTHNets::extractDescriptorsOpenMVG(std::vector<Image<float>>& patches)
{
  const size_t batch_size = 15360;
  const size_t N = patches.size();

}

// Default should be to assume that we're using CUDA and full float.
DeepClassifierTHNets::DeepClassifierTHNets(const std::string& modeldir) {
  std::cout << "From inside DeepClassifierTHNets, opening: " << modeldir << std::endl;
  THInit();

  const int grayscale = 1;
  THNETWORK* net = THLoadNetwork(modeldir.c_str(), grayscale);

  if (!net) {
    THError("Couldn't load network");
  }

  const int lastlayer = 0x7fffffff;
  if (net->net->nelem > lastlayer)
    net->net->nelem = lastlayer;
  
  THUseSpatialConvolutionMM(net, 0);
  
  m_network = THCreateCudaNetwork(net);
  THFreeNetwork(net);
  m_network->grayscale = 1;
  //m_network = net;

  m_descriptors = cv::Mat::zeros(1, 1, CV_32F);

  m_count = 0;

}
// TODO: BROKEN
const std::vector<features::AffinePointFeature>& DeepClassifierTHNets::extractKeypointsOpenMVG(const Image<unsigned char>& image) {
  std::vector<features::AffinePointFeature> feats_tbmr_bright;
  std::vector<features::AffinePointFeature> feats_tbmr_dark;
/*
  tbmr::Extract_tbmr(image, feats_tbmr_bright, std::less<uint8_t>());
  tbmr::Extract_tbmr(image, feats_tbmr_dark, std::greater<uint8_t>());
  feats_tbmr_bright.insert(feats_tbmr_bright.end(), feats_tbmr_dark.begin(), feats_tbmr_dark.end());
  */
  const std::vector<features::AffinePointFeature>& feats_tbmr(feats_tbmr_bright);
  
  return feats_tbmr;
}

float* DeepClassifierTHNets::describeOpenMVG(const openMVG::image::Image<unsigned char>& image, std::vector<cv::KeyPoint>& keypoints) {
  std::vector<cv::Mat> patches;

  cv::Mat imageMat;
  cv::eigen2cv(image, imageMat);
  std::cout << "Converted eigen to cv" << std::endl;
  
  // 1. Run BRISK extractor
  /*cv::Ptr<cv::BRISK> briskDetector = cv::BRISK::create();
  std::vector<cv::KeyPoint> briefKeypoints;
  briskDetector->detect(imageMat, briefKeypoints, cv::noArray());
  briskDetector.release();
  const size_t briskSize = briefKeypoints.size();
  std::cout << "Brisk size: " << briskSize << std::endl;
*/
  // 2. Use BRIEF size to feed to SIFT
  /*cv::Ptr<cv::xfeatures2d::SIFT> siftDetector = cv::xfeatures2d::SIFT::create(briskSize);
  siftDetector->detect(imageMat, keypoints, cv::noArray());
  RetainBestKeypoints(keypoints, briskSize);
  keypoints.resize(m_testingArr[m_count]);
  m_count++;
  std::cout << "Sift keypoints size: " << keypoints.size() << std::endl;*/
  // 3. Run deep descriptors across all images
  
   cv::cuda::GpuMat img1g;
  
   img1g.upload(imageMat);
    
    {
      cv::Ptr<cv::cuda::ORB> m_orbClassifier = cv::cuda::ORB::create(30720);
      m_orbClassifier->setBlurForDescriptor(true);
  
      m_orbClassifier->detect(img1g, keypoints);
  
      keypoints.resize(keypoints.size() - (keypoints.size() % 16));
      m_orbClassifier.release();
      
    }

  extractPatches(imageMat, keypoints, patches);
  
  std::cout << "Extracted all patches" << std::endl;
  extractDescriptors(patches);

  std::cout << "Extracted all descriptors" << std::endl;

  return (float*)(m_descriptors.data);
}

float* DeepClassifierTHNets::describeOpenMVG(const openMVG::image::Image<unsigned char>& image, const std::vector<openMVG::features::AffinePointFeature>& keypoints) {
  std::vector<cv::Mat> patches;

  cv::Mat imageMat;
  cv::eigen2cv(image, imageMat);
  std::cout << "Converted eigen to cv" << std::endl;
  
  extractPatches(imageMat, keypoints, patches);
  
  std::cout << "Extracted all patches" << std::endl;
  extractDescriptors(patches);

  std::cout << "Extracted all descriptors" << std::endl;

  return (float*)(m_descriptors.data);
}
  
const std::vector<DeepClassifierKeypoint>& DeepClassifierTHNets::extractKeypoints(const Image<unsigned char>& image) {
  auto feats = extractKeypointsOpenMVG(image);

  std::vector<DeepClassifierKeypoint> deepClassifierKeypoints;
  for (auto it : feats) {
    deepClassifierKeypoints.push_back(
      DeepClassifierKeypoint(
        it.x(), it.y(),
        it.l1(), it.l2(),
        it.orientation(),
        it.a(), it.b(), it.c()
      )
    );
  }
  const std::vector<DeepClassifierKeypoint>& constDeepClassifierKeypoints(deepClassifierKeypoints);
  return constDeepClassifierKeypoints;
}

float* DeepClassifierTHNets::describe(const openMVG::image::Image<unsigned char>& image, const std::vector<DeepClassifierKeypoint>& keypoints) {
  Image<unsigned char> Icpy (image);
  std::vector<Image<float>> patches;

  //extractDescriptors(patches);

  return m_descriptorsOpenMVG.data();
}

DeepClassifierTHNets::~DeepClassifierTHNets() {
  THFreeNetwork(m_network);
}
