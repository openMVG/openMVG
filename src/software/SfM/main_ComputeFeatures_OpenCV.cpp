
// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


#include "openMVG/image/image.hpp"
#include "openMVG/sfm/sfm.hpp"

/// Feature/Regions & Image describer interfaces
#include "openMVG/features/features.hpp"
#include <cereal/archives/json.hpp>

#include "openMVG/system/timer.hpp"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"
#include "third_party/progress/progress.hpp"

/// OpenCV Includes
#include <opencv2/opencv.hpp>
#include "opencv2/core/eigen.hpp"
#ifdef USE_OCVSIFT
#include "opencv2/xfeatures2d.hpp"
#endif

#include <cstdlib>
#include <fstream>
#include <algorithm>

using namespace openMVG;
using namespace openMVG::image;
using namespace openMVG::features;
using namespace openMVG::sfm;
using namespace std;

enum eGeometricModel
{
  FUNDAMENTAL_MATRIX = 0,
  ESSENTIAL_MATRIX   = 1,
  HOMOGRAPHY_MATRIX  = 2
};

enum ePairMode
{
  PAIR_EXHAUSTIVE = 0,
  PAIR_CONTIGUOUS = 1,
  PAIR_FROM_FILE  = 2
};

///
//- Create an Image_describer interface that use an OpenCV feature extraction method
// i.e. with the AKAZE detector+descriptor
//--/!\ If you use a new Regions type you define and register it in
//   "openMVG/features/regions_factory.hpp" file.
///
// Reuse the existing AKAZE floating point Keypoint.
typedef features::AKAZE_Float_Regions AKAZE_OpenCV_Regions;
// Define the Interface
class AKAZE_OCV_Image_describer : public Image_describer
{
public:
  AKAZE_OCV_Image_describer():Image_describer(){}


  bool Set_configuration_preset(EDESCRIBER_PRESET preset)
  {
    return false;
  }
  /**
  @brief Detect regions on the image and compute their attributes (description)
  @param image Image.
  @param regions The detected regions and attributes (the caller must delete the allocated data)
  @param mask 8-bit gray image for keypoint filtering (optional).
     Non-zero values depict the region of interest.
  */
  bool Describe(const Image<unsigned char>& image,
    std::unique_ptr<Regions> &regions,
    const Image<unsigned char> * mask = NULL)
  {
    cv::Mat img;
    cv::eigen2cv(image.GetMat(), img);

    std::vector< cv::KeyPoint > vec_keypoints;
    cv::Mat m_desc;

    cv::Ptr<cv::Feature2D> extractor = cv::AKAZE::create(cv::AKAZE::DESCRIPTOR_KAZE);
    extractor->detectAndCompute(img, cv::Mat(), vec_keypoints, m_desc);

    if (!vec_keypoints.empty())
    {
      Allocate(regions);

      // Build alias to cached data
      AKAZE_OpenCV_Regions * regionsCasted = dynamic_cast<AKAZE_OpenCV_Regions*>(regions.get());
      // reserve some memory for faster keypoint saving
      regionsCasted->Features().reserve(vec_keypoints.size());
      regionsCasted->Descriptors().reserve(vec_keypoints.size());

      typedef Descriptor<float, 64> DescriptorT;
      DescriptorT descriptor;
      int cpt = 0;
      for(std::vector< cv::KeyPoint >::const_iterator i_keypoint = vec_keypoints.begin();
        i_keypoint != vec_keypoints.end(); ++i_keypoint, ++cpt){

        SIOPointFeature feat((*i_keypoint).pt.x, (*i_keypoint).pt.y, (*i_keypoint).size, (*i_keypoint).angle);
        regionsCasted->Features().push_back(feat);

        memcpy(descriptor.getData(),
               m_desc.ptr<typename DescriptorT::bin_type>(cpt),
               DescriptorT::static_size*sizeof(typename DescriptorT::bin_type));
        regionsCasted->Descriptors().push_back(descriptor);
      }
    }
    return true;
  };

  /// Allocate Regions type depending of the Image_describer
  void Allocate(std::unique_ptr<Regions> &regions) const
  {
    regions.reset( new AKAZE_OpenCV_Regions );
  }

  template<class Archive>
  void serialize( Archive & ar )
  {
  }
};
#include <cereal/cereal.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/archives/json.hpp>
CEREAL_REGISTER_TYPE_WITH_NAME(AKAZE_OCV_Image_describer, "AKAZE_OCV_Image_describer");

#ifdef USE_OCVSIFT

class SIFT_OPENCV_Params
{
public:
  SIFT_OPENCV_Params() {}
  ~SIFT_OPENCV_Params() {}
  
  bool Set_configuration_preset(EDESCRIBER_PRESET preset)
  {
      switch(preset)
      {
        case LOW_PRESET:
          contrastThreshold = 0.01;
          maxTotalKeypoints = 1000;
          break;
        case MEDIUM_PRESET:
          contrastThreshold = 0.005;
          maxTotalKeypoints = 5000;
          break;
        case NORMAL_PRESET:
          contrastThreshold = 0.005;
          edgeThreshold = 15;
          maxTotalKeypoints = 10000;
          break;
        case HIGH_PRESET:
          contrastThreshold = 0.005;
          edgeThreshold = 20;
          maxTotalKeypoints = 20000;
          break;
        case ULTRA_PRESET:
          contrastThreshold = 0.005;
          edgeThreshold = 20;
          maxTotalKeypoints = 40000;
          break;
      }
      return true;
  }

  template<class Archive>
  void serialize( Archive & ar )
  {
    ar(
      cereal::make_nvp("grid_size", gridSize),
      cereal::make_nvp("max_total_keypoints", maxTotalKeypoints),
      cereal::make_nvp("n_octave_layers", nOctaveLayers),
      cereal::make_nvp("contrast_threshold", contrastThreshold),
      cereal::make_nvp("edge_threshold", edgeThreshold),
      cereal::make_nvp("sigma", sigma));
      // cereal::make_nvp("root_sift", root_sift));
  }

  // Parameters
  std::size_t gridSize = 4;
  std::size_t maxTotalKeypoints = 1000;
  int nOctaveLayers = 6;  // default opencv value is 3
  double contrastThreshold = 0.04;  // default opencv value is 0.04
  double edgeThreshold = 10;
  double sigma = 1.6;
  // bool rootSift = true;
};

///
//- Create an Image_describer interface that use and OpenCV extraction method
// i.e. with the SIFT detector+descriptor
// Regions is the same as classic SIFT : 128 unsigned char
class SIFT_OPENCV_Image_describer : public Image_describer
{
public:
  SIFT_OPENCV_Image_describer() : Image_describer() {}

  ~SIFT_OPENCV_Image_describer() {}

  bool Set_configuration_preset(EDESCRIBER_PRESET preset)
  {
    return _params.Set_configuration_preset(preset);
  }
  /**
  @brief Detect regions on the image and compute their attributes (description)
  @param image Image.
  @param regions The detected regions and attributes (the caller must delete the allocated data)
  @param mask 8-bit gray image for keypoint filtering (optional).
     Non-zero values depict the region of interest.
  */
  bool Describe(const image::Image<unsigned char>& image,
    std::unique_ptr<Regions> &regions,
    const image::Image<unsigned char> * mask = NULL)
  {
    // Convert for opencv
    cv::Mat img;
    cv::eigen2cv(image.GetMat(), img);

    // Create a SIFT detector
    std::vector< cv::KeyPoint > v_keypoints;
    cv::Mat m_desc;
    std::size_t maxDetect = 0; // No max value by default
    if(_params.maxTotalKeypoints)
      if(!_params.gridSize)  // If no grid filtering, use opencv to limit the number of features
        maxDetect = _params.maxTotalKeypoints;

    cv::Ptr<cv::Feature2D> siftdetector = cv::xfeatures2d::SIFT::create(maxDetect, _params.nOctaveLayers, _params.contrastThreshold, _params.edgeThreshold, _params.sigma);

    // Detect SIFT keypoints
    auto detect_start = std::chrono::steady_clock::now();
    siftdetector->detect(img, v_keypoints);
    auto detect_end = std::chrono::steady_clock::now();
    auto detect_elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(detect_end - detect_start);

    std::cout << "SIFT: contrastThreshold: " << _params.contrastThreshold << ", edgeThreshold: " << _params.edgeThreshold << std::endl;
    std::cout << "Detect SIFT: " << detect_elapsed.count() << " milliseconds." << std::endl;
    std::cout << "Image size: " << img.cols << " x " << img.rows << std::endl;
    std::cout << "Grid size: " << _params.gridSize << ", maxTotalKeypoints: " << _params.maxTotalKeypoints << std::endl;
    std::cout << "Number of detected features: " << v_keypoints.size() << std::endl;

    // cv::KeyPoint::response: the response by which the most strong keypoints have been selected.
    // Can be used for the further sorting or subsampling.
    std::sort(v_keypoints.begin(), v_keypoints.end(), [](const cv::KeyPoint& a, const cv::KeyPoint& b) { return a.size > b.size; });

    // Grid filtering of the keypoints to ensure a global repartition
    if(_params.gridSize && _params.maxTotalKeypoints)
    {
      // Only filter features if we have more features than the maxTotalKeypoints
      if(v_keypoints.size() > _params.maxTotalKeypoints)
      {
        std::vector< cv::KeyPoint > filtered_keypoints;
        std::vector< cv::KeyPoint > rejected_keypoints;
        filtered_keypoints.reserve(std::min(v_keypoints.size(), _params.maxTotalKeypoints));
        rejected_keypoints.reserve(v_keypoints.size());

        cv::Mat countFeatPerCell(_params.gridSize, _params.gridSize, cv::DataType<std::size_t>::type, cv::Scalar(0));
        const std::size_t keypointsPerCell = _params.maxTotalKeypoints / countFeatPerCell.total();
        const double regionWidth = image.Width() / double(countFeatPerCell.cols);
        const double regionHeight = image.Height() / double(countFeatPerCell.rows);

        std::cout << "Grid filtering -- keypointsPerCell: " << keypointsPerCell
                  << ", regionWidth: " << regionWidth
                  << ", regionHeight: " << regionHeight << std::endl;

        for(const cv::KeyPoint& keypoint: v_keypoints)
        {
          const std::size_t cellX = std::min(std::size_t(keypoint.pt.x / regionWidth), _params.gridSize);
          const std::size_t cellY = std::min(std::size_t(keypoint.pt.y / regionHeight), _params.gridSize);
          // std::cout << "- keypoint.pt.x: " << keypoint.pt.x << ", keypoint.pt.y: " << keypoint.pt.y << std::endl;
          // std::cout << "- cellX: " << cellX << ", cellY: " << cellY << std::endl;
          // std::cout << "- countFeatPerCell: " << countFeatPerCell << std::endl;
          // std::cout << "- gridSize: " << _params.gridSize << std::endl;

          const std::size_t count = countFeatPerCell.at<std::size_t>(cellX, cellY);
          countFeatPerCell.at<std::size_t>(cellX, cellY) = count + 1;
          if(count < keypointsPerCell)
            filtered_keypoints.push_back(keypoint);
          else
            rejected_keypoints.push_back(keypoint);
        }
        // If we don't have enough features (less than maxTotalKeypoints) after the grid filtering (empty regions in the grid for example).
        // We add the best other ones, without repartition constraint.
        if( filtered_keypoints.size() < _params.maxTotalKeypoints )
        {
          const std::size_t remainingElements = std::min(rejected_keypoints.size(), _params.maxTotalKeypoints - filtered_keypoints.size());
          std::cout << "Grid filtering -- Copy remaining points: " << remainingElements << std::endl;
          filtered_keypoints.insert(filtered_keypoints.end(), rejected_keypoints.begin(), rejected_keypoints.begin() + remainingElements);
        }

        v_keypoints.swap(filtered_keypoints);
      }
    }
    std::cout << "Number of features: " << v_keypoints.size() << std::endl;

    // Compute SIFT descriptors
    auto desc_start = std::chrono::steady_clock::now();
    siftdetector->compute(img, v_keypoints, m_desc);
    auto desc_end = std::chrono::steady_clock::now();
    auto desc_elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(desc_end - desc_start);
    std::cout << "Compute descriptors: " << desc_elapsed.count() << " milliseconds." << std::endl;

    Allocate(regions);

    // Build alias to cached data
    SIFT_Regions * regionsCasted = dynamic_cast<SIFT_Regions*>(regions.get());
    // reserve some memory for faster keypoint saving
    regionsCasted->Features().reserve(v_keypoints.size());
    regionsCasted->Descriptors().reserve(v_keypoints.size());

    // Prepare a column vector with the sum of each descriptor
    cv::Mat m_siftsum;
    cv::reduce(m_desc, m_siftsum, 1, cv::REDUCE_SUM);

    // Copy keypoints and descriptors in the regions
    int cpt = 0;
    for(std::vector< cv::KeyPoint >::const_iterator i_kp = v_keypoints.begin();
        i_kp != v_keypoints.end();
        ++i_kp, ++cpt)
    {
      SIOPointFeature feat((*i_kp).pt.x, (*i_kp).pt.y, (*i_kp).size, (*i_kp).angle);
      regionsCasted->Features().push_back(feat);

      Descriptor<unsigned char, 128> desc;
      for(int j = 0; j < 128; j++)
      {
        desc[j] = static_cast<unsigned char>(512.0*sqrt(m_desc.at<float>(cpt, j)/m_siftsum.at<float>(cpt, 0)));
      }
      regionsCasted->Descriptors().push_back(desc);
    }

    return true;
  };

  /// Allocate Regions type depending of the Image_describer
  void Allocate(std::unique_ptr<Regions> &regions) const
  {
    regions.reset( new SIFT_Regions );
  }

  template<class Archive>
  void serialize( Archive & ar )
  {
    ar(cereal::make_nvp("params", _params));
  }

private:
  SIFT_OPENCV_Params _params;
};

CEREAL_REGISTER_TYPE_WITH_NAME(SIFT_OPENCV_Image_describer, "SIFT_OPENCV_Image_describer");
#endif //USE_OCVSIFT

/// Compute between the Views
/// Compute view image description (feature & descriptor extraction using OpenCV)
/// Compute putative local feature matches (descriptor matching)
/// Compute geometric coherent feature matches (robust model estimation from putative matches)
/// Export computed data
int main(int argc, char **argv)
{
  CmdLine cmd;

  std::string sSfM_Data_Filename;
  std::string sOutDir = "";
  bool bForce = false;
#ifdef USE_OCVSIFT
  std::string sImage_Describer_Method = "AKAZE_OPENCV";
#endif
  std::string sFeaturePreset = "NORMAL";
  int rangeStart = -1;
  int rangeSize = 1;

  // required
  cmd.add( make_option('i', sSfM_Data_Filename, "input_file") );
  cmd.add( make_option('o', sOutDir, "outdir") );
  // Optional
  cmd.add( make_option('f', bForce, "force") );
#ifdef USE_OCVSIFT
  cmd.add( make_option('m', sImage_Describer_Method, "describerMethod") );
#endif
  cmd.add( make_option('p', sFeaturePreset, "describerPreset") );
  cmd.add( make_option('s', rangeStart, "range_start") );
  cmd.add( make_option('r', rangeSize, "range_size") );

  try {
      if (argc == 1) throw std::string("Invalid command line parameter.");
      cmd.process(argc, argv);
  } catch(const std::string& s) {
      std::cerr << "Usage: " << argv[0] << '\n'
      << "[-i|--input_file]: a SfM_Data file \n"
      << "[-o|--outdir] path \n"
      << "\n[Optional]\n"
      << "[-f|--force: Force to recompute data]\n"
#ifdef USE_OCVSIFT
      << "[-m|--describerMethod\n"
      << "  (method to use to describe an image):\n"
      << "   AKAZE_OPENCV (default),\n"
      << "   SIFT_OPENCV: SIFT FROM OPENCV,\n"
#endif
      << "[-p|--describerPreset]\n"
      << "  (used to control the Image_describer configuration):\n"
      << "   LOW,\n"
      << "   MEDIUM,\n"
      << "   NORMAL (default),\n"
      << "   HIGH,\n"
      << "   ULTRA: !!Can take long time!!\n"
      << "[-s]--range_start] range image index start\n"
      << "[-r]--range_size] range size\n"
      << std::endl;

      std::cerr << s << std::endl;
      return EXIT_FAILURE;
  }

  std::cout << " You called : " <<std::endl
            << argv[0] << std::endl
            << "--input_file " << sSfM_Data_Filename << std::endl
            << "--outdir " << sOutDir << std::endl
#ifdef USE_OCVSIFT
            << "--describerMethod " << sImage_Describer_Method << std::endl
#endif
            << "--describerPreset " << sFeaturePreset << std::endl
            << "--range_start " << rangeStart << std::endl
            << "--range_size " << rangeSize << std::endl;

  if (sOutDir.empty())  {
    std::cerr << "\nIt is an invalid output directory" << std::endl;
    return EXIT_FAILURE;
  }

  // Create output dir
  if (!stlplus::folder_exists(sOutDir))
  {
    if (!stlplus::folder_create(sOutDir))
    {
      std::cerr << "Cannot create output directory" << std::endl;
      return EXIT_FAILURE;
    }
  }

  //---------------------------------------
  // a. Load input scene
  //---------------------------------------
  SfM_Data sfm_data;
  if (!Load(sfm_data, sSfM_Data_Filename, ESfM_Data(VIEWS|INTRINSICS))) {
    std::cerr << std::endl
      << "The input file \""<< sSfM_Data_Filename << "\" cannot be read" << std::endl;
    return false;
  }

  // Init the image_describer
  // - retrieve the used one in case of pre-computed features
  // - else create the desired one

  using namespace openMVG::features;
  std::unique_ptr<Image_describer> image_describer;

  const std::string sImage_describer = stlplus::create_filespec(sOutDir, "image_describer", "json");
  if (stlplus::is_file(sImage_describer))
  {
    // Dynamically load the image_describer from the file (will restore old used settings)
    std::ifstream stream(sImage_describer.c_str());
    if (!stream.is_open())
      return false;

    cereal::JSONInputArchive archive(stream);
    archive(cereal::make_nvp("image_describer", image_describer));
  }
  else
  {
#ifdef USE_OCVSIFT
    if (sImage_Describer_Method == "AKAZE_OPENCV")
    {
      image_describer.reset(new AKAZE_OCV_Image_describer);
    }
    else
    if (sImage_Describer_Method == "SIFT_OPENCV")
    {
      image_describer.reset(new SIFT_OPENCV_Image_describer());
    }
    else
    {
      std::cerr << "Unknown image describer method." << std::endl;
      return EXIT_FAILURE;
    }
#else
    image_describer.reset(new AKAZE_OCV_Image_describer);
#endif

    image_describer->Set_configuration_preset(sFeaturePreset);

    // Export the used Image_describer and region type for:
    // - dynamic future regions computation and/or loading
    {
      std::ofstream stream(sImage_describer.c_str());
      if (!stream.is_open())
        return false;

      cereal::JSONOutputArchive archive(stream);
      archive(cereal::make_nvp("image_describer", image_describer));
      std::unique_ptr<Regions> regionsType;
      image_describer->Allocate(regionsType);
      archive(cereal::make_nvp("regions_type", regionsType));
    }
  }

  // Feature extraction routines
  // For each View of the SfM_Data container:
  // - if regions file exist continue,
  // - if no file, compute features
  {
    system::Timer timer;
    Image<unsigned char> imageGray;
    C_Progress_display my_progress_bar( sfm_data.GetViews().size(),
      std::cout, "\n- EXTRACT FEATURES -\n" );
    Views::const_iterator iterViews = sfm_data.views.begin();
    Views::const_iterator iterViewsEnd = sfm_data.views.end();
    if(rangeStart != -1)
    {
      if(rangeStart < 0 || rangeStart > sfm_data.views.size())
      {
        std::cerr << "Bad specific index" << std::endl;
        return EXIT_FAILURE;
      }
      if(rangeSize < 0 || rangeSize > sfm_data.views.size())
      {
        std::cerr << "Bad range size. " << std::endl;
        return EXIT_FAILURE;
      }
      if(rangeStart + rangeSize > sfm_data.views.size())
        rangeSize = sfm_data.views.size() - rangeStart;

      std::advance(iterViews, rangeStart);
      iterViewsEnd = iterViews;
      std::advance(iterViewsEnd, rangeSize);
    }
    for(;
        iterViews != iterViewsEnd;
        ++iterViews, ++my_progress_bar)
    {
      const View * view = iterViews->second.get();
      const std::string sView_filename = stlplus::create_filespec(sfm_data.s_root_path,
        view->s_Img_path);
      const std::string sFeat = stlplus::create_filespec(sOutDir,
        stlplus::basename_part(std::to_string(view->id_view)), "feat");
      const std::string sDesc = stlplus::create_filespec(sOutDir,
        stlplus::basename_part(std::to_string(view->id_view)), "desc");
      

      //If features or descriptors file are missing, compute them
      if (bForce || !stlplus::file_exists(sFeat) || !stlplus::file_exists(sDesc))
      {
        if (!ReadImage(sView_filename.c_str(), &imageGray))
        {
          std::cout << "Error: can't read image: " << sView_filename << std::endl;
          continue;
        }

        // Compute features and descriptors and export them to files
        std::cout << "Extracting features from image " << view->id_view << " - (" << sFeat << ")" << std::endl;
        std::unique_ptr<Regions> regions;

        image_describer->Describe(imageGray, regions);
        image_describer->Save(regions.get(), sFeat, sDesc);
      }
    }
    std::cout << "Task done in (s): " << timer.elapsed() << std::endl;
  }
  return EXIT_SUCCESS;
}
