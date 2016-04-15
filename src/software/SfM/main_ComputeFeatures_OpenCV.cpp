
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
    const Image<unsigned char> * mask = nullptr)
  {
    cv::Mat img;
    cv::eigen2cv(image.GetMat(), img);

    cv::Mat m_mask;
    if(mask != nullptr) {
      cv::eigen2cv(mask->GetMat(), m_mask);
    }

    std::vector< cv::KeyPoint > vec_keypoints;
    cv::Mat m_desc;

    cv::Ptr<cv::Feature2D> extractor = cv::AKAZE::create(cv::AKAZE::DESCRIPTOR_KAZE);
    extractor->detectAndCompute(img, m_mask, vec_keypoints, m_desc);

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

        memcpy(descriptor.data(),
               m_desc.ptr<typename DescriptorT::bin_type>(cpt),
               DescriptorT::static_size*sizeof(DescriptorT::bin_type));
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
///
//- Create an Image_describer interface that use and OpenCV extraction method
// i.e. with the SIFT detector+descriptor
// Regions is the same as classic SIFT : 128 unsigned char
class SIFT_OPENCV_Image_describer : public Image_describer
{
public:
  SIFT_OPENCV_Image_describer() : Image_describer() {}

  ~SIFT_OPENCV_Image_describer() {}

  bool Set_configuration_preset(EDESCRIBER_PRESET preset){
    return true;
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
    const image::Image<unsigned char> * mask = nullptr)
  {
    // Convert for opencv
    cv::Mat img;
    cv::eigen2cv(image.GetMat(), img);

    // Convert mask image into cv::Mat
    cv::Mat m_mask;
    if(mask != nullptr) {
      cv::eigen2cv(mask->GetMat(), m_mask);
    }

    // Create a SIFT detector
    std::vector< cv::KeyPoint > v_keypoints;
    cv::Mat m_desc;
    cv::Ptr<cv::Feature2D> siftdetector = cv::xfeatures2d::SIFT::create();

    // Process SIFT computation
    siftdetector->detectAndCompute(img, m_mask, v_keypoints, m_desc);

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
  }
private:
  int _i;
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

  // required
  cmd.add( make_option('i', sSfM_Data_Filename, "input_file") );
  cmd.add( make_option('o', sOutDir, "outdir") );
  // Optional
  cmd.add( make_option('f', bForce, "force") );
#ifdef USE_OCVSIFT
  cmd.add( make_option('m', sImage_Describer_Method, "describerMethod") );
#endif

  try {
      if (argc == 1) throw std::string("Invalid command line parameter.");
      cmd.process(argc, argv);
  } catch(const std::string& s) {
      std::cerr << "Usage: " << argv[0] << '\n'
      << "[-i|--input_file]: a SfM_Data file \n"
      << "[-o|--outdir] path \n"
      << "\n[Optional]\n"
      << "[-f|--force] Force to recompute data\n"
#ifdef USE_OCVSIFT
      << "[-m|--describerMethod]\n"
      << "  (method to use to describe an image):\n"
      << "   AKAZE_OPENCV (default),\n"
      << "   SIFT_OPENCV: SIFT FROM OPENCV\n"
#endif
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
            << "--force " << bForce << std::endl;

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
    Image<unsigned char> imageGray, globalMask, imageMask;

    const std::string sGlobalMask_filename = stlplus::create_filespec(sOutDir, "mask.png");
    if(stlplus::file_exists(sGlobalMask_filename))
      ReadImage(sGlobalMask_filename.c_str(), &globalMask);

    C_Progress_display my_progress_bar( sfm_data.GetViews().size(),
      std::cout, "\n- EXTRACT FEATURES -\n" );
    for(Views::const_iterator iterViews = sfm_data.views.begin();
        iterViews != sfm_data.views.end();
        ++iterViews, ++my_progress_bar)
    {
      const View * view = iterViews->second.get();
      const std::string
        sView_filename = stlplus::create_filespec(sfm_data.s_root_path, view->s_Img_path),
        sFeat = stlplus::create_filespec(sOutDir, stlplus::basename_part(sView_filename), "feat"),
        sDesc = stlplus::create_filespec(sOutDir, stlplus::basename_part(sView_filename), "desc");

      //If features or descriptors file are missing, compute them
      if (bForce || !stlplus::file_exists(sFeat) || !stlplus::file_exists(sDesc))
      {
        if (!ReadImage(sView_filename.c_str(), &imageGray))
          continue;

        Image<unsigned char> * mask = nullptr; // The mask is null by default

        const std::string sImageMask_filename =
          stlplus::create_filespec(sfm_data.s_root_path,
            stlplus::basename_part(sView_filename) + "_mask", "png");

        if(stlplus::file_exists(sImageMask_filename))
          ReadImage(sImageMask_filename.c_str(), &imageMask);

        // The mask point to the globalMask, if a valid one exists for the current image
        if(globalMask.Width() == imageGray.Width() && globalMask.Height() == imageGray.Height())
          mask = &globalMask;
        // The mask point to the imageMask (individual mask) if a valid one exists for the current image
        if(imageMask.Width() == imageGray.Width() && imageMask.Height() == imageGray.Height())
          mask = &imageMask;

        // Compute features and descriptors and export them to files
        std::unique_ptr<Regions> regions;
        image_describer->Describe(imageGray, regions, mask);
        image_describer->Save(regions.get(), sFeat, sDesc);
      }
    }
    std::cout << "Task done in (s): " << timer.elapsed() << std::endl;
  }
  return EXIT_SUCCESS;
}
