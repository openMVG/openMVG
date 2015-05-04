
// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


#include "openMVG/image/image.hpp"
#include "openMVG/sfm/sfm.hpp"
#include "openMVG/features/features.hpp"
#include <cereal/archives/json.hpp>

/// Generic Image Collection image matching
#include "openMVG/matching_image_collection/Matcher_Regions_AllInMemory.hpp"
#include "openMVG/matching_image_collection/GeometricFilter.hpp"
#include "openMVG/matching_image_collection/F_ACRobust.hpp"
#include "openMVG/matching_image_collection/E_ACRobust.hpp"
#include "openMVG/matching_image_collection/H_ACRobust.hpp"
#include "software/SfM/pairwiseAdjacencyDisplay.hpp"
#include "openMVG/matching/indMatch_utils.hpp"
#include "openMVG/system/timer.hpp"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"
#include "third_party/progress/progress.hpp"

#include <opencv2/opencv.hpp>
#include "opencv2/core/eigen.hpp"

#include <cstdlib>
#include <fstream>

using namespace openMVG;
using namespace openMVG::matching;
using namespace openMVG::robust;
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
using namespace openMVG::features;
// Reuse the existing AKAZE floating point Keypoint.
typedef features::AKAZE_Float_Regions AKAZE_OpenCV_Regions;
// Define the Interface
class AKAZE_OCV_Image_describer : public Image_describer
{
public:
  AKAZE_OCV_Image_describer():Image_describer(){}

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
  std::string sGeometricModel = "f";
  float fDistRatio = .6f;
  int iMatchingVideoMode = -1;
  std::string sPredefinedPairList = "";

  cmd.add( make_option('i', sSfM_Data_Filename, "input_file") );
  cmd.add( make_option('o', sOutDir, "outdir") );
  cmd.add( make_option('r', fDistRatio, "ratio") );
  cmd.add( make_option('g', sGeometricModel, "geometricModel") );
  cmd.add( make_option('v', iMatchingVideoMode, "videoModeMatching") );
  cmd.add( make_option('l', sPredefinedPairList, "pairList") );

  try {
      if (argc == 1) throw std::string("Invalid command line parameter.");
      cmd.process(argc, argv);
  } catch(const std::string& s) {
      std::cerr << "Usage: " << argv[0] << '\n'
      << "[-i|--input_file]: a SfM_Data file \n"
      << "[-o|--outdir path] \n"
      << "\n[Optional]\n"
      << "[-r|--ratio Distance ratio to discard non meaningful matches\n"
      << "   0.6 typical value (you can use 0.8 to have more matches)]"
      << "[-g|--geometricModel\n"
      << "  (pairwise correspondences filtering thanks to robust model estimation):\n"
      << "   f: fundamental matrix,\n"
      << "   e: essential matrix,\n"
      << "   h: homography matrix]\n"
      << "[-v|--videoModeMatching\n"
      << "  (sequence matching with an overlap of X images)\n"
      << "   X: with match 0 with (1->X), ...]\n"
      << "   2: will match 0 with (1,2), 1 with (2,3), ...\n"
      << "   3: will match 0 with (1,2,3), 1 with (2,3,4), ...]\n"
      << "[-l]--pairList file"
      << std::endl;

      std::cerr << s << std::endl;
      return EXIT_FAILURE;
  }

  std::cout << " You called : " <<std::endl
            << argv[0] << std::endl
            << "--input_file " << sSfM_Data_Filename << std::endl
            << "--outdir " << sOutDir << std::endl
            << "--ratio " << fDistRatio << std::endl
            << "--geometricModel " << sGeometricModel << std::endl
            << "--videoModeMatching " << iMatchingVideoMode << std::endl;

  ePairMode ePairmode = (iMatchingVideoMode == -1 ) ? PAIR_EXHAUSTIVE : PAIR_CONTIGUOUS;

  if (sPredefinedPairList.length()) {
    std::cout << "--pairList " << sPredefinedPairList << std::endl;
    ePairmode = PAIR_FROM_FILE;
    if (iMatchingVideoMode>0) {
      std::cerr << "\nIncompatible options: --videoModeMatching and --pairList" << std::endl;
      return EXIT_FAILURE;
    }
  }

  if (sOutDir.empty())  {
    std::cerr << "\nIt is an invalid output directory" << std::endl;
    return EXIT_FAILURE;
  }

  eGeometricModel eGeometricModelToCompute = FUNDAMENTAL_MATRIX;
  std::string sGeometricMatchesFilename = "";
  switch(sGeometricModel[0])
  {
    case 'f': case 'F':
      eGeometricModelToCompute = FUNDAMENTAL_MATRIX;
      sGeometricMatchesFilename = "matches.f.txt";
    break;
    case 'e': case 'E':
      eGeometricModelToCompute = ESSENTIAL_MATRIX;
      sGeometricMatchesFilename = "matches.e.txt";
    break;
    case 'h': case 'H':
      eGeometricModelToCompute = HOMOGRAPHY_MATRIX;
      sGeometricMatchesFilename = "matches.h.txt";
    break;
    default:
      std::cerr << "Unknown geometric model" << std::endl;
      return EXIT_FAILURE;
  }

  // -----------------------------
  // a. Load input scene
  // b. Compute features and descriptors
  // c. Compute putative descriptor matches
  // d. Geometric filtering of putative matches
  // e. Export some statistics
  // -----------------------------

  // Create output dir
  if (!stlplus::folder_exists(sOutDir))
    if (!stlplus::folder_create(sOutDir))
    {
      std::cerr << "Cannot create output directory" << std::endl;
      return EXIT_FAILURE;
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

  //---------------------------------------
  // b. Compute features and descriptor
  //    - extract sift features and descriptor
  //    - if keypoints already computed, re-load them
  //    - else save features and descriptors on disk
  //---------------------------------------

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
    image_describer.reset(new AKAZE_OCV_Image_describer);

    // Export the used Image_describer to a file for future regions loading
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


  {
    Timer timer;
    std::cout << "\n\n - EXTRACT FEATURES - " << std::endl;

    Image<unsigned char> imageGray;
    C_Progress_display my_progress_bar( sfm_data.getViews().size() );
    for(Views::const_iterator iterViews = sfm_data.views.begin();
        iterViews != sfm_data.views.end();
        ++iterViews, ++my_progress_bar)
    {
      const View * view = iterViews->second.get();
      const std::string sView_filename = stlplus::create_filespec(sfm_data.s_root_path,
        view->s_Img_path);
      const std::string sFeat = stlplus::create_filespec(sOutDir,
        stlplus::basename_part(sView_filename), "feat");
      const std::string sDesc = stlplus::create_filespec(sOutDir,
        stlplus::basename_part(sView_filename), "desc");

      //If features or descriptors file are missing, compute them
      if (!stlplus::file_exists(sFeat) || !stlplus::file_exists(sDesc))
      {
        if (!ReadImage(sView_filename.c_str(), &imageGray))
          continue;

        // Compute features and descriptors and export them to files
        std::unique_ptr<Regions> regions;
        image_describer->Describe(imageGray, regions);
        image_describer->Save(regions.get(), sFeat, sDesc);
      }
    }
    std::cout << "Task done in (s): " << timer.elapsed() << std::endl;
  }

  //---------------------------------------
  // c. Compute putative descriptor matches
  //    - Descriptor matching (according user method choice)
  //    - Keep correspondences only if NearestNeighbor ratio is ok
  //---------------------------------------
  PairWiseMatches map_PutativesMatches;

  // List views as a vector of filenames & imagesizes (alias)
  std::vector<std::string> vec_fileNames;
  vec_fileNames.reserve(sfm_data.getViews().size());
  std::vector<std::pair<size_t, size_t> > vec_imagesSize;
  vec_imagesSize.reserve(sfm_data.getViews().size());
  for (Views::const_iterator iter = sfm_data.getViews().begin();
    iter != sfm_data.getViews().end();
    ++iter)
  {
    const View * v = iter->second.get();
    vec_fileNames.push_back(stlplus::create_filespec(sfm_data.s_root_path,
        v->s_Img_path));
    vec_imagesSize.push_back( std::make_pair( v->ui_width, v->ui_height) );
  }



  std::cout << std::endl << " - PUTATIVE MATCHES - " << std::endl;
  // If the matches already exists, reload them
  if (stlplus::file_exists(sOutDir + "/matches.putative.txt"))
  {
    PairedIndMatchImport(sOutDir + "/matches.putative.txt", map_PutativesMatches);
    std::cout << "\t PREVIOUS RESULTS LOADED" << std::endl;
  }
  else // Compute the putative matches
  {
    std::cout << "Use: ";
    switch (ePairmode)
    {
      case PAIR_EXHAUSTIVE: std::cout << "exhaustive pairwise matching" << std::endl; break;
      case PAIR_CONTIGUOUS: std::cout << "sequence pairwise matching" << std::endl; break;
      case PAIR_FROM_FILE:  std::cout << "user defined pairwise matching" << std::endl; break;
    }

    Timer timer;
    // TODO: make matching method a parameter
    Matcher_Regions_AllInMemory collectionMatcher(fDistRatio, ANN_L2);
    //Matcher_Regions_AllInMemory collectionMatcher(fDistRatio, BRUTE_FORCE_HAMMING);
    if (collectionMatcher.loadData(*image_describer.get(), vec_fileNames, sOutDir))
    {
      // Get pair to match according the matching mode:
      Pair_Set pairs;
      switch (ePairmode)
      {
        case PAIR_EXHAUSTIVE: pairs = exhaustivePairs(sfm_data.getViews().size()); break;
        case PAIR_CONTIGUOUS: pairs = contiguousWithOverlap(sfm_data.getViews().size(), iMatchingVideoMode); break;
        case PAIR_FROM_FILE:
          if(!loadPairs(sfm_data.getViews().size(), sPredefinedPairList, pairs))
          {
              return EXIT_FAILURE;
          };
          break;
      }
      // Photometric matching of putative pairs
      collectionMatcher.Match(vec_fileNames, pairs, map_PutativesMatches);
      //---------------------------------------
      //-- Export putative matches
      //---------------------------------------
      std::ofstream file (std::string(sOutDir + "/matches.putative.txt").c_str());
      if (file.is_open())
        PairedIndMatchToStream(map_PutativesMatches, file);
      file.close();
    }
    std::cout << "Task done in (s): " << timer.elapsed() << std::endl;
  }
  //-- export putative matches Adjacency matrix
  PairWiseMatchingToAdjacencyMatrixSVG(vec_fileNames.size(),
    map_PutativesMatches,
    stlplus::create_filespec(sOutDir, "PutativeAdjacencyMatrix", "svg"));

  //---------------------------------------
  // d. Geometric filtering of putative matches
  //    - AContrario Estimation of the desired geometric model
  //    - Use an upper bound for the a contrario estimated threshold
  //---------------------------------------

  // Prepare the features and matches provider
  std::shared_ptr<Features_Provider> feats_provider = std::make_shared<Features_Provider>();
  if (!feats_provider->load(sfm_data, sOutDir, image_describer)) {
    std::cerr << std::endl << "Invalid features." << std::endl;
    return EXIT_FAILURE;
  }

  PairWiseMatches map_GeometricMatches;

  ImageCollectionGeometricFilter collectionGeomFilter(feats_provider.get());
  const double maxResidualError = 4.0;
  {
    Timer timer;
    std::cout << std::endl << " - GEOMETRIC FILTERING - " << std::endl;
    switch (eGeometricModelToCompute)
    {
      case FUNDAMENTAL_MATRIX:
      {
       collectionGeomFilter.Filter(
          GeometricFilter_FMatrix_AC(maxResidualError),
          map_PutativesMatches,
          map_GeometricMatches,
          vec_imagesSize);
      }
      break;
      case ESSENTIAL_MATRIX:
      {
        // Build the intrinsic parameter map for each view
        std::map<IndexT, Mat3> map_K;
        size_t cpt = 0;
        for (Views::const_iterator iter = sfm_data.getViews().begin();
          iter != sfm_data.getViews().end();
          ++iter, ++cpt)
        {
          const View * v = iter->second.get();
          if (sfm_data.getIntrinsics().count(v->id_intrinsic))
          {
            const IntrinsicBase * ptrIntrinsic = sfm_data.getIntrinsics().find(v->id_intrinsic)->second.get();
            switch (ptrIntrinsic->getType())
            {
              case PINHOLE_CAMERA:
              case PINHOLE_CAMERA_RADIAL1:
              case PINHOLE_CAMERA_RADIAL3:
                const Pinhole_Intrinsic * ptrPinhole = (const Pinhole_Intrinsic*)(ptrIntrinsic);
                map_K[cpt] = ptrPinhole->K();
              break;
            }
          }
        }

        collectionGeomFilter.Filter(
          GeometricFilter_EMatrix_AC(map_K, maxResidualError),
          map_PutativesMatches,
          map_GeometricMatches,
          vec_imagesSize);

        //-- Perform an additional check to remove pairs with poor overlap
        std::vector<PairWiseMatches::key_type> vec_toRemove;
        for (PairWiseMatches::const_iterator iterMap = map_GeometricMatches.begin();
          iterMap != map_GeometricMatches.end(); ++iterMap)
        {
          const size_t putativePhotometricCount = map_PutativesMatches.find(iterMap->first)->second.size();
          const size_t putativeGeometricCount = iterMap->second.size();
          const float ratio = putativeGeometricCount / (float)putativePhotometricCount;
          if (putativeGeometricCount < 50 || ratio < .3f)  {
            // the pair will be removed
            vec_toRemove.push_back(iterMap->first);
          }
        }
        //-- remove discarded pairs
        for (std::vector<PairWiseMatches::key_type>::const_iterator
          iter =  vec_toRemove.begin(); iter != vec_toRemove.end(); ++iter)
        {
          map_GeometricMatches.erase(*iter);
        }
      }
      break;
      case HOMOGRAPHY_MATRIX:
      {
        collectionGeomFilter.Filter(
          GeometricFilter_HMatrix_AC(maxResidualError),
          map_PutativesMatches,
          map_GeometricMatches,
          vec_imagesSize);
      }
      break;
    }

    //---------------------------------------
    //-- Export geometric filtered matches
    //---------------------------------------
    std::ofstream file (string(sOutDir + "/" + sGeometricMatchesFilename).c_str());
    if (file.is_open())
      PairedIndMatchToStream(map_GeometricMatches, file);
    file.close();

    std::cout << "Task done in (s): " << timer.elapsed() << std::endl;

    //-- export Adjacency matrix
    std::cout << "\n Export Adjacency Matrix of the pairwise's geometric matches"
      << std::endl;
    PairWiseMatchingToAdjacencyMatrixSVG(vec_fileNames.size(),
      map_GeometricMatches,
      stlplus::create_filespec(sOutDir, "GeometricAdjacencyMatrix", "svg"));
  }
  return EXIT_SUCCESS;
}
