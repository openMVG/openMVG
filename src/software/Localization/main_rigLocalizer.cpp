#include <openMVG/localization/VoctreeLocalizer.hpp>
#if HAVE_CCTAG
#include <openMVG/localization/CCTagLocalizer.hpp>
#endif
#include <openMVG/rig/Rig.hpp>
#include <openMVG/image/image_io.hpp>
#include <openMVG/dataio/FeedProvider.hpp>
#include <openMVG/features/image_describer.hpp>
#include <openMVG/robust_estimation/robust_estimators.hpp>
#include <openMVG/logger.hpp>

#include <boost/filesystem.hpp>
#include <boost/progress.hpp>
#include <boost/program_options.hpp> 
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/min.hpp>
#include <boost/accumulators/statistics/max.hpp>
#include <boost/accumulators/statistics/sum.hpp>

#include <iostream>
#include <string>
#include <vector>
#include <chrono>
#include <memory>

#if HAVE_ALEMBIC
#include <openMVG/sfm/AlembicExporter.hpp>
#endif // HAVE_ALEMBIC


namespace bfs = boost::filesystem;
namespace bacc = boost::accumulators;
namespace po = boost::program_options;

using namespace openMVG;

enum DescriberType
{
  SIFT
#if HAVE_CCTAG
  ,CCTAG,
  SIFT_CCTAG
#endif
};

inline DescriberType stringToDescriberType(const std::string& describerType)
{
  if(describerType == "SIFT")
    return DescriberType::SIFT;
#if HAVE_CCTAG
  if (describerType == "CCTAG")
    return DescriberType::CCTAG;
  if(describerType == "SIFT_CCTAG")
    return DescriberType::SIFT_CCTAG;
#endif
  throw std::invalid_argument("Unsupported describer type "+describerType);
}

inline std::string describerTypeToString(DescriberType describerType)
{
  if(describerType == DescriberType::SIFT)
    return "SIFT";
#if HAVE_CCTAG
  if (describerType == DescriberType::CCTAG)
    return "CCTAG";
  if(describerType == DescriberType::SIFT_CCTAG)
    return "SIFT_CCTAG";
#endif
  throw std::invalid_argument("Unrecognized DescriberType "+std::to_string(describerType));
}

std::ostream& operator<<(std::ostream& os, const DescriberType describerType)
{
  os << describerTypeToString(describerType);
  return os;
}

std::istream& operator>>(std::istream &in, DescriberType &describerType)
{
  std::string token;
  in >> token;
  describerType = stringToDescriberType(token);
  return in;
}

std::string myToString(std::size_t i, std::size_t zeroPadding)
{
  std::stringstream ss;
  ss << std::setw(zeroPadding) << std::setfill('0') << i;
  return ss.str();
}

/**
 * @brief It checks if the value for the reprojection error or the matching error
 * is compatible with the given robust estimator. The value cannot be 0 for 
 * LORansac, for ACRansac a value of 0 means to use infinity (ie estimate the 
 * threshold during ransac process)
 * @param e The estimator to be checked.
 * @param value The value for the reprojection or matching error.
 * @return true if the value is compatible
 */
bool checkRobustEstimator(robust::EROBUST_ESTIMATOR e, double &value)
{
  if(e != robust::EROBUST_ESTIMATOR::ROBUST_ESTIMATOR_LORANSAC &&
     e != robust::EROBUST_ESTIMATOR::ROBUST_ESTIMATOR_ACRANSAC)
  {
    POPART_CERR("Only " << robust::EROBUST_ESTIMATOR::ROBUST_ESTIMATOR_ACRANSAC 
            << " and " << robust::EROBUST_ESTIMATOR::ROBUST_ESTIMATOR_LORANSAC 
            << " are supported.");
    return false;
  }
  if(value == 0 && 
     e == robust::EROBUST_ESTIMATOR::ROBUST_ESTIMATOR_ACRANSAC)
  {
    // for acransac set it to infinity
    value = std::numeric_limits<double>::infinity();
  }
  // for loransac we need thresholds > 0
  if(e == robust::EROBUST_ESTIMATOR::ROBUST_ESTIMATOR_LORANSAC)
  {
    const double minThreshold = 1e-6;
    if(value <= minThreshold)
    {
      POPART_CERR("Error: errorMax and matchingError cannot be 0 with " 
              << robust::EROBUST_ESTIMATOR::ROBUST_ESTIMATOR_LORANSAC 
              << " estimator.");
      return false;     
    }
  }

  return true;
}

int main(int argc, char** argv)
{
  // common parameters
  std::string sfmFilePath;                //< the OpenMVG .json data file
  std::string descriptorsFolder;          //< the the folder containing the descriptors
  std::string mediaPath;                  //< the media file to localize
  std::string filelist;                  //< the media file to localize
  std::string rigCalibPath;               //< the file containing the calibration data for the file (subposes)
  //< the preset for the feature extractor
  features::EDESCRIBER_PRESET featurePreset = features::EDESCRIBER_PRESET::NORMAL_PRESET;     
  //< the type of features to use for localization
  DescriberType descriptorType = DescriberType::SIFT;        
  //< the estimator to use for resection
  robust::EROBUST_ESTIMATOR resectionEstimator = robust::EROBUST_ESTIMATOR::ROBUST_ESTIMATOR_ACRANSAC;        
  //< the estimator to use for matching
  robust::EROBUST_ESTIMATOR matchingEstimator = robust::EROBUST_ESTIMATOR::ROBUST_ESTIMATOR_ACRANSAC;        
  //< the possible choices for the estimators as strings
  const std::string str_estimatorChoices = ""+robust::EROBUST_ESTIMATOR_enumToString(robust::EROBUST_ESTIMATOR::ROBUST_ESTIMATOR_ACRANSAC)
                                          +","+robust::EROBUST_ESTIMATOR_enumToString(robust::EROBUST_ESTIMATOR::ROBUST_ESTIMATOR_LORANSAC);
  bool refineIntrinsics = false;
  bool useLocalizeRigNaive = false;
  double resectionErrorMax = 4.0;                    //< the maximum error allowed for resection
  double matchingErrorMax = 4.0;               //< the maximum error allowed for image matching with geometric validation


  // parameters for voctree localizer
  std::string vocTreeFilepath;            //< the vocabulary tree file
  std::string weightsFilepath;            //< the vocabulary tree weights file
  std::string algostring = "AllResults";   //< the localization algorithm to use for the voctree localizer
  std::size_t numResults = 4;              //< number of documents to search when querying the voctree
  std::size_t maxResults = 10;             //< maximum number of matching documents to retain
  // parameters for cctag localizer
  std::size_t nNearestKeyFrames = 5;           //

#if HAVE_ALEMBIC
  std::string exportFile = "trackedcameras.abc"; //!< the export file
#endif
  
  std::size_t numCameras = 3;
  po::options_description desc("This program is used to localize a camera rig composed of internally calibrated cameras");
  desc.add_options()
      ("help,h", "Print this message")
      ("descriptors", po::value<DescriberType>(&descriptorType)->default_value(descriptorType), 
        "Type of descriptors to use {SIFT"
#ifdef HAVE_CCTAG
        ", CCTAG, SIFT_CCTAG"
#endif
        "}")
      ("preset", po::value<features::EDESCRIBER_PRESET>(&featurePreset)->default_value(featurePreset), 
          "Preset for the feature extractor when localizing a new image "
          "{LOW,MEDIUM,NORMAL,HIGH,ULTRA}")
      ("resectionEstimator", po::value<robust::EROBUST_ESTIMATOR>(&resectionEstimator)->default_value(resectionEstimator), 
          std::string("The type of *sac framework to use for resection "
          "{"+str_estimatorChoices+"}").c_str())
      ("matchingEstimator", po::value<robust::EROBUST_ESTIMATOR>(&matchingEstimator)->default_value(matchingEstimator), 
          std::string("The type of *sac framework to use for matching "
          "{"+str_estimatorChoices+"}").c_str())
      ("sfmdata", po::value<std::string>(&sfmFilePath)->required(),
          "The sfm_data.json kind of file generated by OpenMVG.")
      ("descriptorPath", po::value<std::string>(&descriptorsFolder)->required(),
          "Folder containing the .desc.")
      ("mediapath", po::value<std::string>(&mediaPath)->required(),
          "The folder path containing all the synchronised image subfolders "
          "assocated to each camera")
      ("filelist", po::value<std::string>(&filelist),
          "An optional txt file containing the images to use for calibration. "
          "This file must have the same name in each camera folder and contains "
          "the list of images to load.")
      ("refineIntrinsics", po::bool_switch(&refineIntrinsics),
          "Enable/Disable camera intrinsics refinement for each localized image")
      ("nCameras", po::value<size_t>(&numCameras)->default_value(numCameras), 
          "Number of cameras composing the rig")
      ("calibration", po::value<std::string>(&rigCalibPath)->required(), 
          "The file containing the calibration data for the rig (subposes)")
      ("reprojectionError", po::value<double>(&resectionErrorMax)->default_value(resectionErrorMax), 
          "Maximum reprojection error (in pixels) allowed for resectioning. If set "
          "to 0 it lets the ACRansac select an optimal value.")
      ("useLocalizeRigNaive", po::bool_switch(&useLocalizeRigNaive),
          "Enable/Disable the naive method for rig localization: naive method tries "
          "to localize each camera separately. This is enabled by default if the "
          "library has not been built with openGV.")
  // parameters for voctree localizer
      ("voctree", po::value<std::string>(&vocTreeFilepath),
          "[voctree] Filename for the vocabulary tree")
      ("voctreeWeights", po::value<std::string>(&weightsFilepath),
          "[voctree] Filename for the vocabulary tree weights")
      ("algorithm", po::value<std::string>(&algostring)->default_value(algostring),
          "[voctree] Algorithm type: {FirstBest,BestResult,AllResults,Cluster}" )
      ("nbImageMatch", po::value<size_t>(&numResults)->default_value(numResults),
          "[voctree] Number of images to retrieve in the database")
      ("maxResults", po::value<size_t>(&maxResults)->default_value(maxResults), 
          "[voctree] For algorithm AllResults, it stops the image matching when "
          "this number of matched images is reached. If 0 it is ignored.")
      ("matchingError", po::value<double>(&matchingErrorMax)->default_value(matchingErrorMax), 
          "[voctree] Maximum matching error (in pixels) allowed for image matching with "
          "geometric verification. If set to 0 it lets the ACRansac select "
          "an optimal value.")
#if HAVE_CCTAG
  // parameters for cctag localizer
      ("nNearestKeyFrames", po::value<size_t>(&nNearestKeyFrames)->default_value(nNearestKeyFrames),
        "[cctag] Number of images to retrieve in database")
#endif
#if HAVE_ALEMBIC
      ("output", po::value<std::string>(&exportFile)->default_value(exportFile),
        "Filename for the SfM_Data export file (where camera poses will be stored)."
        " Default : trackedcameras.abc.")
#endif
          ;

  po::variables_map vm;

  try
  {
    po::store(po::parse_command_line(argc, argv, desc), vm);

    if(vm.count("help") || (argc == 1))
    {
      POPART_COUT(desc);
      return EXIT_SUCCESS;
    }

    po::notify(vm);
  }
  catch(boost::program_options::required_option& e)
  {
    POPART_CERR("ERROR: " << e.what() << std::endl);
    POPART_COUT("Usage:\n\n" << desc);
    return EXIT_FAILURE;
  }
  catch(boost::program_options::error& e)
  {
    POPART_CERR("ERROR: " << e.what() << std::endl);
    POPART_COUT("Usage:\n\n" << desc);
    return EXIT_FAILURE;
  }

  if(!checkRobustEstimator(matchingEstimator, matchingErrorMax) || 
     !checkRobustEstimator(resectionEstimator, resectionErrorMax))
  {
    return EXIT_FAILURE;
  }

  // just debugging prints, print out all the parameters
  {
    POPART_COUT("Program called with the following parameters:");
    POPART_COUT("\tsfmdata: " << sfmFilePath);
    POPART_COUT("\tpreset: " << featurePreset);
    POPART_COUT("\tmediapath: " << mediaPath);
    POPART_COUT("\tresectionEstimator: " << resectionEstimator);
    POPART_COUT("\tmatchingEstimator: " << matchingEstimator);
    POPART_COUT("\tdescriptorPath: " << descriptorsFolder);
    POPART_COUT("\trefineIntrinsics: " << refineIntrinsics);
    POPART_COUT("\tuseLocalizeRigNaive: " << useLocalizeRigNaive);
    POPART_COUT("\rreprojectionError: " << resectionErrorMax);
    POPART_COUT("\tnCameras: " << numCameras);
    if(!filelist.empty())
      POPART_COUT("\tfilelist: " << filelist);
    POPART_COUT("\tdescriptors: " << descriptorType);
    if((DescriberType::SIFT==descriptorType)
#if HAVE_CCTAG
            ||(DescriberType::SIFT_CCTAG==descriptorType)
#endif
      )
    {
      // parameters for voctree localizer
      POPART_COUT("\tvoctree: " << vocTreeFilepath);
      POPART_COUT("\tweights: " << weightsFilepath);
      POPART_COUT("\tnbImageMatch: " << numResults);
      POPART_COUT("\tmaxResults: " << maxResults);
      POPART_COUT("\talgorithm: " << algostring);
      POPART_COUT("\tmatchingError: " << matchingErrorMax);
    }
#if HAVE_CCTAG
    else
    {
      POPART_COUT("\tnNearestKeyFrames: " << nNearestKeyFrames);
    }
#endif

  }

  std::unique_ptr<localization::LocalizerParameters> param;
  
  std::unique_ptr<localization::ILocalizer> localizer;
  
  // initialize the localizer according to the chosen type of describer
  if((DescriberType::SIFT==descriptorType)
#if HAVE_CCTAG
            ||(DescriberType::SIFT_CCTAG==descriptorType)
#endif
      )
  {
    POPART_COUT("Localizing sequence using the voctree localizer");
    localization::VoctreeLocalizer* tmpLoc = new localization::VoctreeLocalizer(sfmFilePath,
                                                            descriptorsFolder,
                                                            vocTreeFilepath,
                                                            weightsFilepath
#if HAVE_CCTAG
                                                            , DescriberType::SIFT_CCTAG==descriptorType
#endif
                                                            );
    localizer.reset(tmpLoc);
    
    localization::VoctreeLocalizer::Parameters *tmpParam = new localization::VoctreeLocalizer::Parameters();
    param.reset(tmpParam);
    tmpParam->_algorithm = localization::VoctreeLocalizer::initFromString(algostring);;
    tmpParam->_numResults = numResults;
    tmpParam->_maxResults = maxResults;
    tmpParam->_ccTagUseCuda = false;
    tmpParam->_matchingError = matchingErrorMax;
    
  }
#if HAVE_CCTAG
  else
  {
    localization::CCTagLocalizer* tmpLoc = new localization::CCTagLocalizer(sfmFilePath, descriptorsFolder);
    localizer.reset(tmpLoc);
    
    localization::CCTagLocalizer::Parameters *tmpParam = new localization::CCTagLocalizer::Parameters();
    param.reset(tmpParam);
    tmpParam->_nNearestKeyFrames = nNearestKeyFrames;
  }
#endif 

  assert(localizer);
  assert(param);
  
  // set other common parameters
  param->_featurePreset = featurePreset;
  param->_refineIntrinsics = refineIntrinsics;
  param->_errorMax = resectionErrorMax;
  param->_resectionEstimator = resectionEstimator;
  param->_matchingEstimator = matchingEstimator;
  param->_useLocalizeRigNaive = useLocalizeRigNaive;

  if(!localizer->isInit())
  {
    POPART_CERR("ERROR while initializing the localizer!");
    return EXIT_FAILURE;
  }

#if HAVE_ALEMBIC
  dataio::AlembicExporter exporter(exportFile);
  exporter.initAnimatedCamera("rig");
  exporter.addPoints(localizer->getSfMData().GetLandmarks());
  
  std::vector<dataio::AlembicExporter> cameraExporters;
  cameraExporters.reserve(numCameras);
  // this contains the full path and the root name of the file without the extension
  const std::string basename = (bfs::path(exportFile).parent_path() / bfs::path(exportFile).stem()).string();
  for(std::size_t i = 0; i < numCameras; ++i)
  {
    cameraExporters.emplace_back(basename+".cam"+myToString(i, 2)+".abc");
    cameraExporters.back().initAnimatedCamera("cam"+myToString(i, 2));
  }
#endif

  vector<dataio::FeedProvider*> feeders(numCameras);

  // Init the feeder for each camera
  for(std::size_t idCamera = 0; idCamera < numCameras; ++idCamera)
  {
    const std::string subMediaFilepath = mediaPath + "/" + std::to_string(idCamera);
    const std::string calibFile = subMediaFilepath + "/intrinsics.txt";
    const std::string feedPath = subMediaFilepath + "/"+filelist;

    // create the feedProvider
    feeders[idCamera] = new dataio::FeedProvider(feedPath, calibFile);
    if(!feeders[idCamera]->isInit())
    {
      POPART_CERR("ERROR while initializing the FeedProvider for the camera " 
              << idCamera << " " << feedPath);
      return EXIT_FAILURE;
    }
  }
  
  bool haveNext = true;
  std::size_t frameCounter = 0;
  std::size_t numLocalizedFrames = 0;
  
  // load the subposes
  std::vector<geometry::Pose3> vec_subPoses;
  rig::loadRigCalibration(rigCalibPath, vec_subPoses);
  assert(vec_subPoses.size() == numCameras-1);
  geometry::Pose3 rigPose;
  
  // Define an accumulator set for computing the mean and the
  // standard deviation of the time taken for localization
  bacc::accumulator_set<double, bacc::stats<bacc::tag::mean, bacc::tag::min, bacc::tag::max, bacc::tag::sum > > stats;

  // store the result
  std::vector< std::vector<localization::LocalizationResult> > rigResultPerFrame;
  
  while(haveNext)
  {
    // @fixme It's better to have arrays of pointers...
    std::vector<image::Image<unsigned char> > vec_imageGrey;
    std::vector<cameras::Pinhole_Intrinsic_Radial_K3 > vec_queryIntrinsics;
    vec_imageGrey.reserve(numCameras);
    vec_queryIntrinsics.reserve(numCameras);
           
    // for each camera get the image and the associated internal parameters
    for(int idCamera = 0; idCamera < numCameras; ++idCamera)
    {
      image::Image<unsigned char> imageGrey;
      cameras::Pinhole_Intrinsic_Radial_K3 queryIntrinsics;
      bool hasIntrinsics = false;
      std::string currentImgName;
      haveNext = feeders[idCamera]->next(imageGrey, queryIntrinsics, currentImgName, hasIntrinsics);
      
      if(!haveNext)
      {
        if(idCamera > 0)
        {
          // this is quite odd, it means that eg the fist camera has an image but
          // one of the others has not image
          POPART_CERR("This is weird... Camera " << idCamera << " seems not to have any available images while some other cameras do...");
          return EXIT_FAILURE;  // a bit harsh but if we are here it's cheesy to say the less
        }
        break;
      }
      
      // for now let's suppose that the cameras are calibrated internally too
      if(!hasIntrinsics)
      {
        POPART_CERR("For now only internally calibrated cameras are supported!"
                << "\nCamera " << idCamera << " does not have calibration for image " << currentImgName);
        return EXIT_FAILURE;  // a bit harsh but if we are here it's cheesy to say the less
      }
      
      vec_imageGrey.push_back(imageGrey);
      vec_queryIntrinsics.push_back(queryIntrinsics);
    }
    
    if(!haveNext)
    {
      // no more images are available
      break;
    }
    
    POPART_COUT("******************************");
    POPART_COUT("FRAME " << myToString(frameCounter, 4));
    POPART_COUT("******************************");
    auto detect_start = std::chrono::steady_clock::now();
    std::vector<localization::LocalizationResult> localizationResults;
    const bool isLocalized = localizer->localizeRig(vec_imageGrey,
                                                    param.get(),
                                                    vec_queryIntrinsics,
                                                    vec_subPoses,
                                                    rigPose,
                                                    localizationResults);
    auto detect_end = std::chrono::steady_clock::now();
    auto detect_elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(detect_end - detect_start);
    POPART_COUT("\nLocalization took  " << detect_elapsed.count() << " [ms]");
    stats(detect_elapsed.count());
    
    rigResultPerFrame.push_back(localizationResults);
    
    if(isLocalized)
    {
      ++numLocalizedFrames;
#if HAVE_ALEMBIC
      // save the position of the main camera
      exporter.addCameraKeyframe(rigPose, &vec_queryIntrinsics[0], mediaPath, frameCounter, frameCounter);
      assert(cameraExporters.size()==numCameras);
      assert(localizationResults.size()==numCameras);
      assert(vec_queryIntrinsics.size()==numCameras);
      // save the position of all cameras of the rig
      for(std::size_t camIDX = 0; camIDX < numCameras; ++camIDX)
      {
        POPART_COUT("cam pose" << camIDX << "\n" <<  localizationResults[camIDX].getPose().rotation() << "\n" << localizationResults[camIDX].getPose().center());
        if(camIDX > 0)
          POPART_COUT("cam subpose" << camIDX-1 << "\n" <<  vec_subPoses[camIDX-1].rotation() << "\n" << vec_subPoses[camIDX-1].center());
        cameraExporters[camIDX].addCameraKeyframe(localizationResults[camIDX].getPose(), &vec_queryIntrinsics[camIDX], mediaPath, frameCounter, frameCounter);
      }
#endif
    }
    else
    {
     POPART_CERR("Unable to localize frame " << frameCounter);
#if HAVE_ALEMBIC
      exporter.jumpKeyframe();
      assert(cameraExporters.size()==numCameras);
      for(std::size_t camIDX = 0; camIDX < numCameras; ++camIDX)
      {
        cameraExporters[camIDX].jumpKeyframe();
      }
#endif
    }

    ++frameCounter;
  }
  
  // print out some time stats
  POPART_COUT("\n\n******************************");
  POPART_COUT("Localized " << numLocalizedFrames << " / " << frameCounter << " images");
  POPART_COUT("Processing took " << bacc::sum(stats) / 1000 << " [s] overall");
  POPART_COUT("Mean time for localization:   " << bacc::mean(stats) << " [ms]");
  POPART_COUT("Max time for localization:   " << bacc::max(stats) << " [ms]");
  POPART_COUT("Min time for localization:   " << bacc::min(stats) << " [ms]");
}
