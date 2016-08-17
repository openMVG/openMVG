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
  int i;
  in >> i;
  describerType = static_cast<DescriberType>(i);
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
  std::string filelist;                   //< the media file to localize
  std::string outputFile;                 //< the name of the file where to store the calibration data
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
  po::options_description desc("This program is used to calibrate a camera rig composed of internally calibrated cameras");
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
      ("reprojectionError", po::value<double>(&resectionErrorMax)->default_value(resectionErrorMax), 
          "Maximum reprojection error (in pixels) allowed for resectioning. If set "
          "to 0 it lets the ACRansac select an optimal value.")
      ("outfile,o", po::value<std::string>(&outputFile)->required(),
          "The name of the file where to store the calibration data")
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
      ("export,e", po::value<std::string>(&exportFile)->default_value(exportFile),
          "Filename for the alembic file containing the rig poses with the 3D points. "
          "It also saves a file for each camera named 'filename.cam##.abc'.")
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
    POPART_COUT("\reprojectionError: " << resectionErrorMax);
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
    POPART_COUT("Calibrating sequence using the voctree localizer");
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

  if(!localizer->isInit())
  {
    POPART_CERR("ERROR while initializing the localizer!");
    return EXIT_FAILURE;
  }

#if HAVE_ALEMBIC
  dataio::AlembicExporter exporter(exportFile);
  exporter.addPoints(localizer->getSfMData().GetLandmarks());
#endif

  // Create a camera rig
  rig::Rig rig;

  // Loop over all cameras of the rig
  for(std::size_t idCamera = 0; idCamera < numCameras; ++idCamera)
  {
    POPART_COUT("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
    POPART_COUT("CAMERA " << idCamera);
    POPART_COUT("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n");
    const std::string subMediaFilepath = mediaPath + "/" + std::to_string(idCamera);
    const std::string calibFile = subMediaFilepath + "/intrinsics.txt";
    const std::string feedPath = subMediaFilepath + "/"+filelist;

    // create the feedProvider
    dataio::FeedProvider feed(feedPath, calibFile);
    if(!feed.isInit())
    {
      POPART_CERR("ERROR while initializing the FeedProvider!");
      return EXIT_FAILURE;
    }

    //std::string featureFile, cameraResultFile, pointsFile;
    //featureFile = subMediaFilepath + "/cctag" + std::to_string(nRings) + "CC.out";
    //cameraResultFile = inputFolder + "/" + std::to_string(i) + "/cameras.txt";
    //std::ofstream result;
    //result.open(cameraResultFile);
    //pointsFile = inputFolder + "/points.txt";

    image::Image<unsigned char> imageGrey;
    cameras::Pinhole_Intrinsic_Radial_K3 queryIntrinsics;
    bool hasIntrinsics = false;

    size_t frameCounter = 0;
    std::string currentImgName;

    // Define an accumulator set for computing the mean and the
    // standard deviation of the time taken for localization
    bacc::accumulator_set<double, bacc::stats<bacc::tag::mean, bacc::tag::min, bacc::tag::max, bacc::tag::sum > > stats;

    // used to collect the match data result
    std::vector<localization::LocalizationResult> vLocalizationResults;
    while(feed.next(imageGrey, queryIntrinsics, currentImgName, hasIntrinsics))
    {
      POPART_COUT("******************************");
      POPART_COUT("FRAME " << myToString(frameCounter, 4));
      POPART_COUT("******************************");
      auto detect_start = std::chrono::steady_clock::now();
      localization::LocalizationResult localizationResult;
      const bool ok = localizer->localize(imageGrey,
                                          param.get(),
                                          hasIntrinsics/*useInputIntrinsics*/,
                                          queryIntrinsics,
                                          localizationResult);
      assert( ok == localizationResult.isValid() );
      vLocalizationResults.emplace_back(localizationResult);
      auto detect_end = std::chrono::steady_clock::now();
      auto detect_elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(detect_end - detect_start);
      POPART_COUT("\nLocalization took  " << detect_elapsed.count() << " [ms]");
      stats(detect_elapsed.count());
      
#if HAVE_ALEMBIC
      if(localizationResult.isValid())
      {
        exporter.appendCamera("camera"+std::to_string(idCamera)+"."+myToString(frameCounter,4), localizationResult.getPose(), &queryIntrinsics, subMediaFilepath, frameCounter, frameCounter);
      }
      else
      {
        // @fixme for now just add a fake camera so that it still can be see in MAYA
        exporter.appendCamera("camera"+std::to_string(idCamera)+".V."+myToString(frameCounter,4), geometry::Pose3(), &queryIntrinsics, subMediaFilepath, frameCounter, frameCounter);
      }
#endif
      ++frameCounter;
    }

    rig.setTrackingResult(vLocalizationResults, idCamera);

    // print out some time stats
    POPART_COUT("\n\n******************************");
    POPART_COUT("Processed " << frameCounter << " images for camera " << idCamera);
    POPART_COUT("Processing took " << bacc::sum(stats) / 1000 << " [s] overall");
    POPART_COUT("Mean time for localization:   " << bacc::mean(stats) << " [ms]");
    POPART_COUT("Max time for localization:   " << bacc::max(stats) << " [ms]");
    POPART_COUT("Min time for localization:   " << bacc::min(stats) << " [ms]");
  }
  
  {
    // just for statistics purposes
    const std::size_t numRigCam = rig.nCams();
    POPART_COUT("\n\n******************************");
    for(std::size_t idCam = 0; idCam < numRigCam; ++idCam)
    {
      auto & currResult = rig.getLocalizationResults(idCam);
      std::size_t numLocalized = 0;
      for(const auto &curr : currResult)
      {
        if(curr.isValid())
          ++numLocalized;
      }
      POPART_COUT("Camera " << idCam << " localized " 
              << numLocalized << "/" << currResult.size());
    }
    
  }
  
  POPART_COUT("Rig calibration initialization...");
  if(!rig.initializeCalibration())
  {
    POPART_CERR("Unable to find a proper initialization for the relative poses! Aborting...");
    return EXIT_FAILURE;
  }
  POPART_COUT("Rig calibration optimization...");
  if(!rig.optimizeCalibration())
  {
    POPART_CERR("Unable to optimize the relative poses! Aborting...");
    return EXIT_FAILURE;
  }
  
  // save the rig calibration (subposes)
  rig.saveCalibration(outputFile);
  
  
  // just print out the results
  // the first rig pose
  if(rig.getPosesSize() > 0)
  {
    POPART_COUT("First pose of the rig");
    const geometry::Pose3 &pose = rig.getPose(0); 
    POPART_COUT("R\n" << pose.rotation());
    POPART_COUT("center\n" << pose.center());
    POPART_COUT("t\n" << pose.translation());
  }
  
  // get the subposes of the cameras inside the rig
  const std::vector<geometry::Pose3>& subposes = rig.getRelativePoses();
  assert(numCameras-1 == subposes.size());
  for(std::size_t i = 0; i < subposes.size(); ++i)
  {
    const geometry::Pose3 &pose = subposes[i];
    POPART_COUT("--------------------");
    POPART_COUT("Subpose p0" << i+1); // from camera 0 to camera i+1
    POPART_COUT("R\n" << pose.rotation());
    POPART_COUT("center\n" << pose.center());
    POPART_COUT("t\n" << pose.translation());
    POPART_COUT("--------------------\n");
  }
  
  return EXIT_SUCCESS;
}
