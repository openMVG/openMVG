#include <openMVG/localization/ILocalizer.hpp>
#include <openMVG/localization/VoctreeLocalizer.hpp>
#if HAVE_CCTAG
#include <openMVG/localization/CCTagLocalizer.hpp>
#endif
#include <openMVG/localization/LocalizationResult.hpp>
#include <openMVG/localization/optimization.hpp>
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
    if( value <= minThreshold)
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
  std::string calibFile;                    //< the calibration file
  std::string sfmFilePath;                  //< the OpenMVG .json data file
  std::string descriptorsFolder;            //< the OpenMVG .json data file
  std::string mediaFilepath;                //< the media file to localize
  //< the preset for the feature extractor
  features::EDESCRIBER_PRESET featurePreset = features::EDESCRIBER_PRESET::NORMAL_PRESET;     
  //< the preset for the feature extractor
  DescriberType descriptorType = DescriberType::SIFT;        
  //< the estimator to use for resection
  robust::EROBUST_ESTIMATOR resectionEstimator = robust::EROBUST_ESTIMATOR::ROBUST_ESTIMATOR_ACRANSAC;        
  //< the estimator to use for matching
  robust::EROBUST_ESTIMATOR matchingEstimator = robust::EROBUST_ESTIMATOR::ROBUST_ESTIMATOR_ACRANSAC;        
  //< the possible choices for the estimators as strings
  const std::string str_estimatorChoices = ""+robust::EROBUST_ESTIMATOR_enumToString(robust::EROBUST_ESTIMATOR::ROBUST_ESTIMATOR_ACRANSAC)
                                          +","+robust::EROBUST_ESTIMATOR_enumToString(robust::EROBUST_ESTIMATOR::ROBUST_ESTIMATOR_LORANSAC);
  bool refineIntrinsics = false;
  double resectionErrorMax = 4.0;  //< the maximum error allowed for resection
  double matchingErrorMax = 4.0;   //< the maximum error allowed for image matching with geometric validation
  
  // voctree parameters
  std::string algostring = "AllResults";
  std::size_t numResults = 4;       //< number of documents to search when querying the voctree
  std::size_t maxResults = 10;      //< maximum number of matching documents to retain
  std::size_t numCommonViews = 3;
  std::string vocTreeFilepath;      //< the vocabulary tree file
  std::string weightsFilepath;      //< the vocabulary tree weights file
  
#if HAVE_ALEMBIC
  std::string exportFile = "trackedcameras.abc"; //!< the export file
#else
  std::string exportFile = "localizationResult.json"; //!< the export file
#endif
#if HAVE_CCTAG
  // parameters for cctag localizer
  std::size_t nNearestKeyFrames = 5;   
#endif
  bool globalBundle = false;              ///< If !refineIntrinsics it can run a final global budndle to refine the scene
  bool noDistortion = false;              ///< It does not count the distortion
  bool noBArefineIntrinsics = false;      ///< It does not refine intrinsics during BA
  std::size_t minPointVisibility = 0;
  
  std::string visualDebug = "";        ///< whether to save visual debug info
  bool useVoctreeLocalizer = true;        ///< whether to use the voctreeLocalizer or cctagLocalizer
  bool useSIFT_CCTAG = false;             ///< whether to use SIFT_CCTAG

  po::options_description desc(
      "This program takes as input a media (image, image sequence, video) and a database (voctree, 3D structure data) \n"
      "and returns for each frame a pose estimation for the camera.");
  
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
      ("calibration", po::value<std::string>(&calibFile)/*->required( )*/, 
          "Calibration file")
      ("sfmdata", po::value<std::string>(&sfmFilePath)->required(), 
          "The sfm_data.json kind of file generated by OpenMVG.")
      ("descriptorPath", po::value<std::string>(&descriptorsFolder)->required(), 
          "Folder containing the .desc.")
      ("mediafile", po::value<std::string>(&mediaFilepath)->required(), 
          "The folder path or the filename for the media to track")
      ("refineIntrinsics", po::bool_switch(&refineIntrinsics), 
          "Enable/Disable camera intrinsics refinement for each localized image")
      ("reprojectionError", po::value<double>(&resectionErrorMax)->default_value(resectionErrorMax), 
          "Maximum reprojection error (in pixels) allowed for resectioning. If set "
          "to 0 it lets the ACRansac select an optimal value.")
// voctree specific options
      ("nbImageMatch", po::value<std::size_t>(&numResults)->default_value(numResults), 
          "[voctree] Number of images to retrieve in database")
      ("maxResults", po::value<std::size_t>(&maxResults)->default_value(maxResults), 
          "[voctree] For algorithm AllResults, it stops the image matching when "
          "this number of matched images is reached. If 0 it is ignored.")
      ("commonviews", po::value<std::size_t>(&numCommonViews)->default_value(numCommonViews), 
          "[voctree] Number of minimum images in which a point must be seen to "
          "be used in cluster tracking")
      ("voctree", po::value<std::string>(&vocTreeFilepath)->required(), 
          "[voctree] Filename for the vocabulary tree")
      ("voctreeWeights", po::value<std::string>(&weightsFilepath), 
          "[voctree] Filename for the vocabulary tree weights")
      ("algorithm", po::value<std::string>(&algostring)->default_value(algostring), 
          "[voctree] Algorithm type: FirstBest, BestResult, AllResults, Cluster" )
      ("matchingError", po::value<double>(&matchingErrorMax)->default_value(matchingErrorMax), 
          "[voctree] Maximum matching error (in pixels) allowed for image matching with "
          "geometric verification. If set to 0 it lets the ACRansac select "
          "an optimal value.")
// cctag specific options
#if HAVE_CCTAG
      ("nNearestKeyFrames", po::value<size_t>(&nNearestKeyFrames)->default_value(nNearestKeyFrames), 
          "[cctag] Number of images to retrieve in the database")
#endif
// final bundle adjustment options
      ("globalBundle", po::bool_switch(&globalBundle), 
          "[bundle adjustment] If --refineIntrinsics is not set, this option "
          "allows to run a final global budndle adjustment to refine the scene")
      ("noDistortion", po::bool_switch(&noDistortion), 
          "[bundle adjustment] It does not take into account distortion during "
          "the BA, it consider the distortion coefficients all equal to 0")
      ("noBArefineIntrinsics", po::bool_switch(&noBArefineIntrinsics), 
          "[bundle adjustment] It does not refine intrinsics during BA")
      ("minPointVisibility", po::value<size_t>(&minPointVisibility)->default_value(minPointVisibility), 
          "[bundle adjustment] Minimum number of observation that a point must "
          "have in order to be considered for bundle adjustment")
      ("visualDebug", po::value<std::string>(&visualDebug), 
          "If a directory is provided it enables visual debug and saves all the "
          "debugging info in that directory")
#if HAVE_ALEMBIC
      ("output", po::value<std::string>(&exportFile)->default_value(exportFile), 
          "Filename for the SfM_Data export file (where camera poses will be stored). "
          "Default : trackedcameras.abc. It will also save the localization "
          "results (raw data) as .json with the same name")
#else
      ("output", po::value<std::string>(&exportFile)->default_value(exportFile), 
          "Filename for the SfM_Data export file containing the localization "
          "results. Default : localizationResult.json.")
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
  
  // just for debugging purpose, print out all the parameters
  {
    // decide the localizer to use based on the type of feature
    useVoctreeLocalizer = ((DescriberType::SIFT==descriptorType)
#if HAVE_CCTAG
            || (DescriberType::SIFT_CCTAG==descriptorType)
#endif
            );
    
#if HAVE_CCTAG   
    // check whether we have to use SIFT and CCTAG together
    useSIFT_CCTAG = (DescriberType::SIFT_CCTAG==descriptorType);
#endif    
    
    // the bundle adjustment can be run for now only if the refine intrinsics option is not set
    globalBundle = (globalBundle && !refineIntrinsics);
    POPART_COUT("Program called with the following parameters:");
    POPART_COUT("\tsfmdata: " << sfmFilePath);
    POPART_COUT("\tdescriptors: " << descriptorType);
    POPART_COUT("\tpreset: " << featurePreset);
    POPART_COUT("\tresectionEstimator: " << resectionEstimator);
    POPART_COUT("\tmatchingEstimator: " << matchingEstimator);
    POPART_COUT("\tcalibration: " << calibFile);
    POPART_COUT("\tdescriptorPath: " << descriptorsFolder);
    POPART_COUT("\trefineIntrinsics: " << refineIntrinsics);
    POPART_COUT("\treprojectionError: " << resectionErrorMax);
    POPART_COUT("\tmediafile: " << mediaFilepath);
    if(useVoctreeLocalizer)
    {
      POPART_COUT("\tvoctree: " << vocTreeFilepath);
      POPART_COUT("\tweights: " << weightsFilepath);
      POPART_COUT("\tnbImageMatch: " << numResults);
      POPART_COUT("\tmaxResults: " << maxResults);
      POPART_COUT("\tcommon views: " << numCommonViews);
      POPART_COUT("\talgorithm: " << algostring);
      POPART_COUT("\tmatchingError: " << matchingErrorMax);
    }
#if HAVE_CCTAG 
    else
    {
      POPART_COUT("\tnNearestKeyFrames: " << nNearestKeyFrames);
    }
#endif 
    POPART_COUT("\tminPointVisibility: " << minPointVisibility);
    POPART_COUT("\tglobalBundle: " << globalBundle);
    POPART_COUT("\tnoDistortion: " << noDistortion);
    POPART_COUT("\tnoBArefineIntrinsics: " << noBArefineIntrinsics);
    POPART_COUT("\tvisualDebug: " << visualDebug);
  }

  // if the provided directory for visual debugging does not exist create it
  // recursively
  if((!visualDebug.empty()) && (!bfs::exists(visualDebug)))
  {
    bfs::create_directories(visualDebug);
  }
 
  // this contains the full path and the root name of the file without the extension
  const std::string basename = (bfs::path(exportFile).parent_path() / bfs::path(exportFile).stem()).string();
  
  
  //***********************************************************************
  // Localizer initialization
  //***********************************************************************
  
  std::unique_ptr<localization::LocalizerParameters> param;
  
  std::unique_ptr<localization::ILocalizer> localizer;
  
  // initialize the localizer according to the chosen type of describer
  if((DescriberType::SIFT==descriptorType)
#if HAVE_CCTAG
            ||(DescriberType::SIFT_CCTAG==descriptorType)
#endif
      )
  {
    localization::VoctreeLocalizer* tmpLoc = new localization::VoctreeLocalizer(sfmFilePath,
                                                   descriptorsFolder,
                                                   vocTreeFilepath,
                                                   weightsFilepath
#ifdef HAVE_CCTAG
                                                    , useSIFT_CCTAG
#endif
                                                  );
    localizer.reset(tmpLoc);
    
    localization::VoctreeLocalizer::Parameters *tmpParam = new localization::VoctreeLocalizer::Parameters();
    param.reset(tmpParam);
    tmpParam->_algorithm = localization::VoctreeLocalizer::initFromString(algostring);;
    tmpParam->_numResults = numResults;
    tmpParam->_maxResults = maxResults;
    tmpParam->_numCommonViews = numCommonViews;
    tmpParam->_ccTagUseCuda = false;
    tmpParam->_matchingError = matchingErrorMax;
  }
#if HAVE_CCTAG
  else
  {
    localization::CCTagLocalizer* tmpLoc = new localization::CCTagLocalizer(sfmFilePath, descriptorsFolder);
    localizer.reset(tmpLoc);
    
    localization::CCTagLocalizer::Parameters* tmpParam = new localization::CCTagLocalizer::Parameters();
    param.reset(tmpParam);
    tmpParam->_nNearestKeyFrames = nNearestKeyFrames;
  }
#endif 
   
  assert(localizer);
  assert(param);
  
  // set other common parameters
  param->_featurePreset = featurePreset;
  param->_refineIntrinsics = refineIntrinsics;
  param->_visualDebug = visualDebug;
  param->_errorMax = resectionErrorMax;
  param->_resectionEstimator = resectionEstimator;
  param->_matchingEstimator = matchingEstimator;
  
  
  if(!localizer->isInit())
  {
    POPART_CERR("ERROR while initializing the localizer!");
    return EXIT_FAILURE;
  }
  
  // create the feedProvider
  dataio::FeedProvider feed(mediaFilepath, calibFile);
  if(!feed.isInit())
  {
    POPART_CERR("ERROR while initializing the FeedProvider!");
    return EXIT_FAILURE;
  }
  
#if HAVE_ALEMBIC
  // init alembic exporter
  dataio::AlembicExporter exporter( exportFile );
  exporter.addPoints(localizer->getSfMData().GetLandmarks());
  exporter.initAnimatedCamera("camera");
#endif
  
  image::Image<unsigned char> imageGrey;
  cameras::Pinhole_Intrinsic_Radial_K3 queryIntrinsics;
  bool hasIntrinsics = false;
  
  std::size_t frameCounter = 0;
  std::size_t goodFrameCounter = 0;
  std::vector<std::string> goodFrameList;
  std::string currentImgName;
  
  //***********************************************************************
  // Main loop
  //***********************************************************************
  
  // Define an accumulator set for computing the mean and the
  // standard deviation of the time taken for localization
  bacc::accumulator_set<double, bacc::stats<bacc::tag::mean, bacc::tag::min, bacc::tag::max, bacc::tag::sum > > stats;
  
  std::vector<localization::LocalizationResult> vec_localizationResults;
  
  while(feed.readImage(imageGrey, queryIntrinsics, currentImgName, hasIntrinsics))
  {
    POPART_COUT("******************************");
    POPART_COUT("FRAME " << myToString(frameCounter,4));
    POPART_COUT("******************************");
    localization::LocalizationResult localizationResult;
    auto detect_start = std::chrono::steady_clock::now();
    localizer->localize(imageGrey, 
                       param.get(),
                       hasIntrinsics /*useInputIntrinsics*/,
                       queryIntrinsics,
                       localizationResult,
                       currentImgName);
    auto detect_end = std::chrono::steady_clock::now();
    auto detect_elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(detect_end - detect_start);
    POPART_COUT("\nLocalization took  " << detect_elapsed.count() << " [ms]");
    stats(detect_elapsed.count());
    
    vec_localizationResults.emplace_back(localizationResult);

    // save data
    if(localizationResult.isValid())
    {
#if HAVE_ALEMBIC
      exporter.addCameraKeyframe(localizationResult.getPose(), &queryIntrinsics, currentImgName, frameCounter, frameCounter);
#endif
      
      goodFrameCounter++;
      goodFrameList.push_back(currentImgName + " : " + std::to_string(localizationResult.getIndMatch3D2D().size()) );
    }
    else
    {
      POPART_CERR("Unable to localize frame " << frameCounter);
#if HAVE_ALEMBIC
      exporter.jumpKeyframe(currentImgName);
#endif
    }
    ++frameCounter;
    feed.goToNextFrame();
  }
  localization::save(vec_localizationResults, basename+".locres.json");
  
  
  //***********************************************************************
  // Global bundle
  //***********************************************************************
  
  if(globalBundle)
  {
    POPART_COUT("\n\n\n***********************************************");
    POPART_COUT("Bundle Adjustment - Refining the whole sequence");
    POPART_COUT("***********************************************\n\n");
    // run a bundle adjustment
    const bool b_allTheSame = true;
    const bool b_refineStructure = false;
    const bool b_refinePose = true;
    const bool BAresult = localization::refineSequence(vec_localizationResults,
                                                       b_allTheSame, 
                                                       !noBArefineIntrinsics, 
                                                       noDistortion, 
                                                       b_refinePose,
                                                       b_refineStructure,
                                                       basename+".sfmdata.BUNDLE",
                                                       minPointVisibility);
    if(!BAresult)
    {
      POPART_CERR("Bundle Adjustment failed!");
    }
    else
    {
#if HAVE_ALEMBIC
      // now copy back in a new abc with the same name file and BUNDLE appended at the end
      dataio::AlembicExporter exporterBA( basename+".BUNDLE.abc" );
      exporterBA.initAnimatedCamera("camera");
      std::size_t idx = 0;
      for(const localization::LocalizationResult &res : vec_localizationResults)
      {
        if(res.isValid())
        {
          assert(idx < vec_localizationResults.size());
          exporterBA.addCameraKeyframe(res.getPose(), &res.getIntrinsics(), currentImgName, frameCounter, frameCounter);
        }
        else
        {
          exporterBA.jumpKeyframe(currentImgName);
        }
        idx++;
      }
      exporterBA.addPoints(localizer->getSfMData().GetLandmarks());
#endif
      localization::save(vec_localizationResults, basename+".locres.BUNDLE.json");
    }
  }
  
  // print out some time stats
  POPART_COUT("\n\n******************************");
  POPART_COUT("Localized " << goodFrameCounter << "/" << frameCounter << " images");
  POPART_COUT("Images localized with the number of 2D/3D matches during localization :");
  for(std::size_t i = 0; i < goodFrameList.size(); ++i)
    POPART_COUT(goodFrameList[i]);
  POPART_COUT("Processing took " << bacc::sum(stats)/1000 << " [s] overall");
  POPART_COUT("Mean time for localization:   " << bacc::mean(stats) << " [ms]");
  POPART_COUT("Max time for localization:   " << bacc::max(stats) << " [ms]");
  POPART_COUT("Min time for localization:   " << bacc::min(stats) << " [ms]");
}
