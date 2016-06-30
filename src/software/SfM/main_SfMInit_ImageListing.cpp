// Copyright (c) 2012, 2013, 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
#include "openMVG/exif/exif_IO_EasyExif.hpp"

#include "openMVG/exif/sensor_width_database/ParseDatabase.hpp"

#include "openMVG/image/image.hpp"
#include "openMVG/stl/split.hpp"

#include "openMVG/sfm/sfm.hpp"
#include "openMVG/sfm/sfm_view_metadata.hpp"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <memory>
#include <string>
#include <vector>

using namespace openMVG;
using namespace openMVG::cameras;
using namespace openMVG::exif;
using namespace openMVG::image;
using namespace openMVG::sfm;

/// Check that Kmatrix is a string like "f;0;ppx;0;f;ppy;0;0;1"
/// With f,ppx,ppy as valid numerical value
bool checkIntrinsicStringValidity(const std::string & Kmatrix, double & focal, double & ppx, double & ppy)
{
  std::vector<std::string> vec_str;
  stl::split(Kmatrix, ";", vec_str);
  if (vec_str.size() != 9)  {
    std::cerr << "\n Missing ';' character" << std::endl;
    return false;
  }
  // Check that all K matrix value are valid numbers
  for (size_t i = 0; i < vec_str.size(); ++i) {
    double readvalue = 0.0;
    std::stringstream ss;
    ss.str(vec_str[i]);
    if (! (ss >> readvalue) )  {
      std::cerr << "\n Used an invalid not a number character" << std::endl;
      return false;
    }
    if (i==0) focal = readvalue;
    if (i==2) ppx = readvalue;
    if (i==5) ppy = readvalue;
  }
  return true;
}

/// Recursively list all files from a folder with a specific extension
void listFiles(const std::string& folderOrFile, const std::vector<std::string>& extensions, std::vector<std::string>& outFiles)
{
  if(stlplus::is_file(folderOrFile))
  {
    std::string fileExtension = stlplus::extension_part(folderOrFile);
    std::transform(fileExtension.begin(), fileExtension.end(), fileExtension.begin(), ::tolower);
    for(const std::string& extension: extensions)
    {
      if(fileExtension == extension)
      {
        outFiles.push_back(folderOrFile);
        return;
      }
    }
  }
  else if(stlplus::is_folder(folderOrFile))
  {
    // list all files of the folder
    const std::vector<std::string> allFiles = stlplus::folder_all(folderOrFile);
    for(const std::string& item: allFiles)
    {
      const std::string itemPath = stlplus::create_filespec(folderOrFile, item);
      listFiles(itemPath, extensions, outFiles);
    }
  }
}

/// Retrieve resources path from a json file
/// Need a "resource" variable in the json
bool retrieveResources(const std::string& jsonFile, std::vector<std::string>& vec_imagePaths)
{
  if(!stlplus::file_exists(jsonFile))
  {
    std::cerr << "File \"" << jsonFile << "\" does not exists." << std::endl;
    return false;
  }

  // Read file
  std::ifstream jsonStream(jsonFile, std::ifstream::binary);
  
  if(!jsonStream.is_open())
    throw std::runtime_error("Unable to open "+jsonFile);

  // get length of file:
  jsonStream.seekg (0, jsonStream.end);
  const int length = jsonStream.tellg();
  jsonStream.seekg (0, jsonStream.beg);
  // read data as a block:
  std::string jsonString;
  jsonString.resize(length);
  jsonStream.read(&jsonString[0], length);
  jsonStream.close();

  // Parse json
  rapidjson::Document document;
  document.Parse<0>(&jsonString[0]);
  if(!document.IsObject())
  {
    std::cerr << "File \"" << jsonFile << "\" is not in json format." << std::endl;
    return false;
  }
  if(!document.HasMember("resources"))
  {
    std::cerr << "No member \"resources\" in json file" << std::endl;
    return false;
  }
  rapidjson::Value& resourcesValue = document["resources"];
  for(rapidjson::Value::ConstValueIterator it = resourcesValue.Begin(); it != resourcesValue.End(); it++)
    vec_imagePaths.push_back(it->GetString());

  return true;

}

//
// Create the description of an input image dataset for OpenMVG toolsuite
// - Export a SfM_Data file with View & Intrinsic data
//
int main(int argc, char **argv)
{
  CmdLine cmd;

  std::string sImageDir;
  std::string sJsonFile;
  std::string sfileDatabase;
  std::string sOutputDir;
  std::string sKmatrix;

  std::string i_User_camera_model;

  int b_Group_camera_model = 0;
  bool b_use_UID = false;
  bool b_storeMetadata = false;

  double userFocalLengthPixel = -1.0;
  double userSensorWidth = -1.0;

  cmd.add( make_option('i', sImageDir, "imageDirectory") );
  cmd.add( make_option('j', sJsonFile, "jsonFile") );
  cmd.add( make_option('d', sfileDatabase, "sensorWidthDatabase") );
  cmd.add( make_option('o', sOutputDir, "outputDirectory") );
  cmd.add( make_option('f', userFocalLengthPixel, "focal") );
  cmd.add( make_option('s', userSensorWidth, "sensorWidth") );
  cmd.add( make_option('k', sKmatrix, "intrinsics") );
  cmd.add( make_option('c', i_User_camera_model, "camera_model") );
  cmd.add( make_option('g', b_Group_camera_model, "group_camera_model") );
  cmd.add( make_option('u', b_use_UID, "use_UID") );
  cmd.add( make_option('m', b_storeMetadata, "storeMetadata") );

  try {
      if (argc == 1) throw std::string("Invalid command line parameter.");
      cmd.process(argc, argv);
  } catch(const std::string& s)
  {
      std::cerr << "Usage: " << argv[0] << '\n'
      << "[-i|--imageDirectory]\n"
      << "[-j|--jsonFile] Input file with all the user options. It can be used to provide a list of images instead of a directory.\n"
      << "[-d|--sensorWidthDatabase]\n"
      << "[-o|--outputDirectory]\n"
      << "[-f|--focal] (pixels)\n"
      << "[-s|--sensorWidth] (mm)\n"
      << "[-k|--intrinsics] Kmatrix: \"f;0;ppx;0;f;ppy;0;0;1\"\n"
      << "[-c|--camera_model] Camera model type (pinhole, radial1, radial3, brown or fisheye4)\n"
      << "[-g|--group_camera_model]\n"
      << "\t 0-> each view have its own camera intrinsic parameters,\n"
      << "\t 1-> (default) view share camera intrinsic parameters based on metadata, if no metadata each view has its own camera intrinsic parameters\n"
      << "\t 2-> view share camera intrinsic parameters based on metadata, if no metadata they are grouped by folder\n"
      << "[-u|--use_UID] Generate a UID (unique identifier) for each view. By default, the key is the image index.\n"
      << "[-m|--storeMetadata] Store image metadata in the sfm data\n"
      << std::endl;

      std::cerr << s << std::endl;
      return EXIT_FAILURE;
  }

  std::cout << " You called : " <<std::endl
            << argv[0] << std::endl
            << "--imageDirectory " << sImageDir << std::endl
            << "--jsonFile " << sJsonFile << std::endl
            << "--sensorWidthDatabase " << sfileDatabase << std::endl
            << "--outputDirectory " << sOutputDir << std::endl
            << "--focal " << userFocalLengthPixel << std::endl
            << "--sensorWidth " << userSensorWidth << std::endl
            << "--intrinsics " << sKmatrix << std::endl
            << "--camera_model " << i_User_camera_model << std::endl
            << "--group_camera_model " << b_Group_camera_model << std::endl
            << "--use_UID " << b_use_UID << std::endl
            << "--storeMetadata " << b_storeMetadata << std::endl;

  EINTRINSIC e_User_camera_model = PINHOLE_CAMERA_START;
  if(!i_User_camera_model.empty())
    e_User_camera_model = EINTRINSIC_stringToEnum(i_User_camera_model);
  double ppx = -1.0, ppy = -1.0;

  if(!sImageDir.empty() && !sJsonFile.empty())
  {
    std::cerr << "\nCannot combine -i and -j options" << std::endl;
    return EXIT_FAILURE;
  }
  if ( !sImageDir.empty() && !stlplus::folder_exists( sImageDir ) )
  {
    std::cerr << "\nThe input directory doesn't exist" << std::endl;
    return EXIT_FAILURE;
  }
  if (sOutputDir.empty())
  {
    std::cerr << "\nInvalid output directory" << std::endl;
    return EXIT_FAILURE;
  }

  if ( !stlplus::folder_exists( sOutputDir ) )
  {
    if ( !stlplus::folder_create( sOutputDir ))
    {
      std::cerr << "\nCannot create output directory" << std::endl;
      return EXIT_FAILURE;
    }
  }

  if (sKmatrix.size() > 0 && userFocalLengthPixel != -1.0)
  {
    std::cerr << "\nCannot combine -f and -k options" << std::endl;
    return EXIT_FAILURE;
  }

  if (sKmatrix.size() > 0 &&
    !checkIntrinsicStringValidity(sKmatrix, userFocalLengthPixel, ppx, ppy) )
  {
    std::cerr << "\nInvalid K matrix input" << std::endl;
    return EXIT_FAILURE;
  }

  std::vector<Datasheet> vec_database;
  if (!sfileDatabase.empty())
  {
    if ( !parseDatabase( sfileDatabase, vec_database ) )
    {
      std::cerr
       << "\nInvalid input database: " << sfileDatabase
       << ", please specify a valid file." << std::endl;
      return EXIT_FAILURE;
    }
  }

  // Retrieve image paths
  std::vector<std::string> vec_images;
  const std::vector<std::string> supportedExtensions{ "jpg", "jpeg" };

  if (!sJsonFile.empty())
  {
    // Retrieve resources from json file
    std::vector<std::string> vec_resources;
    if (!retrieveResources(sJsonFile, vec_resources))
      return EXIT_FAILURE;
    // Retrieve images from resources
    for(const std::string& item: vec_resources)
    {
      listFiles(item, supportedExtensions, vec_images);
    }
    std::sort(vec_images.begin(), vec_images.end());
  }
  if(!sImageDir.empty())
  {
    vec_images = stlplus::folder_files( sImageDir );
  }
  std::sort(vec_images.begin(), vec_images.end());

  // Configure an empty scene with Views and their corresponding cameras
  SfM_Data sfm_data;
  sfm_data.s_root_path = "";
  if(!sImageDir.empty())
    sfm_data.s_root_path = sImageDir; // Setup main image root_path
  Views & views = sfm_data.views;
  Intrinsics & intrinsics = sfm_data.intrinsics;

  C_Progress_display my_progress_bar( vec_images.size(),
      std::cout, "\n- Image listing -\n" );
  std::ostringstream error_report_stream;
  for ( std::vector<std::string>::const_iterator iter_image = vec_images.begin();
    iter_image != vec_images.end();
    ++iter_image, ++my_progress_bar )
  {
    // Read meta data to fill camera parameter (w,h,focal,ppx,ppy) fields.
    double width = -1.0, height = -1.0, focalPix = -1.0;
    ppx = ppy = -1.0;
    
    const std::string imageFilename = *iter_image;
    std::string imageAbsFilepath;
    if(!sJsonFile.empty())
      imageAbsFilepath = imageFilename;
    else if(!sImageDir.empty())
      imageAbsFilepath = stlplus::create_filespec( sImageDir, imageFilename );
    std::string imageFolder = stlplus::folder_part(imageAbsFilepath);

    // Test if the image format is supported:
    if (openMVG::image::GetFormat(imageAbsFilepath.c_str()) == openMVG::image::Unknown)
    {
      error_report_stream
          << stlplus::filename_part(imageAbsFilepath) << ": Unknown image file format." << "\n";
      continue; // image cannot be opened
    }

    ImageHeader imgHeader;
    if (!openMVG::image::ReadImageHeader(imageAbsFilepath.c_str(), &imgHeader))
      continue; // image cannot be read

    width = imgHeader.width;
    height = imgHeader.height;
    
    if (width <= 0 || height <= 0)
    {
      std::cerr << "Error: Image size is unrecognized for \"" << imageFilename << "\". Skip image.\n"
        << "width, height: " << width << ", " << height << std::endl;
      continue;
    }
    ppx = width / 2.0;
    ppy = height / 2.0;

    Exif_IO_EasyExif exifReader;
    exifReader.open( imageAbsFilepath );

    const bool bHaveValidExifMetadata =
      exifReader.doesHaveExifInfo()
      && !exifReader.getBrand().empty()
      && !exifReader.getModel().empty();

    std::map<std::string, std::string> allExifData;
    if( bHaveValidExifMetadata )
    {
      allExifData = exifReader.getExifData();
    }
    if(!exifReader.doesHaveExifInfo())
    {
      std::cerr << "No metadata for image: " << imageAbsFilepath << std::endl;
    }
    else if(exifReader.getBrand().empty() || exifReader.getModel().empty())
    {
      std::cerr << "No Brand/Model in metadata for image: " << imageAbsFilepath << std::endl;
    }

    int metadataImageWidth = imgHeader.width;
    int metadataImageHeight = imgHeader.height;
    if(allExifData.count("image_width"))
    {
      metadataImageWidth = std::stoi(allExifData.at("image_width"));
      // If the metadata is bad, use the real image size
      if(metadataImageWidth <= 0)
        metadataImageWidth = imgHeader.width;
    }
    if(allExifData.count("image_height"))
    {
      metadataImageHeight = std::stoi(allExifData.at("image_height"));
      // If the metadata is bad, use the real image size
      if(metadataImageHeight <= 0)
        metadataImageHeight = imgHeader.height;
    }
    const bool resizedImage = metadataImageWidth != width || metadataImageHeight != height;
    if(resizedImage)
    {
      std::cout << "Resized image detected:" << std::endl;
      std::cout << " * real image size: " << imgHeader.width << "x" << imgHeader.height << std::endl;
      std::cout << " * image size from metadata is: " << metadataImageWidth << "x" << metadataImageHeight << std::endl;
    }

    // Sensor width
    double ccdw = 0.0;
    if(userSensorWidth != -1.0)
    {
      ccdw = userSensorWidth;
      allExifData.emplace("sensor_width", std::to_string(ccdw));
    }
    else if(bHaveValidExifMetadata)
    {
      const std::string sCamName = exifReader.getBrand();
      const std::string sCamModel = exifReader.getModel();

      Datasheet datasheet;
      if ( getInfo( sCamName, sCamModel, vec_database, datasheet ))
      {
        // The camera model was found in the database so we can compute it's approximated focal length
        ccdw = datasheet._sensorSize;
        allExifData.emplace("sensor_width", std::to_string(ccdw));
      }
      else
      {
        error_report_stream
          << stlplus::basename_part(imageAbsFilepath) << ": Camera \""
          << sCamName << "\" model \"" << sCamModel << "\" doesn't exist in the database" << "\n"
          << "Please consider add your camera model and sensor width in the database." << "\n";
      }
    }
    else
    {
      error_report_stream
        << "No metadata in image:\n" 
        << stlplus::basename_part(imageAbsFilepath);
    }

    // Focal
    float focalLength_mm = -1.0;
    std::string sCamName;
    std::string sCamModel;
    // Consider the case where the focal is provided manually
    if (sKmatrix.size() > 0) // Known user calibration K matrix
    {
      if (!checkIntrinsicStringValidity(sKmatrix, focalPix, ppx, ppy))
        focalPix = -1.0;
    }
    // User provided focal length value
    else if (userFocalLengthPixel != -1.0)
    {
      focalPix = userFocalLengthPixel;
    }
    // If image contains meta data
    else if(bHaveValidExifMetadata)
    {
      sCamName = exifReader.getBrand();
      sCamModel = exifReader.getModel();
      focalLength_mm = exifReader.getFocal();
      // Handle case where focal length is equal to 0
      if (focalLength_mm <= 0.0f)
      {
        error_report_stream
          << stlplus::basename_part(imageAbsFilepath) << ": Focal length is missing in metadata." << "\n";
        focalPix = -1.0;
      }
      else if(ccdw != 0.0)
      {
        // Retrieve the focal from the metadata in mm and convert to pixel.
        focalPix = std::max ( metadataImageWidth, metadataImageHeight ) * focalLength_mm / ccdw;
      }
    }
    else
    {
      error_report_stream
        << stlplus::basename_part(imageAbsFilepath) << ": No metadata. The user needs to provide focal information." << "\n";
      focalPix = -1.0;
    }

    // Build intrinsic parameter related to the view
    EINTRINSIC camera_model = e_User_camera_model;
    // If no user input choose a default camera model
    if(camera_model == PINHOLE_CAMERA_START)
    {
      // Use standard lens with radial distortion by default
      camera_model = PINHOLE_CAMERA_RADIAL3;
      if(resizedImage)
      {
        // If the image has been resized, we assume that it has been undistorted
        // and we use a camera without lens distortion.
        camera_model = PINHOLE_CAMERA;
      }
      else if(focalLength_mm > 0.0 && focalLength_mm < 15)
      {
        // If the focal lens is short, the fisheye model should fit better.
        camera_model = PINHOLE_CAMERA_FISHEYE;
      }
    }

    if (focalPix <= 0 || ppx <= 0 || ppy <= 0)
    {
      std::cerr << "Warning: No intrinsic information for \"" << imageFilename << "\".\n"
        << "focal: " << focalPix << "\n"
        << "ppx,ppy: " << ppx << ", " << ppy << "\n"
        << "width,height: " << width << ", " << height
        << std::endl;
    }
    // Create the desired camera type
    std::shared_ptr<IntrinsicBase> intrinsic = createPinholeIntrinsic(camera_model, width, height, focalPix, ppx, ppy);

    // Initialize distortion parameters
    switch(camera_model)
    {
      case PINHOLE_CAMERA_FISHEYE:
      {
        if(sCamName == "GoPro")
          intrinsic->updateFromParams({focalPix, ppx, ppy, 0.0524, 0.0094, -0.0037, -0.0004});
        break;
      }
      case PINHOLE_CAMERA_FISHEYE1:
      {
        if(sCamName == "GoPro")
          intrinsic->updateFromParams({focalPix, ppx, ppy, 1.04});
        break;
      }
      default:
        break;
    }

    // Add serial number
    if( bHaveValidExifMetadata )
    {
      // Create serial number
      std::string serialNumber = "";
      serialNumber += exifReader.getSerialNumber();
      serialNumber += exifReader.getLensSerialNumber();
      intrinsic->setSerialNumber(serialNumber);
    }
    else
    {
      std::cerr
        << "Error: No instrinsics for \"" << imageFilename << "\".\n"
        << "   - focal: " << focalPix << "\n"
        << "   - ppx,ppy: " << ppx << ", " << ppy << "\n"
        << "   - width,height: " << width << ", " << height << "\n"
        << "   - camera: \"" << sCamName << "\"\n"
        << "   - model \"" << sCamModel << "\""
        << std::endl;
      if(b_Group_camera_model == 2)
      {
        // When we have no metadata at all, we create one intrinsic group per folder.
        // The use case is images extracted from a video without metadata and assumes fixed intrinsics in the video.
        intrinsic->setSerialNumber(imageFolder);
      }
    }

    IndexT id_view = views.size();
    if( b_use_UID )
    {
      const std::size_t uid = computeUID(exifReader, imageFilename);
      id_view = (IndexT)uid;
    }

    // Build the view corresponding to the image
    std::shared_ptr<View> currentView;
    if(!b_storeMetadata)
    {
      currentView.reset(new View(*iter_image, id_view, views.size(), views.size(), width, height));
    }
    else
    {
      currentView.reset(new View_Metadata(*iter_image, id_view, views.size(), views.size(), width, height, allExifData));
    }

    // Add intrinsic related to the image (if any)
    if (intrinsic == NULL)
    {
      //Since the view have invalid intrinsic data
      // (export the view, with an invalid intrinsic field value)
      currentView->id_intrinsic = UndefinedIndexT;
    }
    else
    {
      // Add the intrinsic to the sfm_container
      intrinsics[currentView->id_intrinsic] = intrinsic;
    }

    // Add the view to the sfm_container
    views[currentView->id_view] = currentView;
  }

  // Display saved warning & error messages if any.
  if (!error_report_stream.str().empty())
  {
    std::cerr
      << "\nWarning & Error messages:" << std::endl
      << error_report_stream.str() << std::endl;
  }

  // Group camera that share common properties if desired (leads to more faster & stable BA).
  if (b_Group_camera_model)
  {
    GroupSharedIntrinsics(sfm_data);
  }

  // Store SfM_Data views & intrinsic data
  if (!Save(
    sfm_data,
    stlplus::create_filespec( sOutputDir, "sfm_data.json" ).c_str(),
    ESfM_Data(VIEWS|INTRINSICS)))
  {
    return EXIT_FAILURE;
  }

  std::cout << std::endl
    << "SfMInit_ImageListing report:\n"
    << "listed #File(s): " << vec_images.size() << "\n"
    << "usable #File(s) listed in sfm_data: " << sfm_data.GetViews().size() << std::endl;

  std::size_t viewsWithoutMetatada = 0;
  for(const auto& viewValue: sfm_data.GetViews())
  {
    if(viewValue.second->id_intrinsic == UndefinedIndexT)
      ++viewsWithoutMetatada;
  }
  if(viewsWithoutMetatada == sfm_data.GetViews().size())
  {
    std::cerr << "ERROR: No metadata in images." << std::endl;
    return EXIT_FAILURE;
  }
  else if(viewsWithoutMetatada)
  {
    std::cerr << "WARNING: " << viewsWithoutMetatada << " views without metadata. It may fail the reconstruction." << std::endl;
  }
  return EXIT_SUCCESS;
}
