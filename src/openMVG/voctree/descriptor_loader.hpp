#pragma once

#include <openMVG/features/descriptor.hpp>
#include <openMVG/sfm/sfm_data_io.hpp>
#include "database.hpp"
#include "vocabulary_tree.hpp"

#include <boost/filesystem.hpp>
#include <boost/algorithm/string/case_conv.hpp>
#include <boost/progress.hpp>

#include <iostream>
#include <fstream>
#include <string>

//using namespace std;

namespace openMVG {
namespace voctree {

/**
 * Get the number of descriptors contained inside a .desc file and the number of bytes
 * used to store each descriptor elements
 *
 * @param[in] path The .desc filename
 * @param[in] dim The number of elements per descriptor
 * @param[out] numDescriptors The number of descriptors stored in the file
 * @param[out] bytesPerElement The number of bytes used to store each element of the descriptor
 */
void getInfoBinFile(const std::string &path, int dim, size_t &numDescriptors, int &bytesPerElement);

/**
 * Read a set of descriptors from a file containing the path to the descriptor files.
 * Two format are supported:
 * 1. a simple text file containing a list of filenames of images or descriptors, one for
 * each line, like:
 * imagename1.jpg
 * imagename2.jpg
 * imagename3.jpg
 * or
 * imagename1.desc
 * imagename2.desc
 * imagename3.desc
 * In any case the filename for the descriptors will be inferred by removing the extension
 * and keeping the name. Normally this is the format used by Bundler
 *
 * 2. a json file containing the sfm_data using the OpenMVG data container. The function will
 * parse the view section to retrieve the image name and it will infer the descriptor
 * filename from that
 *
 * @param[in] fileFullPath the input filename containing the list of file to read
 * @param[in,out] descriptors the vector to which append all the read descriptors
 * @param[in,out] numFeatures a vector collecting for each file read the number of features read
 * @return the total number of features read
 *
 */
template<class DescriptorT>
size_t readDescFromFiles(const std::string &fileFullPath, std::vector<DescriptorT>& descriptors, std::vector<size_t> &numFeatures)
{
  namespace boostfs = boost::filesystem;
  std::ifstream fs;
  boostfs::path pathToFiles;

  size_t numDescriptors = 0;
  boostfs::path bp(fileFullPath);

  if(!bp.has_extension())
  {
    std::cerr << "File without extension not recognized! " << fileFullPath << std::endl;
    std::cerr << "The file  " + fileFullPath + " is neither a JSON nor a txt file" << std::endl;
    return numDescriptors;
  }

  // get the extension of the file and put it lowercase
  std::string ext = bp.extension().string();
  boost::to_lower(ext);

  // two cases, either the input file is a text file with the relative paths or
  // it is a JSON file from OpenMVG
  // in the two cases we fill a vector with paths to the descriptors files

  std::vector<std::string> descriptorsFiles;

  // if it is a JSON file
  if(ext == ".json")
  {
    // processing a JSON file containing sfm_data

    // open the sfm_data file
    openMVG::sfm::SfM_Data sfmdata;
    openMVG::sfm::Load(sfmdata, fileFullPath, openMVG::sfm::ESfM_Data::VIEWS);

    // get the number of files to load
    size_t numberOfFiles = sfmdata.GetViews().size();

    // Reserve memory for the file path vector
    descriptorsFiles.reserve(numberOfFiles);

    if(numberOfFiles == 0)
    {
      std::cout << "It seems like there are no views in " << fileFullPath << std::endl;
      return 0;
    }

    // get the base path for the files
    pathToFiles = boostfs::path(fileFullPath).parent_path();

    // explore the sfm_data container to get the files path
    for(openMVG::sfm::Views::const_iterator it = sfmdata.GetViews().begin(); it != sfmdata.GetViews().end(); ++it)
    {
      // get just the image name, remove the extension
      std::string filepath = boostfs::path(it->second->s_Img_path).stem().string();

      // generate the equivalent .desc file path
      filepath = boostfs::path(pathToFiles / boostfs::path(filepath + ".desc")).string();

      // add the filepath in the vector
      descriptorsFiles.push_back(filepath);
    }
  }
  else if(ext == ".txt")
  {
    // processing a file .txt containing the relative paths

    // Extract the folder path from the list file path
    pathToFiles = boostfs::path(fileFullPath).parent_path();

    // Open file
    fs.open(fileFullPath, std::ios::in);
    if(!fs.is_open())
    {
      std::cerr << "Error while opening " << fileFullPath << std::endl;
      return numDescriptors;
    }

    // read the file line by line and store in the vector the descriptors paths
    std::string line;
    while(getline(fs, line))
    {
      // we have to do that because OMVG does not really output a clean list.txt, it also
      // contains other stuff, so we look at the first '.' to get the extension (not robust at all)
      std::string filepath = line.substr(0, line.find_first_of("."));
      filepath = boostfs::path(pathToFiles / boostfs::path(filepath + ".desc")).string();

      // add the filepath in the vector
      descriptorsFiles.push_back(filepath);
    }
  }
  else
  {
    std::cerr << "File not recognized! " << fileFullPath << std::endl;
    std::cerr << "The file  " + fileFullPath + " is neither a JSON nor a txt file" << std::endl;
    return numDescriptors;
  }

  // Allocate the memory by reading in a first time the files to get the number
  // of descriptors
  int bytesPerElement = 0;

  // Display infos and progress bar
  std::cout << "Pre-computing the memory needed..." << std::endl;
  boost::progress_display display(descriptorsFiles.size());

  // Read all files and get the number of descriptors to load
  for(std::vector<std::string>::const_iterator it = descriptorsFiles.begin(); it != descriptorsFiles.end(); ++it)
  {
    // if it is the first one read the number of descriptors and the type of data (we are assuming the the feat are all the same...)
    // bytesPerElement could be 0 even after the first element (eg it has 0 descriptors...), so do it until we get the correct info
    if(bytesPerElement == 0)
    {
      getInfoBinFile(*it, DescriptorT::static_size, numDescriptors, bytesPerElement);
    }
    else
    {
      // get the file size in byte and estimate the number of features without opening the file
      numDescriptors += (boostfs::file_size(*it) / bytesPerElement) / DescriptorT::static_size;
    }
  }
  BOOST_ASSERT(bytesPerElement > 0);
  std::cout << "Found " << numDescriptors << " descriptors overall, allocating memory..." << std::endl;

  // Allocate the memory
  descriptors.reserve(numDescriptors);
  size_t numDescriptorsCheck = numDescriptors; // for later check
  numDescriptors = 0;

  // Read the descriptors
  std::cout << "Reading the descriptors..." << std::endl;
  display.restart(descriptorsFiles.size());

  // Run through the path vector and read the descriptors
  for(std::vector<std::string>::const_iterator it = descriptorsFiles.begin(); it != descriptorsFiles.end(); ++it, ++display)
  {
    // Read the descriptors and append them in the vector
    loadDescsFromBinFile(*it, descriptors, true);
    size_t result = descriptors.size();

    // Add the number of descriptors from this file
    numFeatures.push_back(result - numDescriptors);

    // Update the overall counter
    numDescriptors = result;
  }
  BOOST_ASSERT(numDescriptors == numDescriptorsCheck);

  // Return the result
  return numDescriptors;
}

/**
 * @brief Given a vocabulary tree and a set of features it builds a database
 *
 * @param[in] fileFullPath A file containing the path the features to load, it could be a .txt or an OpenMVG .json
 * @param[in] tree The vocabulary tree to be used for feature quantization
 * @param[out] db The built database
 * @param[out] documents A map containing for each image the list of associated visual words
 * @param[in,out] numFeatures a vector collecting for each file read the number of features read
 * @return the number of overall features read
 */
template<class DescriptorT>
std::size_t populateDatabase(const std::string &fileFullPath,
                             const VocabularyTree<DescriptorT> &tree,
                             Database &db,
                             std::map<size_t, Document> &documents,
                             std::vector<size_t> &numFeatures)
{
  namespace bfs = boost::filesystem;
  std::ifstream fs;
  bfs::path pathToFiles;

  size_t numDescriptors = 0;
  bfs::path bp(fileFullPath);

  if(!bp.has_extension())
  {
    std::cerr << "File without extension not recognized! " << fileFullPath << std::endl;
    std::cerr << "The file  " + fileFullPath + " is neither a JSON nor a txt file" << std::endl;
    return numDescriptors;
  }

  // get the extension of the file and put it lowercase
  std::string ext = bp.extension().string();
  boost::to_lower(ext);

  // two cases, either the input file is a text file with the relative paths or
  // it is a JSON file from OpenMVG
  // in the two cases we fill a vector with paths to the descriptors files

  std::vector<std::string> descriptorsFiles;

  // if it is a JSON file
  if(ext == ".json")
  {
    // processing a JSON file containing sfm_data

    // open the sfm_data file
    openMVG::sfm::SfM_Data sfmdata;
    openMVG::sfm::Load(sfmdata, fileFullPath, openMVG::sfm::ESfM_Data::VIEWS);

    // get the number of files to load
    std::size_t numberOfFiles = sfmdata.GetViews().size();

    // Reserve memory for the file path vector
    descriptorsFiles.reserve(numberOfFiles);

    if(numberOfFiles == 0)
    {
      std::cout << "It seems like there are no views in " << fileFullPath << std::endl;
      return 0;
    }

    // get the base path for the files
    pathToFiles = bfs::path(fileFullPath).parent_path();

    // explore the sfm_data container to get the files path
    for(const auto &view : sfmdata.GetViews())
    {
      // get just the image name, remove the extension
      std::string filepath = bfs::path(view.second->s_Img_path).stem().string();

      // generate the equivalent .desc file path
      filepath = bfs::path(pathToFiles / (filepath + ".desc")).string();

      // add the filepath in the vector
      descriptorsFiles.push_back(filepath);
    }
  }
  else if(ext == ".txt")
  {
    // processing a file .txt containing the relative paths

    // Extract the folder path from the list file path
    pathToFiles = bfs::path(fileFullPath).parent_path();

    // Open file
    fs.open(fileFullPath, std::ios::in);
    if(!fs.is_open())
    {
      std::cerr << "Could not found any file to load in " << fileFullPath <<"..." << std::endl;
      throw std::invalid_argument("Error while opening " + fileFullPath);
    }

    // read the file line by line and store in the vector the descriptors paths
    std::string line;
    while(getline(fs, line))
    {
      // we have to do that because OMVG does not really output a clean list.txt, it also
      // contains other stuff, so we look at the first '.' to get the extension (not robust at all)
      std::string filepath = line.substr(0, line.find_first_of("."));
      filepath = bfs::path(pathToFiles / bfs::path(filepath + ".desc")).string();

      // add the filepath in the vector
      descriptorsFiles.push_back(filepath);
    }
  }
  else
  {
    std::cerr << "File not recognized!" << fileFullPath << std::endl;
    std::cerr << "The file  " + fileFullPath + " is neither a JSON nor a txt file" << std::endl;
    throw std::invalid_argument("Unrecognized file format " + fileFullPath);
  }

  // Read the descriptors
  std::cout << "Reading the descriptors..." << std::endl;
  boost::progress_display display(descriptorsFiles.size());
  size_t docId = 0;
  numFeatures.resize(descriptorsFiles.size());

  // Run through the path vector and read the descriptors
  for(const auto &currentFile : descriptorsFiles)
  {
    std::vector<DescriptorT> descriptors;

    // Read the descriptors
    loadDescsFromBinFile(currentFile, descriptors, false);
    size_t result = descriptors.size();
    
    std::vector<openMVG::voctree::Word> imgVisualWords = tree.quantize(descriptors);

    // Add the vector to the documents
    documents[docId] = imgVisualWords;

    // Insert document in database
    db.insert(imgVisualWords);

    // Update the overall counter
    numDescriptors += result;

    // Save the number of features of this image
    numFeatures[docId] = result;
  }

  // Return the result
  return numDescriptors;
}

}
}
