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
 * @brief Get the number of descriptors contained inside a .desc file and the number of bytes
 * used to store each descriptor elements
 *
 * @param[in] path The .desc filename
 * @param[in] dim The number of elements per descriptor
 * @param[out] numDescriptors The number of descriptors stored in the file
 * @param[out] bytesPerElement The number of bytes used to store each element of the descriptor
 */
void getInfoBinFile(const std::string &path, int dim, size_t &numDescriptors, int &bytesPerElement);

/**
 * @brief Extract a list of decriptor files to read from a file containing the path to the descriptor files.
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
 * @param[in] listFile The input filename containing the list of files to read
 * @param[out] descriptorsFiles A list of descriptor files 
 */
void getListOfDescriptorFiles(const std::string &listFile, std::vector<std::string> &descriptorsFiles);

/**
 * @brief Read a set of descriptors from a file containing the path to the descriptor files.
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
 * @param[in] fileFullPath the input filename containing the list of files to read
 * @param[in,out] descriptors the vector to which append all the read descriptors
 * @param[in,out] numFeatures a vector collecting for each file read the number of features read
 * @return the total number of features read
 *
 */
template<class DescriptorT>
size_t readDescFromFiles(const std::string &fileFullPath, std::vector<DescriptorT>& descriptors, std::vector<size_t> &numFeatures)
{
  namespace bfs = boost::filesystem;
  std::vector<std::string> descriptorsFiles;
  getListOfDescriptorFiles(fileFullPath, descriptorsFiles);
  std::size_t numDescriptors = 0;

  // Allocate the memory by reading in a first time the files to get the number
  // of descriptors
  int bytesPerElement = 0;

  // Display infos and progress bar
  std::cout << "Pre-computing the memory needed..." << std::endl;
  boost::progress_display display(descriptorsFiles.size());

  // Read all files and get the number of descriptors to load
  for(std::vector<std::string>::const_iterator it = descriptorsFiles.begin(); it != descriptorsFiles.end(); ++it)
  for(const auto &currentFile : descriptorsFiles)
  {
    // if it is the first one read the number of descriptors and the type of data (we are assuming the the feat are all the same...)
    // bytesPerElement could be 0 even after the first element (eg it has 0 descriptors...), so do it until we get the correct info
    if(bytesPerElement == 0)
    {
      getInfoBinFile(currentFile, DescriptorT::static_size, numDescriptors, bytesPerElement);
    }
    else
    {
      // get the file size in byte and estimate the number of features without opening the file
      numDescriptors += (bfs::file_size(*it) / bytesPerElement) / DescriptorT::static_size;
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
  std::vector<std::string> descriptorsFiles;
  getListOfDescriptorFiles(fileFullPath, descriptorsFiles);
  std::size_t numDescriptors = 0;
  
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
