#include "descriptor_loader.hpp"

#include <openMVG/sfm/sfm_data_io.hpp>

#include <boost/filesystem.hpp>
#include <boost/algorithm/string/case_conv.hpp>
#include <boost/progress.hpp>

#include <exception>
#include <iostream>
#include <fstream>

namespace openMVG {
namespace voctree {

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
  std::cout << "Reading the descriptors from " << descriptorsFiles.size() <<" files..." << std::endl;
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
    
    ++docId;
    ++display;
  }

  // Return the result
  return numDescriptors;
}

template<class DescriptorT>
void queryDatabase(const std::string &fileFullPath,
                   const openMVG::voctree::VocabularyTree<DescriptorT> &tree,
                   const openMVG::voctree::Database &db,
                   size_t numResults,
                   std::vector<openMVG::voctree::Matches> &allMatches)
{
  std::map<size_t, openMVG::voctree::Document> documents;
  queryDatabase(fileFullPath, tree, db, numResults, allMatches, documents);
}


template<class DescriptorT>
void queryDatabase(const std::string &fileFullPath,
                   const openMVG::voctree::VocabularyTree<DescriptorT> &tree,
                   const openMVG::voctree::Database &db,
                   size_t numResults,
                   std::vector<openMVG::voctree::Matches> &allMatches,
                   std::map<size_t, openMVG::voctree::Document> &documents)
{  
  std::vector<std::string> descriptorsFiles;
  getListOfDescriptorFiles(fileFullPath, descriptorsFiles);
  
  // Read the descriptors
  std::cout << "Reading the descriptors from " << descriptorsFiles.size() <<" files..." << std::endl;
  boost::progress_display display(descriptorsFiles.size());
  size_t docId = 0;

  // Run through the path vector and read the descriptors
  for(const auto &currentFile : descriptorsFiles)
  {
    std::vector<DescriptorT> descriptors;
    openMVG::voctree::Matches matches;

    // Read the descriptors
    loadDescsFromBinFile(currentFile, descriptors, false);
    
    // quantize the descriptors
    std::vector<openMVG::voctree::Word> imgVisualWords = tree.quantize(descriptors);

    // add the vector to the documents
    documents[docId] = imgVisualWords;

    // query the database
    db.find(imgVisualWords, numResults, matches);

    // add the matches to the result vector
    allMatches.emplace_back(matches);
    
    ++display;
    ++docId;
  }

}


} //namespace voctree
} //namespace openMVG