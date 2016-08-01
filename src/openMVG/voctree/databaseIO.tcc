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
 * @param[in] fileFullPath A file containing the path the features to load, 
 * it could be a .txt or an OpenMVG .json or the path to a directory containing
 * the feature. In the case of a json, it is assumed that the descriptors are
 * in the same folder as the file. 
 * @param[in] tree The vocabulary tree to be used for feature quantization
 * @param[out] db The built database
 * @param[out] documents A map containing for each image the list of associated visual words
 * @param[in] Nmax The maximum number of features loaded in each desc file. For Nmax = 0 (default), all the descriptors are loaded.
 * @return the number of overall features read
 */
template<class DescriptorT, class VocDescriptorT>
std::size_t populateDatabase(const std::string &fileFullPath,
                             const VocabularyTree<VocDescriptorT> &tree,
                             Database &db,
                             const int Nmax)
{
  std::map<IndexT, std::string> descriptorsFiles;
  getListOfDescriptorFiles(fileFullPath, descriptorsFiles);
  std::size_t numDescriptors = 0;
  
  // Read the descriptors
  std::cout << "Reading the descriptors from " << descriptorsFiles.size() <<" files..." << std::endl;
  boost::progress_display display(descriptorsFiles.size());

  // Run through the path vector and read the descriptors
  for(const auto &currentFile : descriptorsFiles)
  {
    
    std::vector<DescriptorT> descriptors;

    // Read the descriptors
    loadDescsFromBinFile(currentFile.second, descriptors, false, Nmax);
    size_t result = descriptors.size();
    
    SparseHistogram newDoc =  tree.quantizeToSparse(descriptors);

    // Insert document in database
    db.insert(currentFile.first, newDoc);

    // Update the overall counter
    numDescriptors += result;
    
    ++display;
  }

  // Return the result
  return numDescriptors;
}

template<class DescriptorT, class VocDescriptorT>
std::size_t populateDatabase(const std::string &fileFullPath,
                             const VocabularyTree<VocDescriptorT> &tree,
                             Database &db,
                             std::map<size_t, std::vector<DescriptorT>> &allDescriptors,
                             const int Nmax)
{
  std::map<IndexT, std::string> descriptorsFiles;
  getListOfDescriptorFiles(fileFullPath, descriptorsFiles);
  std::size_t numDescriptors = 0;

  // Read the descriptors
  std::cout << "Reading the descriptors from " << descriptorsFiles.size() <<" files..." << std::endl;
  boost::progress_display display(descriptorsFiles.size());

  // Run through the path vector and read the descriptors
  for(const auto &currentFile : descriptorsFiles)
  {
    std::vector<DescriptorT> descriptors;

    // Read the descriptors
    loadDescsFromBinFile(currentFile.second, descriptors, false, Nmax);
    size_t result = descriptors.size();
    
    allDescriptors[currentFile.first] = descriptors;
    
    SparseHistogram newDoc = tree.quantizeToSparse(descriptors);
    
    // Insert document in database
    db.insert(currentFile.first, newDoc);

    // Update the overall counter
    numDescriptors += result;

    ++display;
  }

  // Return the result
  return numDescriptors;
}

/**
 * @brief Given an non empty database, it queries the database with a set of images
 * and their associated features and returns, for each image, the first \p numResults best
 * matching documents in the database
 * 
 * @param[in] fileFullPath A file containing the path the features to load, it could be a .txt or an OpenMVG .json
 * @param[in] tree The vocabulary tree to be used for feature quantization
 * @param[in] db The built database
 * @param[in] numResults The number of results to retrieve for each image
 * @param[out] allMatches The matches for all the images
 * @param[in] distanceMethod the method used to compute distance between histograms.
 * @param[in] Nmax The maximum number of features loaded in each desc file. For Nmax = 0 (default), all the descriptors are loaded. 
 * @see queryDatabase()
 */
template<class DescriptorT, class VocDescriptorT>
void queryDatabase(const std::string &fileFullPath,
                   const VocabularyTree<VocDescriptorT> &tree,
                   const Database &db,
                   size_t numResults,
                   std::map<size_t, DocMatches> &allDocMatches,
                   const std::string &distanceMethod,
                   const int Nmax)
{
  std::map<size_t, Document> documents;
  queryDatabase<DescriptorT>(fileFullPath, tree, db, numResults, allDocMatches, documents, distanceMethod, Nmax);
}

/**
 * @brief Given an non empty database, it queries the database with a set of images
 * and their associated features and returns, for each image, the first \p numResults best
 * matching documents in the database
 * 
 * @param[in] fileFullPath A file containing the path the features to load, it could be a .txt or an OpenMVG .json
 * @param[in] tree The vocabulary tree to be used for feature quantization
 * @param[in] db The built database
 * @param[in] numResults The number of results to retrieve for each image
 * @param[out] allMatches The matches for all the images
 * @param[out] documents For each document, it contains the list of associated visual words 
 * @param[in] distanceMethod The method used to compute distance between histograms.
 * @param[in] Nmax The maximum number of features loaded in each desc file. For Nmax = 0 (default), all the descriptors are loaded.
 */
template<class DescriptorT, class VocDescriptorT>
void queryDatabase(const std::string &fileFullPath,
                   const VocabularyTree<VocDescriptorT> &tree,
                   const Database &db,
                   size_t numResults,
                   std::map<size_t, DocMatches> &allDocMatches,
                   std::map<size_t, SparseHistogram> &documents,
                   const std::string &distanceMethod,
                   const int Nmax)
{
  std::map<IndexT, std::string> descriptorsFiles;
  getListOfDescriptorFiles(fileFullPath, descriptorsFiles);
  
  // Read the descriptors
  std::cout << "queryDatabase: Reading the descriptors from " << descriptorsFiles.size() <<" files..." << std::endl;
  boost::progress_display display(descriptorsFiles.size());

  #ifdef OPENMVG_USE_OPENMP
    #pragma omp parallel for
  #endif
  // Run through the path vector and read the descriptors
  for(std::size_t i = 0; i < descriptorsFiles.size(); ++i)
  {
    std::map<IndexT, std::string>::const_iterator currentFileIt = descriptorsFiles.cbegin();
    std::advance(currentFileIt, i);
    std::vector<DescriptorT> descriptors;

    // Read the descriptors
    loadDescsFromBinFile(currentFileIt->second, descriptors, false, Nmax);

    // quantize the descriptors
    SparseHistogram query = tree.quantizeToSparse(descriptors);

    openMVG::voctree::DocMatches docMatches;
    // query the database
    db.find(query, numResults, docMatches, distanceMethod);

    #ifdef OPENMVG_USE_OPENMP
      #pragma omp critical
    #endif
    {
      // add the vector to the documents
      documents[currentFileIt->first] = query;

      // add the matches to the result vector
      allDocMatches[currentFileIt->first] = docMatches;

      ++display;
    }
  }
}

/**
 * @brief Returns some statistics (histogram) 
 * 
 * @param[in] fileFullPath A file containing the path the features to load, it could be a .txt or an OpenMVG .json
 * @param[in] tree The vocabulary tree to be used for feature quantization
 * @param[in] db The built database
 * @param[in] distanceMethod The distance method used for create the pair list
 * @param[in/out] globalHistogram The histogram of the "population" of voctree leaves. 
 * @see queryDatabase()
 */
template<class DescriptorT, class VocDescriptorT>
void voctreeStatistics(
    const std::string &fileFullPath,
    const VocabularyTree<VocDescriptorT> &tree,
    const Database &db,
    const std::string &distanceMethod,
    std::map<int, int> &globalHistogram)
{
  std::map<IndexT, std::string> descriptorsFiles;
  getListOfDescriptorFiles(fileFullPath, descriptorsFiles);
  
  // Read the descriptors
  std::cout << "Reading the descriptors from " << descriptorsFiles.size() << " files." << std::endl;

  // Run through the path vector and read the descriptors
  for(const auto &currentFile : descriptorsFiles)
  {
    std::vector<DescriptorT> descriptors;

    // Read the descriptors
    loadDescsFromBinFile(currentFile.second, descriptors, false);

    // query the database
    SparseHistogram query = tree.quantizeToSparse(descriptors);
    std::map<int,int> localHisto;
    
    for(auto q: query)
    {
      int nb = (int)q.second.size();
      if(globalHistogram.find(nb) == globalHistogram.end())
        globalHistogram[nb] = 1;
      else
        globalHistogram[nb] += 1;
      
      if(localHisto.find(nb) == localHisto.end())
        localHisto[nb] = 1;
      else
        localHisto[nb] += 1;   
    }
    
    std::cout << "Histogram of " << currentFile.first << std::endl;
    
    for(auto itHisto = localHisto.begin(); itHisto != localHisto.end(); itHisto++)
    {
      std::cout << itHisto->first << ": " << itHisto->second  << ", ";
    }
    
    localHisto.clear();
    
    std::cout << std::endl;
    
  }
}


} //namespace voctree
} //namespace openMVG