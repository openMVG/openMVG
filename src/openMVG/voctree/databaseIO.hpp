#pragma once

#include "database.hpp"
#include "vocabulary_tree.hpp"

#include <openMVG/features/descriptor.hpp>
#include <openMVG/sfm/sfm_data_io.hpp>

#include <string>

namespace openMVG {
namespace voctree {

/**
 * @brief Given a vocabulary tree and a set of features it builds a database
 *
 * @param[in] fileFullPath A file containing the path the features to load, it could be a .txt or an OpenMVG .json
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
                             std::map<size_t, Document> &documents,
                             const int Nmax = 0);

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
 * @param[in] distanceMethod The distance method used for create the pair list
 * @param[in] Nmax The maximum number of features loaded in each desc file. For Nmax = 0 (default), all the descriptors are loaded. 
 * @see queryDatabase()
 */
template<class DescriptorT, class VocDescriptorT>
void queryDatabase(const std::string &fileFullPath,
                   const openMVG::voctree::VocabularyTree<VocDescriptorT> &tree,
                   const openMVG::voctree::Database &db,
                   size_t numResults,
                   std::map<size_t, DocMatches> &allMatches,
                   const std::string &distanceMethod,
                   const int Nmax = 0);

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
 * @param[in] distanceMethod The distance method used for create the pair list
 * @param[in] Nmax The maximum number of features loaded in each desc file. For Nmax = 0 (default), all the descriptors are loaded.
 */
template<class DescriptorT, class VocDescriptorT>
void queryDatabase(const std::string &fileFullPath,
                   const VocabularyTree<VocDescriptorT> &tree,
                   const Database &db,
                   size_t numResults,
                   std::map<size_t, DocMatches> &allMatches,
                   std::map<size_t, Document> &documents,
                   const std::string &distanceMethod,
                   const int Nmax = 0);

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
    std::map<int, int> &globalHistogram);

} //namespace voctree
} //namespace openMVG

#include "databaseIO.tcc"