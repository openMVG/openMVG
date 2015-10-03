#include "descriptor_loader.hpp"

#include <boost/filesystem.hpp>
#include <boost/algorithm/string/case_conv.hpp>
#include <boost/progress.hpp>

#include <exception>
#include <iostream>
#include <fstream>

namespace openMVG {
namespace voctree {

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
  std::cout << "Reading the descriptors..." << std::endl;
  boost::progress_display display(descriptorsFiles.size());
  size_t docId = 0;

  // Run through the path vector and read the descriptors
  //for(vector<string>::const_iterator it = descriptorsFiles.begin(); it != descriptorsFiles.end(); ++it, ++display, ++docId)
  for(const auto &currentFile : descriptorsFiles)
  {
    std::vector<DescriptorT> descriptors;
    openMVG::voctree::Matches matches;

    // Read the descriptors
    loadDescsFromBinFile(currentFile, descriptors, false);

    std::vector<openMVG::voctree::Word> imgVisualWords=tree.quantize(descriptors);

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