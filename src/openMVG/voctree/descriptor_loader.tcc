#include <openMVG/sfm/sfm_data_io.hpp>

#include <boost/filesystem.hpp>
#include <boost/algorithm/string/case_conv.hpp>
#include <boost/progress.hpp>

#include <iostream>
#include <fstream>

namespace openMVG {
namespace voctree {

template<class DescriptorT>
size_t readDescFromFiles(const std::string &fileFullPath, std::vector<DescriptorT>& descriptors, std::vector<size_t> &numFeatures)
{
  namespace bfs = boost::filesystem;
  std::map<IndexT, std::string> descriptorsFiles;
  getListOfDescriptorFiles(fileFullPath, descriptorsFiles);
  std::size_t numDescriptors = 0;

  // Allocate the memory by reading in a first time the files to get the number
  // of descriptors
  int bytesPerElement = 0;

  // Display infos and progress bar
  std::cout << "Pre-computing the memory needed..." << std::endl;
  boost::progress_display display(descriptorsFiles.size());

  // Read all files and get the number of descriptors to load
  for(const auto &currentFile : descriptorsFiles)
  {
    // if it is the first one read the number of descriptors and the type of data (we are assuming the the feat are all the same...)
    // bytesPerElement could be 0 even after the first element (eg it has 0 descriptors...), so do it until we get the correct info
    if(bytesPerElement == 0)
    {
      getInfoBinFile(currentFile.second, DescriptorT::static_size, numDescriptors, bytesPerElement);
    }
    else
    {
      // get the file size in byte and estimate the number of features without opening the file
      numDescriptors += (bfs::file_size(currentFile.second) / bytesPerElement) / DescriptorT::static_size;
    }
    ++display;
  }
  std::cout << "Found " << numDescriptors << " descriptors overall, allocating memory..." << std::endl;
  if(bytesPerElement == 0)
  {
    std::cout << "WARNING: Empty descriptor file: " << fileFullPath << std::endl;
    return 0;
  }
  
  else
  {
    // Allocate the memory
    descriptors.reserve(numDescriptors);
    size_t numDescriptorsCheck = numDescriptors; // for later check
    numDescriptors = 0;

    // Read the descriptors
    std::cout << "Reading the descriptors..." << std::endl;
    display.restart(descriptorsFiles.size());

    // Run through the path vector and read the descriptors
    for(const auto &currentFile : descriptorsFiles)
    {
      // Read the descriptors and append them in the vector
      loadDescsFromBinFile(currentFile.second, descriptors, true);
      size_t result = descriptors.size();

      // Add the number of descriptors from this file
      numFeatures.push_back(result - numDescriptors);

      // Update the overall counter
      numDescriptors = result;

      ++display;
    }
    BOOST_ASSERT(numDescriptors == numDescriptorsCheck);

  }
  // Return the result
  return numDescriptors;
}

} // namespace voctree
} // namespace openMVG