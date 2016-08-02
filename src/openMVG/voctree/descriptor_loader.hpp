#pragma once

#include "database.hpp"
#include "vocabulary_tree.hpp"
#include <openMVG/types.hpp>
#include <openMVG/features/descriptor.hpp>

#include <string>


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
 * Three formats are supported:
 * 
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
 * 3. a path to a directory containing the features and the descriptors. The function
 * looks for all .desc files there inside.
 *  
 * @param[in] listFile The input filename containing the list of files to read
 * @param[out] descriptorsFiles A list of descriptor files 
 */
void getListOfDescriptorFiles(const std::string &fileFullPath, std::map<IndexT, std::string> &descriptorsFiles);

/**
 * @brief Read a set of descriptors from a file containing the path to the descriptor files.
 * Three formats are supported:
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
 * filename from that. It assumes that the descriptors are in the same folder as
 * the .json file.
 * 
 * 3. a path to a directory containing the features and the descriptors. The function
 * looks for all .desc files there inside.
 * 
 * \p FileDescriptorT is the type of descriptors that are stored in the file. Usually 
 * the two types should be the same, but it could be the case in which the descriptor 
 * stored in the file  * has different type representation: for example the file 
 * could contain SIFT descriptors stored as uchar (the default type) and we want 
 * to cast these into SIFT descriptors stored in memory as floats.
 *
 * @param[in] fileFullPath the input filename containing the list of files to read
 * @param[in,out] descriptors the vector to which append all the read descriptors
 * @param[in,out] numFeatures a vector collecting for each file read the number of features read
 * @return the total number of features read
 *
 */
template<class DescriptorT, class FileDescriptorT>
size_t readDescFromFiles(const std::string &fileFullPath, std::vector<DescriptorT>& descriptors, std::vector<size_t> &numFeatures);

} // namespace voctree
} // namespace openMVG

#include "descriptor_loader.tcc"
