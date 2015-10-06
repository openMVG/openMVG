#include "descriptor_loader.hpp"

namespace openMVG {
namespace voctree {

void getInfoBinFile(const std::string &path, int dim, size_t &numDescriptors, int &bytesPerElement)
{
  std::fstream fs;

  // the file is supposed to have the number of descriptors as first element and then
  // the set of descriptors of dimension dim either as chars or floats

  // Open file and get the number of descriptors
  fs.open(path, std::ios::in | std::ios::binary);

  if(!fs.is_open())
  {
    std::cerr << "Error while opening " << path << std::endl;
    std::cerr << "Error while opening " + path << std::endl;
  }

  // go to the end of the file
  fs.seekg(0, fs.end);

  // get the length in byte of the file
  //@fixeme we are ignoring the first element of the file which is the number of
  // feature. However given the large amount of data of the feature this is mitigate
  // by the integer division in bytepd later
  int length = fs.tellg();

  // go back to the beginning of the file
  fs.seekg(0, fs.beg);

  // get the number of descriptors
  fs.read((char*) &numDescriptors, sizeof (size_t));

  if(numDescriptors > 0)
  {
    // get the number of bytes per descriptor element
    bytesPerElement = (length / numDescriptors) / dim;
  }
  else
  {
    bytesPerElement = 0;
  }
}

}
}
