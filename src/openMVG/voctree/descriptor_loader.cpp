#include "descriptor_loader.hpp"
#include <boost/algorithm/string/predicate.hpp>

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

void getListOfDescriptorFiles(const std::string &fileFullPath, std::map<IndexT, std::string> &descriptorsFiles)
{
  namespace bfs = boost::filesystem;
  std::ifstream fs;
  bfs::path pathToFiles;
  
  descriptorsFiles.clear();

  bfs::path bp(fileFullPath);

  // If the input is a directory, list all .desc files recursively.
  if(bfs::is_directory(bp))
  {
    std::size_t viewId = 0;
    for(bfs::recursive_directory_iterator it(bp), end; it != end; ++it)
    {
        if(!bfs::is_directory(*it) && boost::algorithm::ends_with(it->path().string(), ".desc"))
        {
          descriptorsFiles[viewId++] = it->path().string();
        }
    }
    return;
  }
  
  if(!bp.has_extension())
  {
    std::cerr << "File without extension not recognized! " << fileFullPath << std::endl;
    std::cerr << "The file  " + fileFullPath + " is neither a JSON nor a txt file" << std::endl;
    throw std::invalid_argument("Unrecognized extension for " + fileFullPath);
  }

  // get the extension of the file and put it lowercase
  std::string ext = bp.extension().string();
  boost::to_lower(ext);

  // two cases, either the input file is a text file with the relative paths or
  // it is a JSON file from OpenMVG
  // in the two cases we fill a vector with paths to the descriptors files

  if(ext == ".txt")
  {
    // processing a file .txt containing the relative paths

    // Extract the folder path from the list file path
    pathToFiles = bfs::path(fileFullPath).parent_path();

    // Open file
    fs.open(fileFullPath, std::ios::in);
    if(!fs.is_open())
    {
      std::cerr << "Error while opening " << fileFullPath << std::endl;
      throw std::invalid_argument("Error while opening " + fileFullPath);
    }

    // read the file line by line and store in the vector the descriptors paths
    std::string line;
    IndexT viewId = 0;
    while(getline(fs, line))
    {
      // extract the filename without extension and create on with .desc as extension
      const std::string filename = bfs::path(line).stem().string() + ".desc";
      const std::string filepath = (pathToFiles / filename).string();

      // add the filepath in the vector
      descriptorsFiles[viewId++] = filepath;
    }
  }
  else
  {
    // processing a JSON file containing sfm_data

    // open the sfm_data file
    openMVG::sfm::SfM_Data sfmdata;
    openMVG::sfm::Load(sfmdata, fileFullPath, openMVG::sfm::ESfM_Data::VIEWS);

    // get the number of files to load
    size_t numberOfFiles = sfmdata.GetViews().size();

    if(numberOfFiles == 0)
    {
      std::cout << "It seems like there are no views in " << fileFullPath << std::endl;
      return;
    }

    // get the base path for the files
    pathToFiles = bfs::path(fileFullPath).parent_path();

    // explore the sfm_data container to get the files path
    for(const auto &view : sfmdata.GetViews())
    {
      // generate the equivalent .desc file path
      const std::string filepath = bfs::path(pathToFiles / (std::to_string(view.first) + ".desc")).string();

      // add the filepath in the vector
      descriptorsFiles[view.first] = filepath;
    }
  }
}

}
}
