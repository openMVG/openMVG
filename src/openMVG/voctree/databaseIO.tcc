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
  namespace bfs = boost::filesystem;
  std::ifstream fs;
  bfs::path pathToFiles;
  std::string line;

  bfs::path bp(fileFullPath);

  if(!bp.has_extension())
  {
    std::cerr << "File without extension not recognized! " << fileFullPath << std::endl;
    std::cerr << "The file  " + fileFullPath + " is neither a JSON nor a txt file" << std::endl;
  }

  // get the extension of the file and put it lowercase
  std::string ext = bp.extension().string();
  boost::to_lower(ext);

  // two cases, either the input file is a text file with the relative paths or
  // it is a JSON file from OpenMVG
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
      return;
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

    // Open file and fill the vector
    fs.open(fileFullPath, std::ios::in);

    if(!fs.is_open())
    {
      std::cerr << "Error while opening " << fileFullPath << std::endl;
      throw std::invalid_argument("Error while opening " + fileFullPath);
    }

    // count the name of files to load (ie the number of lines)
    auto numberOfFiles = std::count(std::istreambuf_iterator<char>(fs),
                                    std::istreambuf_iterator<char>(), '\n');

    if(numberOfFiles == 0)
    {
      std::cerr << "Could not found any file to load in " << fileFullPath <<"..." << std::endl;
      return;
    }

    // get back to the beginning of the file
    fs.seekg(0, std::ios::beg);
    descriptorsFiles.reserve(numberOfFiles);
    
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
    // Close and return
    fs.close();

  }
  else
  {
    std::cerr << "File not recognized! " << fileFullPath << std::endl;
    std::cerr << "The file  " + fileFullPath + " is neither a JSON nor a txt file" << std::endl;
    throw std::invalid_argument("Unrecognized file format " + fileFullPath);
  }
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