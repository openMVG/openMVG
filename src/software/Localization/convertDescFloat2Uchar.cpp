/* 
 * File:   convertDescFloat2Uchar.cpp
 * Author: sgaspari
 *
 * Created on November 28, 2015, 1:30 PM
 */
#include <openMVG/features/descriptor.hpp>
#include <openMVG/logger.hpp>

#include <boost/progress.hpp>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp> 
#include <boost/algorithm/string/case_conv.hpp> 

#include <cstdlib>


namespace bfs = boost::filesystem;
namespace po = boost::program_options;

using namespace openMVG;
/*
 * 
 */
int main( int argc, char** argv )
{
  std::string outputFolder;
  std::string inputFolder;
  bool doSanityCheck = false;
  const int siftSize = 128;

  po::options_description desc("This program is used to convert SIFT features from float representation to unsigned char representation");
  desc.add_options()
        ("inputFolder,i", po::value<std::string>(&inputFolder)->required(), "Input folder containing the sift in float format.")
        ("outputFolder,o", po::value<std::string>(&outputFolder)->required(), "Output folder that stores the sift in uchar format.")
        ("sanityCheck,s", po::bool_switch(&doSanityCheck)->default_value(doSanityCheck), "Perform a sanity check to check that the conversion and the genrated files are the same.");

  po::variables_map vm;

  try
  {
    po::store(po::parse_command_line(argc, argv, desc), vm);

    if(vm.count("help") || (argc == 1))
    {
      POPART_COUT(desc);
      return EXIT_SUCCESS;
    }

    po::notify(vm);
  }
  catch(boost::program_options::required_option& e)
  {
    POPART_CERR("ERROR: " << e.what() << std::endl);
    POPART_COUT("Usage:\n\n" << desc);
    return EXIT_FAILURE;
  }
  catch(boost::program_options::error& e)
  {
    POPART_CERR("ERROR: " << e.what() << std::endl);
    POPART_COUT("Usage:\n\n" << desc);
    return EXIT_FAILURE;
  }

  if(!(bfs::exists(inputFolder) && bfs::is_directory(inputFolder)))
  {
    POPART_CERR("ERROR: " << inputFolder << " does not exists or it is not a directory" << std::endl);
    return EXIT_FAILURE;
  }

  // if the directory does not exist create it (recursively)
  if(!bfs::exists(outputFolder))
  {
    bfs::create_directories(outputFolder);
  }
  
  size_t countFeat = 0;
  size_t countDesc = 0;

  bfs::directory_iterator iterator(inputFolder);
  for(; iterator != bfs::directory_iterator(); ++iterator)
  {
    // get the extension of the current file to check whether it is an image
    std::string ext = iterator->path().extension().string();
    const std::string &filename = iterator->path().filename().string();
    boost::to_lower(ext);

    if(ext == ".feat")
    {
      // just copy the file into the output directory
      bfs::copy_file(iterator->path(), bfs::path(outputFolder)/bfs::path(filename), bfs::copy_option::overwrite_if_exists);
      
      ++countFeat;
    }
    else if(ext == ".desc")
    {
      const std::string outpath = (bfs::path(outputFolder)/bfs::path(filename)).string(); 
      std::vector<features::Descriptor<float, siftSize> > floatDescriptors;
      
      // load the float descriptors
      features::loadDescsFromBinFile(iterator->path().string(), floatDescriptors, false);
      
      const size_t numDesc = floatDescriptors.size();
      
      std::vector<features::Descriptor<unsigned char, siftSize> > charDescriptors(numDesc);
 
      for(std::size_t i = 0; i < numDesc; ++i)
      {
        float* fptr = floatDescriptors[i].getData();
        assert(fptr!=nullptr);
        unsigned char* uptr = charDescriptors[i].getData();
        assert(uptr!=nullptr);
      
        std::copy(fptr, fptr+siftSize, uptr);
      
        if(!doSanityCheck)
          continue;    
        // check that they are actually the same
        for(std::size_t j = 0; j < siftSize; ++j)
        {      
          const unsigned char compare = (unsigned char) fptr[j];
          assert(compare == uptr[j]);
        }
      }    
      
      assert(charDescriptors.size() == floatDescriptors.size());
      
      // save the unsigned char
      features::saveDescsToBinFile(outpath, charDescriptors);
      
      if(doSanityCheck)
      {
        // sanity check 
        // reload everything and compare
        floatDescriptors.clear();
        charDescriptors.clear();
        features::loadDescsFromBinFile(iterator->path().string(), floatDescriptors, false);
        features::loadDescsFromBinFile(outpath, charDescriptors, false);

        assert(charDescriptors.size() == numDesc);
        assert(charDescriptors.size() == floatDescriptors.size());

        for(std::size_t i = 0; i < numDesc; ++i)
        {
          const features::Descriptor<float, siftSize> &currFloat = floatDescriptors[i];
          const features::Descriptor<unsigned char, siftSize> &currUchar = charDescriptors[i];
          for(std::size_t j = 0; j < siftSize; ++j)
          {
            const unsigned char compare = (unsigned char) currFloat[j];
            assert(compare == currUchar[j]);
          }
        }
      }
      ++countDesc;
    }
  }
  POPART_COUT("Converted " << countDesc << " files .desc and copied " << countFeat << " files .feat");
}
