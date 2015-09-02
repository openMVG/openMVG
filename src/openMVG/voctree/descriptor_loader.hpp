#pragma once

#include <openMVG/features/descriptor.hpp>
#include <openMVG/sfm/sfm_data_io.hpp>

#include <boost/filesystem.hpp>
#include <boost/algorithm/string/case_conv.hpp>
#include <boost/progress.hpp>

#include <iostream>
#include <fstream>
#include <string>

using namespace std;

namespace openMVG {
namespace voctree {

void getInfoBinFile( const std::string &path, int dim, size_t &numDescriptors, int &bytesPerElement );

template<class DescriptorT>
size_t readDescFromFiles( const std::string &fileFullPath, std::vector<DescriptorT>& descriptors, std::vector<size_t> &numFeatures )
{
	namespace boostfs = boost::filesystem;
	std::ifstream fs;
	boostfs::path pathToFiles;
	std::string line;

	size_t numDescriptors = 0;
	boostfs::path bp(fileFullPath);

	if( !bp.has_extension() )
	{
		cerr << "File without extension not recognized! " << fileFullPath << endl;
		cerr << "The file  " + fileFullPath + " is neither a JSON nor a txt file" << endl;
	}

	// get the extension of the file and put it lowercase
	std::string ext = bp.extension().string();
	boost::to_lower( ext );

	// two cases, either the input file is a text file with the relative paths or
	// it is a JSON file from OpenMVG

	// if it is a JSON file
	if( ext == ".json")
	{
		// processing a JSON file containing sfm_data

		// open the sfm_data file
		openMVG::sfm::SfM_Data sfmdata;
		openMVG::sfm::Load( sfmdata, fileFullPath, openMVG::sfm::ESfM_Data::VIEWS );

		// get the number of files to load
		size_t numberOfFiles = sfmdata.GetViews().size();

		if( numberOfFiles == 0 )
		{
			cout << "It seems like there are no views in " << fileFullPath << endl;
			return 0;
		}

		// get the base path for the files
		pathToFiles = boostfs::path( fileFullPath ).parent_path();

		// contains a the size in byte for each descriptor element
		// could be 1 for char/uchar, 4 for float
		int bytesPerElement = 0;

		cout << "Pre-computing the memory needed..." << endl;
		boost::progress_display display( numberOfFiles );
		// pass all the files to get the number of features to load
		for(openMVG::sfm::Views::const_iterator it = sfmdata.GetViews().begin(); it != sfmdata.GetViews().end(); ++it, ++display)
		{
			// get just the image name, remove the extension
			std::string filename = boostfs::path( it->second->s_Img_path ).stem().string();
			// and generate the equivalent .desc filename

			filename = boostfs::path( pathToFiles / boostfs::path(filename+".desc")).string();

			// if it is the first one read the number of descriptors and the type of data (we are assuming the the feat are all the same...)
			// bytesPerElement could be 0 even after the first element (eg it has 0 descriptors...), so do it until
			// we get the correct info
			if( bytesPerElement == 0 )
			{
				getInfoBinFile( filename, DescriptorT::static_size, numDescriptors, bytesPerElement );
			}
			else
			{
				// get the file size in byte and estimate the number of features without opening the file
				numDescriptors += ( boostfs::file_size( filename ) / bytesPerElement ) / DescriptorT::static_size;
			}
		}
		BOOST_ASSERT( bytesPerElement > 0 );
		cout << "Found " << numDescriptors << " descriptors overall, allocating memory..." << endl;


		descriptors.reserve( numDescriptors );
		size_t numDesc =  numDescriptors;  // this is just to check everything is ok after
		numDescriptors = 0;

		cout << "Reading the descriptors..." << endl;
		display.restart( numberOfFiles );
		// run through the poses and read the descriptors
		for(openMVG::sfm::Views::const_iterator it = sfmdata.GetViews().begin(); it != sfmdata.GetViews().end(); ++it, ++display )
		{
			// for the first one, read the descriptor and estimate the size of memory to reserve

			// get just the image name, remove the extension
			std::string filename = boostfs::path( it->second->s_Img_path ).stem().string();
			// and generate the equivalent .desc filename

			filename = boostfs::path( pathToFiles / boostfs::path(filename+".desc")).string();

			// read the descriptor and append them in the vector
			loadDescsFromBinFile(filename, descriptors, true);
			size_t result = descriptors.size();

			// add the number of extracted features for this file
      // so it is the new size minus the previous size
			numFeatures.push_back( result - numDescriptors );
			// update the overall counter
			numDescriptors = result;
		}

		BOOST_ASSERT( numDesc == numDescriptors );
		return numDescriptors;

	}
	else
	{
		cerr << "File not recognized! " << fileFullPath << endl;
		cerr << "The file  " + fileFullPath + " is neither a JSON nor a txt file" << endl;
	}
}


}
}
