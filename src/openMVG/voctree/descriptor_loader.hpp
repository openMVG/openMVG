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

/**
 * Get the number of descriptors contained inside a .desc file and the number of bytes
 * used to store each descriptor elements
 *
 * @param[in] path The .desc filename
 * @param[in] dim The number of elements per descriptor
 * @param[out] numDescriptors The number of descriptors stored in the file
 * @param[out] bytesPerElement The number of bytes used to store each element of the descriptor
 */
void getInfoBinFile( const std::string &path, int dim, size_t &numDescriptors, int &bytesPerElement );

/**
 * Read a set of descriptors from a file containing the path to the descriptor files.
 * Two format are supported:
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
 * @param[in] fileFullPath the input filename containing the list of file to read
 * @param[in,out] descriptors the vector to which append all the read descriptors
 * @param[in,out] numFeatures a vector collecting for each file read the number of features read
 * @return the total number of features read
 *
 */
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
	else if( ext == ".txt" )
	{

		// processing a file .txt containing the relative paths

		// Extract the folder path from the list file path
		pathToFiles = boostfs::path( fileFullPath ).parent_path();

		// Open file and fill the vector
		fs.open(fileFullPath, std::ios::in);

		if( !fs.is_open() )
		{
			cerr << "Error while opening " << fileFullPath << endl;
			cerr << "Error while opening " + fileFullPath << endl;
		}

		// count the name of files to load (ie the number of lines)
		auto numberOfFiles = std::count( std::istreambuf_iterator<char>( fs ),
									std::istreambuf_iterator<char>(), '\n' );

		if( numberOfFiles == 0 )
		{
			cout << "Could not found any file to load..." << endl;
			return 0;
		}

		// reserve space to append the number of extracted features
		numFeatures.reserve( numFeatures.size() + numberOfFiles );

		// get back to the beginning of the file
		fs.seekg (0, std::ios::beg);

		// contains a the size in byte for each descriptor element
		// could be 1 for char/uchar, 4 for float
		int bytesPerElement = 0;

		// pass all the files to get the number of features to load
		while( getline( fs, line ) )
		{
			// we have to do that because OMVG does not really output a clean list.txt, it also
			// contains other stuff, so we look at the first '.' to get the extension (not robust at all)
			std::string filename = line.substr( 0, line.find_first_of(".") );
			filename = boostfs::path( pathToFiles / boostfs::path(filename+".desc")).string();
			std::cout << filename << std::endl;

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

		descriptors.reserve( numDescriptors );
		size_t numDesc =  numDescriptors;  // this is just to check everything is ok after
		numDescriptors = 0;

		// get back to the beginning of the file
		fs.clear();
		fs.seekg (0, std::ios::beg);

		boost::progress_display display( numberOfFiles );
		while( getline( fs, line ) )
		{
			std::string filename = line.substr( 0, line.find_first_of(".") );
			filename = boostfs::path( pathToFiles / boostfs::path(filename+".desc")).string();
			loadDescsFromBinFile(filename, descriptors, true);
			size_t result = descriptors.size();
			// add the number of extracted features for this file
			numFeatures.push_back( result - numDescriptors );
			// update the overall counter
			numDescriptors = result;
			++display;
		}

		// Close and return
		fs.close();

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
