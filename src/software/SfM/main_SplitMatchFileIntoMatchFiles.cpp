// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013, 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/sfm/pipelines/sfm_matches_provider.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"
#include "openMVG/sfm/sfm_data_graph_utils.hpp"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <cstdlib>
#include <memory>
#include <string>

using namespace openMVG;
using namespace openMVG::sfm;

int main(int argc, char **argv)
{
	using namespace std;
	std::cout << std::endl;

	CmdLine cmd;

	std::string sSfM_Data_Filename;
	std::string sMatchsFile;
	std::string sMatchFileComponents;
	bool bBiEdge = false;
	int nMinNodes = 3;

	cmd.add(make_option('i', sSfM_Data_Filename, "input_file"));
	cmd.add(make_option('m', sMatchsFile, "match_file"));
	cmd.add(make_option('o', sMatchFileComponents, "match_file_components"));
	cmd.add(make_option('b', bBiEdge, "biedge"));
	cmd.add(make_option('n', nMinNodes, "min_nodes"));

	try {
		if (argc == 1) throw std::string("Invalid parameter.");
		cmd.process(argc, argv);
	}
	catch (const std::string& s) {
		std::cerr << "Usage: " << argv[0] << '\n'
			<< "[-i|--input_file] path to a SfM_Data scene\n"
			<< "[-m|--match_file] path to the matches that corresponds to the provided SfM_Data scene\n"
			<< "[-o|--match_file_components] path to the matches components that corresponds to the provided SfM_Data scene\n"
			<< "\n[Optional]\n"
			<< "[-b|--biedge]\n"
			<< "[-n|--min_nodes]\n"
			<< "value of n should larger than 3\n" 
			<< std::endl;
		std::cerr << s << std::endl;
		return EXIT_FAILURE;
	}

	if (!stlplus::file_exists(sSfM_Data_Filename))
	{
		return EXIT_FAILURE;
	}

	if (!stlplus::file_exists(sMatchsFile))
	{
		return EXIT_FAILURE;
	}

	const std::string &sMatchFileComponentsDir = stlplus::folder_part(sMatchFileComponents);
	if (!stlplus::folder_exists(sMatchFileComponentsDir))
	{
		stlplus::folder_create(sMatchFileComponentsDir);
	}

	// Load input SfM_Data scene
	SfM_Data sfm_data;
	if (!Load(sfm_data, sSfM_Data_Filename, ESfM_Data(VIEWS | INTRINSICS)))
	{
		std::cerr << std::endl
			<< "The input SfM_Data file \"" << sSfM_Data_Filename << "\" cannot be read." << std::endl;
		return EXIT_FAILURE;
	}

	// Matches reading
	std::shared_ptr<Matches_Provider> matches_provider = std::make_shared<Matches_Provider>();
	if (
		!(matches_provider->load(sfm_data, sMatchsFile))
		)
	{
		std::cerr << std::endl
			<< "Invalid matches file." << std::endl;
		return EXIT_FAILURE;
	}

	// split match_file by connected compents;
	bool bFlag = SplitMatchFileIntoMatchFiles(sfm_data, sMatchsFile, sMatchFileComponents, bBiEdge, nMinNodes);
	if (bFlag)
	{
		return EXIT_SUCCESS;
	}
	else
	{
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}