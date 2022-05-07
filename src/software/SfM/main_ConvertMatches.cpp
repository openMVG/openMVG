// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013, 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/matching/indMatch.hpp"
#include "openMVG/matching/indMatch_utils.hpp"

#include "third_party/cmdLine/cmdLine.h"

#include <cstdlib>
#include <string>

using namespace openMVG;
using namespace openMVG::matching;

int main(int argc, char ** argv)
{
  CmdLine cmd;

  std::string sMatchFile;
  std::string sOutMatchFile;

  cmd.add( make_option('m', sMatchFile, "matchfile") );
  cmd.add( make_option('o', sOutMatchFile, "outmatchfile") );

  try {
      if (argc == 1) throw std::string("Invalid command line parameter.");
      cmd.process(argc, argv);
  } catch (const std::string& s) {
      std::cerr << "Convert matches from bin to txt or txt to bin.\nUsage: " << argv[0] << "\n"
      << "[-m|--sMatchFile filename]\n"
      << "[-o|--outmatchfile filename]\n"
      << std::endl;

      std::cerr << s << std::endl;
      return EXIT_FAILURE;
  }

  if (sMatchFile.empty()) {
    std::cerr << "\nmatchfile cannot be an empty option" << std::endl;
    return EXIT_FAILURE;
  }
  if (sOutMatchFile.empty())  {
    std::cerr << "\noutmatchfile cannot be an empty option" << std::endl;
    return EXIT_FAILURE;
  }

  // Read the matches
  matching::PairWiseMatches pairWise_matches;
  if (!Load(pairWise_matches, sMatchFile)) {
    std::cerr << "\nInvalid matches file." << std::endl;
    return EXIT_FAILURE;
  }

  // Write the matches
  if (!Save(pairWise_matches, sOutMatchFile))
  {
    std::cerr
        << "Cannot save matches to: "
        << sOutMatchFile;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
