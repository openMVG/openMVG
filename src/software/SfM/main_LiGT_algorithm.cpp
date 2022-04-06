// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2022, Qi Cai and Yuanxin Wu

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/multiview/LiGT/LiGT_algorithm_converter.hpp"
#include "third_party/cmdLine/cmdLine.h"

int main(int argc, char **argv)
{
    using namespace std;
    OPENMVG_LOG_INFO
            << "\n-----------------------------------------------------------"
            << "\n The LiGT algorithm (version 1.1):"
            << "\n-----------------------------------------------------------";
    CmdLine cmd;

    // Common options:
    std::string
            gR_file,
            tracks_file,
            output_file,
            time_file;

    IndexT fixed_id = 0;
#ifdef OPENMVG_USE_OPENMP
  int iNumThreads = 1;
#endif

    // Common options
    cmd.add( make_option('r', gR_file, "gR_file") );
    cmd.add( make_option('i', tracks_file, "tracks_file") );
    cmd.add( make_option('o', output_file, "output_file") );
    cmd.add( make_option('t', time_file, "time_file") );
    cmd.add( make_option('f', fixed_id, "fixed_id") );
#ifdef OPENMVG_USE_OPENMP
  cmd.add( make_option('n', iNumThreads, "numThreads") );
#endif

    try {
        if (argc == 1) throw std::string("Invalid parameter.");
        cmd.process(argc, argv);
    } catch (const std::string& s) {

        OPENMVG_LOG_ERROR << s;
        return EXIT_FAILURE;
    }

#ifdef OPENMVG_USE_OPENMP
    const unsigned int nb_max_thread = omp_get_max_threads();

    if (iNumThreads > 0) {
        omp_set_num_threads(iNumThreads);
    } else {
        omp_set_num_threads(nb_max_thread);
    }
OPENMVG_LOG_ERROR << "num threads = " << iNumThreads;
#endif

   {
        LiGT::LiGTProblem LXR_solve(gR_file,
                               tracks_file,
                               output_file,
                               time_file,
                               fixed_id);

        LXR_solve.Solution();
        LXR_solve.WriteTranslation();
        LXR_solve.WriteTime();
        return EXIT_SUCCESS;
    }
    return EXIT_FAILURE;
}
