

// Copyright (c) 2012, 2013, 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/sfm/sfm.hpp"
#include "openMVG/image/image.hpp"

using namespace openMVG;

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"
#include "openMVG/numeric/numeric.h"

#include <fstream>

int main(int argc, char **argv)
{
  CmdLine cmd;

  std::string sSfM_Data_Filename;
  std::string sOutDir = "";

  cmd.add( make_option('i', sSfM_Data_Filename, "sfmdata") );
  cmd.add( make_option('o', sOutDir, "outdir") );

  try {
      if (argc == 1) throw std::string("Invalid command line parameter.");
      cmd.process(argc, argv);
  } catch(const std::string& s) {
      std::cerr << "Usage: " << argv[0] << '\n'
      << "[-i|--sfmdata filename, the SfM_Data file to convert]\n"
      << "[-o|--outdir path]\n"
      << std::endl;

      std::cerr << s << std::endl;
      return EXIT_FAILURE;
  }

  std::cout << " You called : " <<std::endl
            << argv[0] << std::endl
            << "--sfmdata " << sSfM_Data_Filename << std::endl
            << "--outdir " << sOutDir << std::endl;

  // Create output dir
  if (!stlplus::folder_exists(sOutDir))
    stlplus::folder_create( sOutDir );

  // Read the SfM scene
  SfM_Data sfm_data;
  if (!Load(sfm_data, sSfM_Data_Filename, ESfM_Data(VIEWS|INTRINSICS|EXTRINSICS))) {
    std::cerr << std::endl
      << "The input SfM_Data file \""<< sSfM_Data_Filename << "\" cannot be read." << std::endl;
    return EXIT_FAILURE;
  }


  for(Views::const_iterator iter = sfm_data.getViews().begin();
      iter != sfm_data.getViews().end(); ++iter)
  {
    const View * view = iter->second.get();
    std::cout << view->s_Img_path << "\n";
    std::cout << stlplus::basename_part(view->s_Img_path) << std::endl;
    std::ofstream outfile( stlplus::create_filespec(
                sOutDir, stlplus::basename_part(view->s_Img_path), "cam" ).c_str() );
    Poses::const_iterator iterPose = sfm_data.getPoses().find(view->id_pose);
    Intrinsics::const_iterator iterIntrinsic = sfm_data.getIntrinsics().find(view->id_intrinsic);

    if (iterPose == sfm_data.getPoses().end() ||
      iterIntrinsic == sfm_data.getIntrinsics().end())
        continue;

    const std::string srcImage = stlplus::create_filespec(sfm_data.s_root_path, view->s_Img_path);
    const IntrinsicBase * cam = iterIntrinsic->second.get();
    Mat34 P = cam->get_projective_equivalent(iterPose->second);

    Mat3 R, K;
    Vec3 t;
    KRt_From_P(P, &K, &R, &t);

    int largerDim = cam->w() > cam->h() ? cam->w() : cam->h();

    // see https://github.com/nmoehrle/mvs-texturing/blob/master/Arguments.cpp
    // for full specs
    outfile << t(0) << " " << t(1) << " " << t(2) << " "
        << R(0,0) << " " << R(0,1) << " " << R(0,2) << " "
        << R(1,0) << " " << R(1,1) << " " << R(1,2) << " "
        << R(2,0) << " " << R(2,1) << " " << R(2,2) << "\n"
        << K(0,0) / largerDim << " 0 0 1 " << K(0,2) / cam->w() << " " << K(1,2) / cam->h();
    //TODO : undist
    outfile.close();

  }
    return EXIT_SUCCESS;
}
