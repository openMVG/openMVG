// \author Ricardo FABBRI
// \date Tue Jun  1 09:04:21 -03 2021
// \author Gabriel ANDRADE \author Pierre MOULON

#include "trifocal.h"
#include "trifocal_app.h"

void TrifocalSampleApp::
ProcessCmdLine(int argc, char **argv)
{
  CmdLine cmd;
  cmd.add( make_option('a', image_filenames_[0], "image_a") );
  cmd.add( make_option('b', image_filenames_[1], "image_b") );
  cmd.add( make_option('c', image_filenames_[2], "image_c") );
  cmd.add( make_option('K', intrinsics_filename_, "K matrix") );
  
  try {
    if (argc == 1) throw string("Invalid command line parameter.");
    cmd.process(argc, argv);
  } catch (const string& s) {
    cerr << "Usage: " << argv[0] << '\n' << endl;
    cerr << s << endl;
    exit(EXIT_FAILURE);
  }
}


//-------------------------------------------------------------------------------
int main(int argc, char **argv) {
  TrifocalSampleApp T;
  
  T.ProcessCmdLine(argc, argv);
  T.ExtractKeypoints();
  T.MatchKeypoints();
  T.ComputeTracks();
  T.Stats();
  T.ExtractXYOrientation();
  T.Display();
  T.DisplayDesiredIds();
  T.RobustSolve();
  // T.DisplayInliers();
  T.DisplayInliersCamerasAndPoints();
  // T.DisplayInliersCamerasAndPointsSIFT();

  return EXIT_SUCCESS;
}
