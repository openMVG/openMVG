//:\file
//\author Ricardo Fabbri, Brown & Rio de Janeiro State U. (rfabbri.github.io) 
//\date Tue Jun  1 09:04:21 -03 2021
//\author Gabriel ANDRADE Rio de Janeiro State U.
//\author Pierre MOULON
#include "trifocal-app.h"

using namespace trifocal3pt;
//-------------------------------------------------------------------------------
int 
main(int argc, char **argv) 
{
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
