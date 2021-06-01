// \author Ricardo FABBRI
// \date Tue Jun  1 09:04:21 -03 2021
// \author Gabriel ANDRADE \author Pierre MOULON

#include "trifocal.h"
#include "trifocal_app.h"

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
