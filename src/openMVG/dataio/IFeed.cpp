#include "IFeed.hpp"

#include <fstream>
#include <exception>

namespace openMVG{
namespace dataio{

// the structure of the file is
// int #image width
// int #image height
// double #focal
// double #ppx principal point x-coord
// double #ppy principal point y-coord
// double #k0
// double #k1
// double #k2
void readCalibrationFromFile(const std::string &filename, cameras::Pinhole_Intrinsic_Radial_K3 &camIntrinsics)
{
  std::ifstream fs(filename, std::ios::in);
  if(!fs.is_open())
  {
    std::cerr << "Unable to open the calibration file " << filename << std::endl;
    throw std::invalid_argument("Unable to open the calibration file "+filename);
  }
  int width = 0;
  int height = 0;
  const size_t numParam = 6;
  std::vector<double> params(numParam, 0);
  
  fs >> width;
  fs >> height;
  for(size_t i = 0; i < numParam; ++i)
  {
    fs >> params[i];
  }
  camIntrinsics = cameras::Pinhole_Intrinsic_Radial_K3(width, height, 
                                  params[0], params[1], params[2],
                                  params[3], params[4], params[5]);
  
  fs.close();
}

}//namespace dataio 
}//namespace openMVG