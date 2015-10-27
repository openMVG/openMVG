#include "SIFT_CCTAG_describer.hpp"

//#include <openMVG/features/descriptor.hpp>
//#include <openMVG/features/image_describer.hpp>
//#include <openMVG/features/regions_factory.hpp>
#include <openMVG/features/cctag/CCTAG_describer.hpp>

namespace openMVG {
namespace features {

  bool SIFT_CCTAG_Image_describer::Describe(const image::Image<unsigned char>& image,
    std::unique_ptr<Regions> &regions,
    const image::Image<unsigned char> * mask)
  {
    // Sift description
    std::cout << "SIFT description" << std::endl;
    features::SIFT_Image_describer siftDescriber(_paramsSift, _bOrientation);
    siftDescriber.Describe(image, regions);
    
    // CCTag description
    std::cout << "CCTag description" << std::endl;
    features::CCTAG_Image_describer cctagDescriber(_paramsCCTag._nCrowns, true);
    cctagDescriber.Describe(image, regions); 

    return true;
  };

} // features
} // openMVG