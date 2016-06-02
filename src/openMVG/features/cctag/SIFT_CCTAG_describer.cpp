#ifdef HAVE_CCTAG

#include "SIFT_CCTAG_describer.hpp"

//#include <openMVG/features/descriptor.hpp>
#include <openMVG/features/image_describer.hpp>
#include <openMVG/features/regions_factory.hpp>
#include <openMVG/features/regions.hpp>
#include <openMVG/features/cctag/CCTAG_describer.hpp>

namespace openMVG {
namespace features {

SIFT_CCTAG_Image_describer::SIFT_CCTAG_Image_describer(const SiftParams & params, bool bOrientation, std::size_t nRings)
  :Image_describer(), _siftDescriber(params, bOrientation), _cctagDescriber(nRings)
{
}

bool SIFT_CCTAG_Image_describer::Set_configuration_preset(EDESCRIBER_PRESET preset)
{
  bool res = _siftDescriber.Set_configuration_preset(preset);
  res &= _cctagDescriber.Set_configuration_preset(preset);
  return res;
}

bool SIFT_CCTAG_Image_describer::Describe(const image::Image<unsigned char>& image,
  std::unique_ptr<Regions> &regions,
  const image::Image<unsigned char> * mask)
{
  std::unique_ptr<Regions> siftRegionsPtr;
  // Sift description
  std::cout << "SIFT description" << std::endl;
  _siftDescriber.Describe(image, siftRegionsPtr);
  SIFT_Regions* siftRegions = (SIFT_Regions*) siftRegionsPtr.get();

  // CCTag description
  std::cout << "CCTag description" << std::endl;
  _cctagDescriber.Describe(image, regions);
  
  SIFT_Regions* outRegions = (SIFT_Regions*) regions.get();
  std::size_t nbCCTags = outRegions->Features().size();

  outRegions->Features().resize(nbCCTags + siftRegions->Features().size());
  std::copy(siftRegions->Features().begin(), siftRegions->Features().end(), outRegions->Features().begin()+nbCCTags);
  
  outRegions->Descriptors().resize(nbCCTags + siftRegions->Descriptors().size());
  std::copy(siftRegions->Descriptors().begin(), siftRegions->Descriptors().end(), outRegions->Descriptors().begin()+nbCCTags);

  return true;
};

} // features
} // openMVG

#endif //HAVE_CCTAG