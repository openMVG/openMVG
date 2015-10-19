#ifndef OPENMVG_CCTAG_DESCRIBER_H
#define OPENMVG_CCTAG_DESCRIBER_H

//#include <cctag/view.hpp>
#include <cctag/ICCTag.hpp>

#include <cereal/cereal.hpp>
#include <iostream>
#include <numeric>

namespace openMVG {
namespace features {

/** @brief CCTag filter pixel type */

class CCTAG_Image_describer : public Image_describer
{
public:
  CCTAG_Image_describer()
    :Image_describer(), _nRings(3) {}
    
  CCTAG_Image_describer(const std::size_t nRings)
    :Image_describer(), _nRings(nRings) {}   

  ~CCTAG_Image_describer() {}

  bool Set_configuration_preset(EDESCRIBER_PRESET preset)
  {
//    switch(preset)
//    {
//    case NORMAL_PRESET:
//      _params._peak_threshold = 0.04f;
//    break;
//    case HIGH_PRESET:
//      _params._peak_threshold = 0.01f;
//    break;
//    case ULTRA_PRESET:
//      _params._peak_threshold = 0.01f;
//      _params._first_octave = -1;
//    break;
//    default:
//      return false;
//    }
    return true;
  }

  /**
  @brief Detect regions on the image and compute their attributes (description)
  @param image Image.
  @param regions The detected regions and attributes (the caller must delete the allocated data)
  @param mask 8-bit gray image for keypoint filtering (optional).
     Non-zero values depict the region of interest.
  */
  bool Describe(const image::Image<unsigned char>& image,
    std::unique_ptr<Regions> &regions,
    const image::Image<unsigned char> * mask = NULL)
  {
    const int w = image.Width(), h = image.Height();

    Allocate(regions);

    // Build alias to cached data
    CCTAG_Regions * regionsCasted = dynamic_cast<CCTAG_Regions*>(regions.get());
    // reserve some memory for faster keypoint saving
    regionsCasted->Features().reserve(50);
    regionsCasted->Descriptors().reserve(50);
    
    const cv::Mat graySrc(cv::Size(image.Width(), image.Height()), CV_8UC1, (unsigned char *) image.data(), cv::Mat::AUTO_STEP);
    
    boost::ptr_list<cctag::ICCTag> cctags;
    
    cctag::cctagDetection(cctags,1,graySrc,_nRings);
    
    for (const auto & cctag : cctags)
    {
      if ( cctag.getStatus() > 0 )
      {
        std::cout << " New CCTag: " << cctag.id() << " ( " << cctag.x() << " , " << cctag.y() << " ) " << std::endl;

        // Add its associated descriptor
        Descriptor<unsigned char,128> desc;
        for(int i=0; i< desc.size(); ++i)
        {
          desc[i] = 0;
        }
        desc[cctag.id()] = 255;
        regionsCasted->Descriptors().push_back(desc);
        regionsCasted->Features().push_back(SIOPointFeature(cctag.x(), cctag.y()));
      }
    }

    cctags.clear();

    return true;
  };

  /// Allocate Regions type depending of the Image_describer
  void Allocate(std::unique_ptr<Regions> &regions) const
  {
    regions.reset( new CCTAG_Regions );
  }

  template<class Archive>
  void serialize( Archive & ar )
  {
    //ar(
     //cereal::make_nvp("params", _params),
     //cereal::make_nvp("bOrientation", _bOrientation));
  }

private:
  //SiftParams _params;
  //bool _bOrientation;
  std::size_t _nRings;
};

} // namespace features
} // namespace openMVG

#include <cereal/types/polymorphic.hpp>
#include <cereal/archives/json.hpp>
CEREAL_REGISTER_TYPE_WITH_NAME(openMVG::features::CCTAG_Image_describer, "CCTAG_Image_describer");

#endif // OPENMVG_CCTAG_DESCRIBER_H
