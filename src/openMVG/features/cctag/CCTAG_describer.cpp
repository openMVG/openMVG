#include "CCTAG_describer.hpp"

//#define CPU_ADAPT_OF_GPU_PART //todo: #ifdef depreciated
#ifdef CPU_ADAPT_OF_GPU_PART    
  #include "cctag/progBase/MemoryPool.hpp"
#endif

namespace openMVG {
namespace features {

bool CCTAG_Image_describer::Describe(const image::Image<unsigned char>& image,
    std::unique_ptr<Regions> &regions,
    const image::Image<unsigned char> * mask)
  {
    const int w = image.Width(), h = image.Height();
    
    Allocate(regions);

    // Build alias to cached data
    CCTAG_Regions * regionsCasted = dynamic_cast<CCTAG_Regions*>(regions.get());
    // reserve some memory for faster keypoint saving
    regionsCasted->Features().reserve(50);
    regionsCasted->Descriptors().reserve(50);
    
    boost::ptr_list<cctag::ICCTag> cctags;
    
#ifndef CPU_ADAPT_OF_GPU_PART
    const cv::Mat graySrc(cv::Size(image.Width(), image.Height()), CV_8UC1, (unsigned char *) image.data(), cv::Mat::AUTO_STEP);
    //// Invert the image
    //cv::Mat invertImg;
    //cv::bitwise_not(graySrc,invertImg);
    cctag::cctagDetection(cctags,1,graySrc,_params);
#else //todo: #ifdef depreciated
    cctag::MemoryPool::instance().updateMemoryAuthorizedWithRAM();
    cctag::View cctagView((const unsigned char *) image.data(), image.Width(), image.Height(), image.Depth()*image.Width());
    boost::ptr_list<cctag::ICCTag> cctags;
    cctag::cctagDetection(cctags, 1 ,cctagView._grayView ,_params);
#endif
    
    for (const auto & cctag : cctags)
    {
      if ( cctag.getStatus() > 0 )
      {
        std::cout << " New CCTag: " << cctag.id() << " ( " << cctag.x() << " , " << cctag.y() << " ) " << std::endl;

        // Add its associated descriptor
        Descriptor<unsigned char,128> desc;
        for(int i=0; i< desc.size(); ++i)
        {
          desc[i] = 0.f;
        }
        desc[cctag.id()] = 1.f;
        regionsCasted->Descriptors().push_back(desc);
        regionsCasted->Features().push_back(SIOPointFeature(cctag.x(), cctag.y()));
      }
    }

    cctags.clear();

    return true;
  };

} // namespace features
} // namespace openMVG