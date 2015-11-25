/* 
 * File:   svgVisualization.cpp
 * Author: sgaspari
 * 
 * Created on October 19, 2015, 9:46 AM
 */

#include "svgVisualization.hpp"
#if HAVE_CCTAG
#include "CCTagLocalizer.hpp"
#include <openMVG/features/cctag/CCTAG_describer.hpp>
#endif
#include "third_party/vectorGraphics/svgDrawer.hpp"

namespace openMVG {
namespace localization {


void saveMatches2SVG(const std::string &imagePathLeft,
                     const std::pair<size_t,size_t> & imageSizeLeft,
                     const std::vector<features::PointFeature> &keypointsLeft,
                     const std::string &imagePathRight,
                     const std::pair<size_t,size_t> & imageSizeRight,
                     const std::vector<features::PointFeature> &keypointsRight,
                     const matching::IndMatches &matches,
                     const std::string &outputSVGPath)
{
  svg::svgDrawer svgStream( imageSizeLeft.first + imageSizeRight.first, std::max(imageSizeLeft.second, imageSizeRight.second));
  svgStream.drawImage(imagePathLeft, imageSizeLeft.first, imageSizeLeft.second);
  svgStream.drawImage(imagePathRight, imageSizeRight.first, imageSizeRight.second, imageSizeLeft.first);
  
  for(const matching::IndMatch &m : matches) 
  {
    //Get back linked feature, draw a circle and link them by a line
    const features::PointFeature & L = keypointsLeft[m._i];
    const features::PointFeature & R = keypointsRight[m._j];
    svgStream.drawLine(L.x(), L.y(), R.x()+imageSizeLeft.first, R.y(), svg::svgStyle().stroke("green", 2.0));
    svgStream.drawCircle(L.x(), L.y(), 3.0f, svg::svgStyle().stroke("yellow", 2.0));
    svgStream.drawCircle(R.x()+imageSizeLeft.first, R.y(), 3.0f, svg::svgStyle().stroke("yellow", 2.0));
  }
 
  std::ofstream svgFile( outputSVGPath.c_str() );
  svgFile << svgStream.closeSvgFile().str();
  svgFile.close();
}


void saveFeatures2SVG(const std::string &inputImagePath,
                      const std::pair<size_t,size_t> & imageSize,
                      const std::vector<features::PointFeature> &keypoints,
                      const std::string &outputSVGPath)
{
  svg::svgDrawer svgStream( imageSize.first, imageSize.second);
  svgStream.drawImage(inputImagePath, imageSize.first, imageSize.second);
  
  for(const features::PointFeature &kpt : keypoints) 
  {
    svgStream.drawCircle(kpt.x(), kpt.y(), 3.0f, svg::svgStyle().stroke("yellow", 2.0));
  }
 
  std::ofstream svgFile( outputSVGPath.c_str() );
  svgFile << svgStream.closeSvgFile().str();
  svgFile.close();
}

#if HAVE_CCTAG
void saveCCTag2SVG(const std::string &inputImagePath,
                      const std::pair<size_t,size_t> & imageSize,
                      const features::CCTAG_Regions &cctags,
                      const std::string &outputSVGPath)
{
  svg::svgDrawer svgStream( imageSize.first, imageSize.second);
  svgStream.drawImage(inputImagePath, imageSize.first, imageSize.second);
  

  
  const auto &feat = cctags.Features();
  const std::vector<CCTagDescriptor > &desc = cctags.Descriptors();
  
  for(std::size_t i = 0; i < desc.size(); ++i) 
  {
    const IndexT cctagId = features::getCCTagId(desc[i]);
    if ( cctagId == UndefinedIndexT)
    {
      continue;
    }
    const CCTagKeypoint &kpt = feat[i];
    svgStream.drawCircle(kpt.x(), kpt.y(), 3.0f, svg::svgStyle().stroke("yellow", 2.0));
    svgStream.drawText(kpt.x(), kpt.y(), 20, std::to_string(cctagId), "yellow");
  }
 
  std::ofstream svgFile( outputSVGPath.c_str() );
  svgFile << svgStream.closeSvgFile().str();
  svgFile.close();
}
#endif

} // namespace localization
} // namespace openMVG

