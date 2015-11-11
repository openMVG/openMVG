/* 
 * File:   svgVisualization.cpp
 * Author: sgaspari
 * 
 * Created on October 19, 2015, 9:46 AM
 */

#include "svgVisualization.hpp"
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

} // namespace localization
} // namespace openMVG

