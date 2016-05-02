/* 
 * File:   svgVisualization.cpp
 * Author: sgaspari
 * 
 * Created on October 19, 2015, 9:46 AM
 */

#include "svgVisualization.hpp"
#if HAVE_CCTAG
#include <openMVG/localization/CCTagLocalizer.hpp>
#include "cctag/CCTAG_describer.hpp"
#endif
#include "third_party/vectorGraphics/svgDrawer.hpp"

namespace openMVG {
namespace features {


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
 
  std::ofstream svgFile( outputSVGPath );
  svgFile << svgStream.closeSvgFile().str();
  svgFile.close();
}

/**
 * @brief
 * 
 * @param[in] inputImagePath The full path to the image file. The image is only 
 * saved as a link, no image data is stored in the svg.
 * @param[in] imageSize The size of the image <width,height>.
 * @param[in] points A vector containing the points to draw.
 * @param[in] outputSVGPath The name of the svg file to generate.
 */
void saveFeatures2SVG(const std::string &inputImagePath,
                      const std::pair<size_t,size_t> & imageSize,
                      const Mat &points,
                      const std::string &outputSVGPath,
                      const std::vector<size_t> *inliers /*=nullptr*/)
{
  assert(points.rows()>=2);
  svg::svgDrawer svgStream( imageSize.first, imageSize.second);
  svgStream.drawImage(inputImagePath, imageSize.first, imageSize.second);
  
  if(!inliers)
  {
    for(std::size_t i=0; i < points.cols(); ++i) 
    {
      svgStream.drawCircle(points(0,i), points(1,i), 3.0f, svg::svgStyle().stroke("yellow", 2.0));
    }
  }
  else
  {
    for(std::size_t i=0; i < points.cols(); ++i) 
    {
      if(std::find(inliers->begin(), inliers->end(), i) != inliers->end())
        svgStream.drawCircle(points(0,i), points(1,i), 3.0f, svg::svgStyle().stroke("green", 3.0));
      else
        svgStream.drawCircle(points(0,i), points(1,i), 3.0f, svg::svgStyle().stroke("yellow", 2.0));
    }   
  }
 
  std::ofstream svgFile( outputSVGPath );
  svgFile << svgStream.closeSvgFile().str();
  svgFile.close();
}

#if HAVE_CCTAG

void saveCCTag2SVG(const std::string &inputImagePath,
                      const std::pair<size_t,size_t> & imageSize,
                      const features::CCTAG_Regions &cctags,
                      const std::string &outputSVGPath)
{
  // set the text size to 5% if the image heigth
  const float textSize = 0.05*imageSize.second;
  
  svg::svgDrawer svgStream( imageSize.first, imageSize.second);
  svgStream.drawImage(inputImagePath, imageSize.first, imageSize.second);
    
  const auto &feat = cctags.Features();
  const std::vector<localization::CCTagDescriptor > &desc = cctags.Descriptors();
  
  for(std::size_t i = 0; i < desc.size(); ++i) 
  {
    const IndexT cctagId = features::getCCTagId(desc[i]);
    if ( cctagId == UndefinedIndexT)
    {
      continue;
    }
    const localization::CCTagKeypoint &kpt = feat[i];
    svgStream.drawCircle(kpt.x(), kpt.y(), 3.0f, svg::svgStyle().stroke("yellow", 2.0));
    svgStream.drawText(kpt.x(), kpt.y(), textSize, std::to_string(cctagId), "yellow");
  }
 
  std::ofstream svgFile( outputSVGPath );
  svgFile << svgStream.closeSvgFile().str();
  svgFile.close();
}

void saveCCTagMatches2SVG(const std::string &imagePathLeft,
                     const std::pair<size_t,size_t> & imageSizeLeft,
                     const features::CCTAG_Regions &cctagLeft,
                     const std::string &imagePathRight,
                     const std::pair<size_t,size_t> & imageSizeRight,
                     const features::CCTAG_Regions &cctagRight,
                     const matching::IndMatches &matches,
                     const std::string &outputSVGPath, 
                     bool showNotMatched)
{
  // set the text size to 5% if the image heigth
  const float textSize = 0.05*std::min(imageSizeRight.second, imageSizeLeft.second);
  
  svg::svgDrawer svgStream(imageSizeLeft.first + imageSizeRight.first, std::max(imageSizeLeft.second, imageSizeRight.second));
  svgStream.drawImage(imagePathLeft, imageSizeLeft.first, imageSizeLeft.second);
  svgStream.drawImage(imagePathRight, imageSizeRight.first, imageSizeRight.second, imageSizeLeft.first);

  const auto &keypointsLeft = cctagLeft.Features();
  const auto &keypointsRight = cctagRight.Features();
  const std::vector<localization::CCTagDescriptor > &descLeft = cctagLeft.Descriptors();
  const std::vector<localization::CCTagDescriptor > &descRight = cctagRight.Descriptors();
  
  //just to be sure...
  assert(keypointsLeft.size() == descLeft.size());
  assert(keypointsRight.size() == descRight.size());
  
  for(const matching::IndMatch &m : matches)
  {
    //Get back linked feature, draw a circle and link them by a line
    const features::PointFeature & L = keypointsLeft[m._i];
    const features::PointFeature & R = keypointsRight[m._j];
    const IndexT cctagIdLeft = features::getCCTagId(descLeft[m._i]);
    const IndexT cctagIdRight = features::getCCTagId(descRight[m._j]);
    if ( cctagIdLeft == UndefinedIndexT || cctagIdRight == UndefinedIndexT )
    {
      std::cerr << "[svg]\tWarning! cctagIdLeft " << cctagIdLeft << " " << "cctagIdRight " << cctagIdRight << std::endl;
      continue;
    }
    
    svgStream.drawLine(L.x(), L.y(), R.x() + imageSizeLeft.first, R.y(), svg::svgStyle().stroke("green", 2.0));
    svgStream.drawCircle(L.x(), L.y(), 3.0f, svg::svgStyle().stroke("yellow", 2.0));
    svgStream.drawCircle(R.x() + imageSizeLeft.first, R.y(), 3.0f, svg::svgStyle().stroke("yellow", 2.0));
    
    svgStream.drawText(L.x(), L.y(), textSize, std::to_string(cctagIdLeft), "yellow");
    svgStream.drawText(R.x() + imageSizeLeft.first, R.y(), textSize, std::to_string(cctagIdRight), "yellow");
  }
  
  if(showNotMatched)
  {
    const float textSizeSmaller = 0.75*textSize;
    
    // for the left side
    for(std::size_t i = 0; i < keypointsLeft.size(); ++i)
    {
      // look if the index is not already in matches
      bool found = false;
      for(const matching::IndMatch &m : matches)
      {
        found = m._i == i;
        if(found)
          break;
      }
      if(found)
        continue;
      
      assert(i < descLeft.size());
      // find the cctag id
      const IndexT cctagIdLeft = features::getCCTagId(descLeft[i]);
      if(cctagIdLeft == UndefinedIndexT)
        continue;
      
      // draw the center
      const features::PointFeature & L = keypointsLeft[i];
      svgStream.drawCircle(L.x(), L.y(), 3.0f, svg::svgStyle().stroke("red", 2.0));
      // print the id
      svgStream.drawText(L.x(), L.y(), textSizeSmaller, std::to_string(cctagIdLeft), "red");
    }
    
    // for the right side
    for(std::size_t i = 0; i < keypointsRight.size(); ++i)
    {
      // look if the index is not already in matches
      bool found = false;
      for(const matching::IndMatch &m : matches)
      {
        found = m._j == i;
        if(found)
          break;
      }
      if(found)
        continue;
      
      assert(i < descRight.size());
      // find the cctag id
      const IndexT cctagIdRight = features::getCCTagId(descRight[i]);
      if(cctagIdRight == UndefinedIndexT)
        continue;
      
      // draw the center
      const features::PointFeature & R = keypointsRight[i];
      svgStream.drawCircle(R.x() + imageSizeLeft.first, R.y(), 3.0f, svg::svgStyle().stroke("red", 2.0));
      // print the id
      svgStream.drawText(R.x() + imageSizeLeft.first, R.y(), textSizeSmaller, std::to_string(cctagIdRight), "red");

    }
  }

  std::ofstream svgFile(outputSVGPath.c_str());
  svgFile << svgStream.closeSvgFile().str();
  svgFile.close();
}
#endif

} // namespace features
} // namespace openMVG

