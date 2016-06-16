#include "svgVisualization.hpp"
#if HAVE_CCTAG
#include "cctag/CCTAG_describer.hpp"
#endif
#include "third_party/vectorGraphics/svgDrawer.hpp"

namespace openMVG {
namespace features {

float getRadiusEstimate(const std::pair<size_t,size_t> & imgSize)
{
  // heuristic for the radius of the feature according to the size of the image
  // the larger the distance the larger the minimum radius should be in order to be visible
  // We consider a minimum of 2 pixel and we increment it linearly according to 
  // the image size
  return std::max(std::max(imgSize.first, imgSize.second) / float(600), 2.0f);
}

float getStrokeEstimate(const std::pair<size_t,size_t> & imgSize)
{
  return std::max(std::max(imgSize.first, imgSize.second) / float(2200), 2.0f);
}

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
  
  // heuristic for the radius of the feature according to the size of the image
  // the larger the distance the larger the minimum radius should be in order to be visible
  // We consider a minimum of 2 pixel and we increment it linearly according to 
  // the image size
  const float radiusLeft = getRadiusEstimate(imageSizeLeft);
  const float radiusRight = getRadiusEstimate(imageSizeRight);
  const float strokeLeft = getStrokeEstimate(imageSizeLeft);
  const float strokeRight = getStrokeEstimate(imageSizeRight);
 
  
  for(const matching::IndMatch &m : matches) 
  {
    //Get back linked feature, draw a circle and link them by a line
    const features::PointFeature & L = keypointsLeft[m._i];
    const features::PointFeature & R = keypointsRight[m._j];
    svgStream.drawLine(L.x(), L.y(), R.x()+imageSizeLeft.first, R.y(), svg::svgStyle().stroke("green", std::min(strokeRight,strokeLeft)));
    svgStream.drawCircle(L.x(), L.y(), radiusLeft, svg::svgStyle().stroke("yellow", strokeLeft));
    svgStream.drawCircle(R.x()+imageSizeLeft.first, R.y(), radiusRight, svg::svgStyle().stroke("yellow", strokeRight));
  }
 
  std::ofstream svgFile( outputSVGPath.c_str() );
  svgFile << svgStream.closeSvgFile().str();
  svgFile.close();
}


void saveKeypoints2SVG(const std::string &inputImagePath,
                       const std::pair<size_t,size_t> & imageSize,
                       const std::vector<features::SIOPointFeature> &keypoints,
                       const std::string &outputSVGPath,
                       bool richKeypoint /*= true*/)
{
  svg::svgDrawer svgStream( imageSize.first, imageSize.second);
  svgStream.drawImage(inputImagePath, imageSize.first, imageSize.second);
  // heuristic for the radius of the feature according to the size of the image
  // the larger the distance the larger the minimum radius should be in order to be visible
  // We consider a minimum of 2 pixel and we increment it linearly according to 
  // the image size
  const float radius = getRadiusEstimate(imageSize);
  const float strokeWidth = getStrokeEstimate(imageSize);
  
  for(const features::SIOPointFeature &kpt : keypoints) 
  {
    svgStream.drawCircle(kpt.x(), kpt.y(), (richKeypoint) ? kpt.scale()*radius : radius, svg::svgStyle().stroke("yellow", strokeWidth));
    if(richKeypoint)
    {
      // compute the coordinate of the point on the circle used to draw the orientation line
      const float pointX = kpt.x()+std::cos(kpt.orientation())*kpt.scale()*radius;
      const float pointY = kpt.y()+std::sin(kpt.orientation())*kpt.scale()*radius;
      svgStream.drawLine(kpt.x(), kpt.y(), 
                         pointX, pointY,
                         svg::svgStyle().stroke("yellow", strokeWidth));
    }
  }
 
  std::ofstream svgFile( outputSVGPath );
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
  // heuristic for the radius of the feature according to the size of the image
  // the larger the distance the larger the minimum radius should be in order to be visible
  // We consider a minimum of 2 pixel and we increment it linearly according to 
  // the image size
  const float radius = getRadiusEstimate(imageSize);
  const float strokeWidth = getStrokeEstimate(imageSize);
  
  for(const features::PointFeature &kpt : keypoints) 
  {
    svgStream.drawCircle(kpt.x(), kpt.y(), radius, svg::svgStyle().stroke("yellow", strokeWidth));
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
  // heuristic for the radius of the feature according to the size of the image
  // the larger the distance the larger the minimum radius should be in order to be visible
  // We consider a minimum of 2 pixel and we increment it linearly according to 
  // the image size
  const float radius = getRadiusEstimate(imageSize);
  const float strokeWidth = getStrokeEstimate(imageSize);
  
  if(!inliers)
  {
    for(std::size_t i=0; i < points.cols(); ++i) 
    {
      svgStream.drawCircle(points(0,i), points(1,i), radius, svg::svgStyle().stroke("yellow", strokeWidth));
    }
  }
  else
  {
    for(std::size_t i=0; i < points.cols(); ++i) 
    {
      if(std::find(inliers->begin(), inliers->end(), i) != inliers->end())
        svgStream.drawCircle(points(0,i), points(1,i), radius, svg::svgStyle().stroke("green", strokeWidth));
      else
        svgStream.drawCircle(points(0,i), points(1,i), radius, svg::svgStyle().stroke("yellow", strokeWidth));
    }   
  }
 
  std::ofstream svgFile( outputSVGPath );
  svgFile << svgStream.closeSvgFile().str();
  svgFile.close();
}

bool lineToBorderPoints(const Vec3 &epiLine, const std::size_t imgW, const std::size_t imgH, std::vector<Vec2> &intersectionPts)
{
  intersectionPts.clear();
  intersectionPts.reserve(2);
  // @TODO check special case of epiline coincident with the border lines
  
  // intersect epiline with x=0
  //y = -(a*0+c)/b
  double p = - epiLine(2)/epiLine(1);
  if(p >= 0 && p <= imgH)
    intersectionPts.emplace_back(0, p);
  
  // intersect epiline with x=imgW
  //y = -(a*imgW+c)/b
  p = - (imgW*epiLine(0) + epiLine(2))/epiLine(1);
  if(p >= 0 && p <= imgH)
    intersectionPts.emplace_back(imgW, p);
  
  if(intersectionPts.size()==2)
    return true;
  
  // intersect epiline with y=0
  //x = -(b*0+c)/a
  p = - epiLine(2)/epiLine(0);
  if(p >= 0 && p <= imgW)
    intersectionPts.emplace_back(p, 0);
  
  if(intersectionPts.size()==2)
    return true;
  
  // intersect epiline with y=imgH
  //y = -(a*imgW+c)/b
  p = - (imgH*epiLine(1) + epiLine(2))/epiLine(0);
  if(p >= 0 && p <= imgW)
    intersectionPts.emplace_back(p, imgH);
  
  return (intersectionPts.size()==2);
  
}

void saveEpipolarGeometry2SVG(const std::string &imagePath,
                              const std::pair<size_t, size_t> & imageSize,
                              const std::vector<features::PointFeature> &keypoints,
                              const std::vector<features::PointFeature> &otherKeypoints,
                              const matching::IndMatches &matches,
                              const Mat3 &Fmat,
                              const std::string &outputSVGPath,
                              bool left)
{
  svg::svgDrawer svgStream(imageSize.first, imageSize.second);
  svgStream.drawImage(imagePath, imageSize.first, imageSize.second);
  std::size_t count = 0;
  // heuristic for the radius of the feature according to the size of the image
  // the larger the distance the larger the minimum radius should be in order to be visible
  // We consider a minimum of 2 pixel and we increment it linearly according to 
  // the image size
  const float radius = getRadiusEstimate(imageSize);
  const float strokeWidth = getStrokeEstimate(imageSize);
  for(const matching::IndMatch &m : matches)
  {
    //Get back linked feature, draw a circle and link them by a line
    features::PointFeature p;
    features::PointFeature other;
    if(left)
    {
      p = keypoints[m._i];
      other = otherKeypoints[m._j];
    }
    else
    {
      p = keypoints[m._j];
      other = otherKeypoints[m._i];
    }
    if(count > 7)
      svgStream.drawCircle(p.x(), p.y(), radius, svg::svgStyle().stroke("yellow", strokeWidth));
    else
      svgStream.drawCircle(p.x(), p.y(), radius, svg::svgStyle().stroke("red", strokeWidth).fill("red"));

    Vec3 epiLine;
    if(left)
    {
      epiLine = Fmat.transpose() * Vec3(other.x(), other.y(), 1.0);
    }
    else
    {
      epiLine = Fmat * Vec3(other.x(), other.y(), 1.0);
    }

    //    std::cout << "test 1 o*F*p " << (Fmat*Vec3(p.x(), p.y(), 1.0)).transpose()*Vec3(other.x(), other.y(), 1.0) << std::endl;
    //    std::cout << "test 2 p*F*o " << (Fmat.transpose()*Vec3(p.x(), p.y(), 1.0)).transpose()*Vec3(other.x(), other.y(), 1.0) << std::endl;
    //    std::cout << "epiline\n" << epiLine << " dotprod " << (epiLine.dot(Vec3(p.x(), p.y(), 1.0))) << std::endl;
    std::vector<Vec2> pts;
    if(lineToBorderPoints(epiLine, imageSize.first, imageSize.second, pts))
    {
      //      std::cout << "pt1*epiline " << (epiLine.transpose()*Vec3(pts[0](0), pts[0](1), 1)) << std::endl;
      //      std::cout << "pt1 " << pts[0] << std::endl;
      //      std::cout << "pt2*epiline " << (epiLine.transpose()*Vec3(pts[1](0), pts[1](1), 1)) << std::endl;
      //      std::cout << "pt2 " << pts[1] << std::endl;
      if(count > 7)
        svgStream.drawLine(pts[0](0), pts[0](1), pts[1](0), pts[1](1), svg::svgStyle().stroke("green", strokeWidth));
      else
        svgStream.drawLine(pts[0](0), pts[0](1), pts[1](0), pts[1](1), svg::svgStyle().stroke("red", strokeWidth));
    }
    else
    {
      std::cerr << "********** pts size: " << pts.size() << " epiline " << epiLine << " out of image" << std::endl;
      if(pts.size() > 0)
      {
        svgStream.drawLine(pts[0](0), pts[0](1), p.x(), p.y(), svg::svgStyle().stroke("red", 10 * strokeWidth));
        std::cerr << "********** pts: " << pts[0].transpose() << std::endl;
      }
    }
    ++count;
    //    if(count > 7) break;

  }

  //draw the epipole
  Mat epipole;
  if(left)
    epipole = Fmat.fullPivLu().kernel();
  else
    epipole = Fmat.transpose().fullPivLu().kernel();

  if(epipole.cols() > 1)
  {
    std::cerr << "F has kernel of size " << epipole.cols() << std::endl << epipole << std::endl;
  }
  else
  {
    // normalize coordinates
    Vec3 point = epipole.col(0);
    std::cout << "epipole:\n" << point << std::endl;
    //@todo check 0
    point /= point(2);
    std::cout << "epipole normalized:\n" << point << std::endl;
    // check if the point is inside the image
    if(!((point(0) > 0 && point(0) < imageSize.first) &&
            (point(1) > 0 && point(1) < imageSize.second)))
    {
      std::cout << "epipole outside the image:\n" << point << std::endl;
      // point outside the image, clamp it to the borders
      if(point(0) < 0) point(0) = 0;
      if(point(0) > imageSize.first) point(0) = imageSize.first;
      if(point(1) < 0) point(1) = 0;
      if(point(1) > imageSize.second) point(0) = imageSize.second;
      std::cout << "clamped epipole:\n" << point << std::endl;
    }
    svgStream.drawCircle(point(0), point(1), 3 * radius, svg::svgStyle().stroke("red", strokeWidth).fill("red"));
  }

  std::ofstream svgFile(outputSVGPath.c_str());
  svgFile << svgStream.closeSvgFile().str();
  svgFile.close();
}

void saveMatchesAsMotion(const std::string &imagePath,
                         const std::pair<size_t, size_t> & imageSize,
                         const std::vector<features::SIOPointFeature> &keypoints,
                         const std::vector<features::SIOPointFeature> &otherKeypoints,
                         const matching::IndMatches &matches,
                         const std::string &outputSVGPath,
                         bool left,
                         bool richKeypoint)
{
  svg::svgDrawer svgStream(imageSize.first, imageSize.second);
  svgStream.drawImage(imagePath, imageSize.first, imageSize.second);

  const float radius = getRadiusEstimate(imageSize);
  const float strokeWidth = getStrokeEstimate(imageSize);
  for(size_t i = 0; i < matches.size(); ++i)
  {
    //Get back linked feature, draw a circle and link them by a line
    const auto L = keypoints[matches[i]._i];
    const auto R = otherKeypoints[matches[i]._j];
    if(left)
    {
      svgStream.drawLine(L.x(), L.y(), R.x(), R.y(), svg::svgStyle().stroke("green", strokeWidth));
      svgStream.drawCircle(L.x(), L.y(), (richKeypoint) ? L.scale()*radius : radius, svg::svgStyle().stroke("yellow", 2.0));
      svgStream.drawCircle(R.x(), R.y(), (richKeypoint) ? R.scale()*radius : radius, svg::svgStyle().stroke("red", strokeWidth));
    }
    else
    {
      svgStream.drawLine(L.x(), L.y(), R.x(), R.y(), svg::svgStyle().stroke("green", strokeWidth));
      svgStream.drawCircle(R.x(), R.y(), (richKeypoint) ? R.scale()*radius : radius, svg::svgStyle().stroke("yellow", 2.0));
      svgStream.drawCircle(L.x(), L.y(), (richKeypoint) ? L.scale()*radius : radius, svg::svgStyle().stroke("red", strokeWidth));

    }
  }
  std::ofstream svgFile(outputSVGPath.c_str());
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
  
  const float radius = getRadiusEstimate(imageSize);
  const float strokeWidth = getStrokeEstimate(imageSize);
    
  const auto &feat = cctags.Features();
  const auto &desc = cctags.Descriptors();
  
  for(std::size_t i = 0; i < desc.size(); ++i) 
  {
    const IndexT cctagId = features::getCCTagId(desc[i]);
    if ( cctagId == UndefinedIndexT)
    {
      continue;
    }
    const auto &kpt = feat[i];
    svgStream.drawCircle(kpt.x(), kpt.y(), radius, svg::svgStyle().stroke("yellow", strokeWidth));
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
  const auto &descLeft = cctagLeft.Descriptors();
  const auto &descRight = cctagRight.Descriptors();
  
  //just to be sure...
  assert(keypointsLeft.size() == descLeft.size());
  assert(keypointsRight.size() == descRight.size());
  
  const float radiusLeft = getRadiusEstimate(imageSizeLeft);
  const float radiusRight = getRadiusEstimate(imageSizeRight);
  const float strokeLeft = getStrokeEstimate(imageSizeLeft);
  const float strokeRight = getStrokeEstimate(imageSizeRight);
  
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
    
    svgStream.drawLine(L.x(), L.y(), R.x() + imageSizeLeft.first, R.y(), svg::svgStyle().stroke("green", std::min(imageSizeLeft,imageSizeRight)));
    svgStream.drawCircle(L.x(), L.y(), radiusLeft, svg::svgStyle().stroke("yellow", imageSizeLeft));
    svgStream.drawCircle(R.x() + imageSizeLeft.first, R.y(), radiusRight, svg::svgStyle().stroke("yellow", imageSizeRight));
    
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
      svgStream.drawCircle(L.x(), L.y(), radiusLeft, svg::svgStyle().stroke("red", strokeLeft));
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
      svgStream.drawCircle(R.x() + imageSizeLeft.first, R.y(), radiusRight, svg::svgStyle().stroke("red", strokeRight));
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

