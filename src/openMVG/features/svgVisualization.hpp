/* 
 * File:   svgVisualization.hpp
 * Author: sgaspari
 *
 * Created on October 19, 2015, 9:46 AM
 */

#pragma once

#include "image_describer.hpp"
#if HAVE_CCTAG
#include "regions_factory.hpp"
#endif
#include <openMVG/matching/indMatch.hpp>

#include <vector>

namespace openMVG {
namespace features {

/**
 * @brief It saves a svg file containing two images (as linked images) and their
 * feature matches: the two images are showed side by side and each feature in each
 * image (depicted as a circle) is connected to the corresponding feature on the 
 * other image through a line.
 * 
 * @param[in] imagePathLeft The full path to the left iamge. The image is only 
 * saved as a link, no image data is stored in the svg.
 * @param[in] imageSizeLeft The size of the image <width,height>.
 * @param[in] keypointsLeft The keypoints of the left image.
 * @param[in] imagePathRight The full path to the left iamge. The image is only 
 * saved as a link, no image data is stored in the svg.
 * @param[in] imageSizeRight The size of the image <width,height>.
 * @param[in] keypointsRight The keypoints of the right image.
 * @param[in] matches The vector containing the indices of matching cctags.
 * @param[in] outputSVGPath The name of the svg file to generate.
 */
void saveMatches2SVG(const std::string &imagePathLeft,
                     const std::pair<size_t,size_t> & imageSizeLeft,
                     const std::vector<features::PointFeature> &keypointsLeft,
                     const std::string &imagePathRight,
                     const std::pair<size_t,size_t> & imageSizeRight,
                     const std::vector<features::PointFeature> &keypointsRight,
                     const matching::IndMatches &matches,
                     const std::string &outputSVGPath);

/**
 * @brief It saves a svg file containing an image (as linked image) and its detected
 * features.
 * 
 * @param[in] inputImagePath The full path to the image file. The image is only 
 * saved as a link, no image data is stored in the svg.
 * @param[in] imageSize The size of the image <width,height>.
 * @param[in] keypoints The keypoints of the right image.
 * @param[in] outputSVGPath The name of the svg file to generate.
 * @param[in] richKeypoint Draw rich keypoints with a circle proportional to the 
 * octave in which the point has been detected. 
 */
void saveKeypoints2SVG(const std::string &inputImagePath,
                       const std::pair<size_t,size_t> & imageSize,
                       const std::vector<features::SIOPointFeature> &keypoints,
                       const std::string &outputSVGPath,
                       bool richKeypoint = true);
/**
 * @brief It saves a svg file containing an image (as linked image) and its detected
 * features.
 * 
 * @param[in] inputImagePath The full path to the image file. The image is only 
 * saved as a link, no image data is stored in the svg.
 * @param[in] imageSize The size of the image <width,height>.
 * @param[in] keypoints The points of the right image.
 * @param[in] outputSVGPath The name of the svg file to generate.
 **/
void saveFeatures2SVG(const std::string &inputImagePath,
                      const std::pair<size_t,size_t> & imageSize,
                      const std::vector<features::PointFeature> &keypoints,
                      const std::string &outputSVGPath);

/**
 * @brief It saves a svg file containing an image (as linked image) and its detected
 * features.
 * 
 * @param[in] inputImagePath The full path to the image file. The image is only 
 * saved as a link, no image data is stored in the svg.
 * @param[in] imageSize The size of the image <width,height>.
 * @param[in] points A vector containing the points to draw.
 * @param[in] outputSVGPath The name of the svg file to generate.
 * @param[in] inliers [optional] The indices of the features to draw.
 */
void saveFeatures2SVG(const std::string &inputImagePath,
                      const std::pair<size_t,size_t> & imageSize,
                      const Mat &points,
                      const std::string &outputSVGPath,
                      const std::vector<size_t> *inliers = nullptr);

/**
 * @brief It saves an image as a svg file (linked image) with an overlay of the 
 * epipolar geometry. For each feature of the image it displays the features itself
 * as a circle and the corresponding epipolar line. The epipolar line is computed
 * by using the provided F and match of a second image.
 * 
 * @param[in] imagePath The full path to the image file to display.
 * @param[in] imageSize The size of the image <width, height>.
 * @param[in] keypoints The list of keypoints associated to the image.
 * @param[in] otherKeypoints The list of keypoints associated to the other image, 
 * they are used to draw the epipolar lines.
 * @param[in] matches The correspondences between the keypoints of the current 
 * image and the other image.
 * @param[in] Fmat The fundamental matrix between the current image and the other 
 * image, in particular it is assumed the following right' * F * left
 * @param[in] outputSVGPath The filename for the svg file to generate
 * @param[in] left If true it will consider the current image as the left image
 */
void saveEpipolarGeometry2SVG(const std::string &imagePath,
                              const std::pair<size_t, size_t> & imageSize,
                              const std::vector<features::PointFeature> &keypoints,
                              const std::vector<features::PointFeature> &otherKeypoints,
                              const matching::IndMatches &matches,
                              const Mat3 &Fmat,
                              const std::string &outputSVGPath,
                              bool left);

/**
 * @brief It saves a svg file containing an image (as linked image) and its
 * feature matches: the matches are depicted as trails, ie both feature and its 
 * corresponding match are drawn on the same image and connected by a line. The 
 * keypoint of the current image is drawn in yellow, the other in red.
 * 
 * @param[in] imagePath The full path to the image file to display.
 * @param[in] imageSize The size of the image <width, height>.
 * @param[in] keypoints The list of keypoints associated to the image.
 * @param[in] otherKeypoints The list of keypoints associated to the other image, they are used to draw the epipolar lines.
 * @param[in] matches The correspondences between the keypoints of the current image and the other image.
 * @param[in] outputSVGPath The filename for the svg file to generate.
 * @param[in] left If true it will consider the current image as the left image.
 * @param[in] richKeypoint Draw rich keypoints with a circle proportional to the 
 * octave in which the point has been detected.
 */
void saveMatchesAsMotion(const std::string &imagePath,
                         const std::pair<size_t, size_t> & imageSize,
                         const std::vector<features::SIOPointFeature> &keypoints,
                         const std::vector<features::SIOPointFeature> &otherKeypoints,
                         const matching::IndMatches &matches,
                         const std::string &outputSVGPath,
                         bool left,
                         bool richKeypoint = true);

/**
 * @brief Given a 2d line and an image size it returns the intersection points 
 * of the line with the image borders. It returns false if there is no intersection
 * points
 * @param[in] epiLine The 2D line parameters, [a, b, c] for ax+by+c=0.
 * @param[in] imgW The image width.
 * @param[in] imgH The image height
 * @param[out] intersectionPts A list of intersection points
 * @return False if there is no intersection points.
 */
bool lineToBorderPoints(const Vec3 &epiLine, 
                        const std::size_t imgW, 
                        const std::size_t imgH, 
                        std::vector<Vec2> &intersectionPts);

#if HAVE_CCTAG

/**
 * @brief It generates a svg file containing the image and its extracted cctags.
 * The center of the cctag is marked with a small circle and the id of the cctag
 * is rendered as text close to the center.
 * 
 * @param[in] inputImagePath The full path to the image file. The image is only 
 * saved as a link, no image data is stored in the svg.
 * @param[in] imageSize The size of the image <width,height>.
 * @param[in] cctags The CCtag regions (keypoints+descriptors).
 * @param[in] outputSVGPath The name of the svg file to generate.
 */
void saveCCTag2SVG(const std::string &inputImagePath,
                      const std::pair<size_t,size_t> & imageSize,
                      const features::CCTAG_Regions &cctags,
                      const std::string &outputSVGPath);

/**
 * @brief It generates a svg file containing the two images, its extracted cctags, 
 * and their correspondences.
 * Each correspondence is drawn as a line connecting the two centers drawn as a small circle.
 * The ids of the cctags are rendered as text close to their center.
 * If \p showNotMatched is enable also the non matching cctags are drawn.
 * 
 * @param[in] imagePathLeft The full path to the left iamge. The image is only 
 * saved as a link, no image data is stored in the svg.
 * @param[in] imageSizeLeft The size of the image <width,height>.
 * @param[in] cctagLeft The CCtag regions (keypoints+descriptors).
 * @param[in] imagePathRight The full path to the left iamge. The image is only 
 * saved as a link, no image data is stored in the svg.
 * @param[in] imageSizeRight The size of the image <width,height>.
 * @param[in] cctagRight The CCtag regions (keypoints+descriptors).
 * @param[in] matches The vector containing the indices of matching cctags.
 * @param[in] outputSVGPath The name of the svg file to generate.
 * @param[in] showNotMatched If enabled, even the non matched cctags are drawn.
 */
void saveCCTagMatches2SVG(const std::string &imagePathLeft,
                     const std::pair<size_t,size_t> & imageSizeLeft,
                     const features::CCTAG_Regions &cctagLeft,
                     const std::string &imagePathRight,
                     const std::pair<size_t,size_t> & imageSizeRight,
                     const features::CCTAG_Regions &cctagRight,
                     const matching::IndMatches &matches,
                     const std::string &outputSVGPath, 
                     bool showNotMatched);
#endif

} // namespace features
} // namespace openMVG

