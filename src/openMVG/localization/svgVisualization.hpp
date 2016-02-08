/* 
 * File:   svgVisualization.hpp
 * Author: sgaspari
 *
 * Created on October 19, 2015, 9:46 AM
 */

#pragma once

#include <openMVG/features/image_describer.hpp>
#if HAVE_CCTAG
#include <openMVG/features/regions_factory.hpp>
#endif
#include <openMVG/matching/indMatch.hpp>

#include <vector>

namespace openMVG {
namespace localization {

/**
 * @brief
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
 * @brief
 * 
 * @param[in] inputImagePath The full path to the image file. The image is only 
 * saved as a link, no image data is stored in the svg.
 * @param[in] imageSize The size of the image <width,height>.
 * @param[in] keypoints The keypoints of the right image.
 * @param[in] outputSVGPath The name of the svg file to generate.
 */
void saveFeatures2SVG(const std::string &inputImagePath,
                      const std::pair<size_t,size_t> & imageSize,
                      const std::vector<features::PointFeature> &keypoints,
                      const std::string &outputSVGPath);

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
                      const std::vector<size_t> *inliers = nullptr);

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

} // namespace localization
} // namespace openMVG

