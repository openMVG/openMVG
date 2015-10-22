/* 
 * File:   svgVisualization.hpp
 * Author: sgaspari
 *
 * Created on October 19, 2015, 9:46 AM
 */

#pragma once

#include <openMVG/features/image_describer.hpp>
#include <openMVG/matching/indMatch.hpp>

#include <vector>

namespace openMVG {
namespace localization {

/**
 * @brief
 * 
 * @param[in] imagePathLeft
 * @param[in] imageSizeLeft
 * @param[in] keypointsLeft
 * @param[in] imagePathRight
 * @param[in] imageSizeRight
 * @param[in] keypointsRight
 * @param[in] matches
 * @param[in] outputSVGPath
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
 * @param[in] inputImagePath
 * @param[in] imageSize
 * @param[in] keypoints
 * @param[in] outputSVGPath
 */
void saveFeatures2SVG(const std::string &inputImagePath,
                      const std::pair<size_t,size_t> & imageSize,
                      const std::vector<features::PointFeature> &keypoints,
                      const std::string &outputSVGPath);

} // namespace localization
} // namespace openMVG

