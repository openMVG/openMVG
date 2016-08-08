#pragma once

#include <vector>
#include <map>
#include <opencv2/opencv.hpp>

namespace openMVG {
namespace bestImages {

/**
 * @brief This function computes cell indexes per image.
 *
 * @param[in] imagePoints
 * @param[out] cellIndexesPerImage
 * @param[in] imageSize
 * @param[in] calibGridSize
 */
void precomputeCellIndexes(const std::vector<std::vector<cv::Point2f> >& imagePoints,
                           std::vector<std::vector<std::size_t> >& cellIndexesPerImage,
                           const cv::Size& imageSize, const std::size_t calibGridSize);

/**
 * @brief This function counts the number of points in each cell of the grid.
 *
 * @param[in] imagesIndexes
 * @param[in] cellIndexesPerImage
 * @param[out] cellsWeight
 * @param[in] calibGridSize
 */
void computeCellsWeight(const std::vector<std::size_t>& imagesIndexes,
                        const std::vector<std::vector<std::size_t> >& cellIndexesPerImage,
                        std::map<std::size_t, std::size_t>& cellsWeight, const std::size_t calibGridSize);

/**
 * @brief This function computes the score of each image.
 *
 * @param[out] imageScores
 * @param[in] remainingImagesIndexes
 * @param[in] cellIndexesPerImage
 * @param[in] cellsWeight
 */
void computeImageScores(std::vector<std::pair<float, std::size_t> >& imageScores,
                        const std::vector<std::size_t>& remainingImagesIndexes,
                        const std::vector<std::vector<std::size_t> >& cellIndexesPerImage,
                        const std::map<std::size_t, std::size_t>& cellsWeight);

/**
 * @brief This function selects the best images based on repartition in images of the calibration landmarks.
 *
 * @param[in] imagePoints
 * @param[in] imageSize
 * @param[in,out] remainingImagesIndexes
 * @param[in] maxCalibFrames
 * @param[in] validFrames
 * @param[out] calibImageScore
 * @param[out] calibInputFrames
 * @param[out] calibImagePoints
 * @param[in] calibGridSize
 */
void selectBestImages(const std::vector<std::vector<cv::Point2f> >& imagePoints,
                      const cv::Size& imageSize,
                      std::vector<std::size_t>& remainingImagesIndexes,
                      const std::size_t& maxCalibFrames,
                      const std::vector<std::size_t>& validFrames,
                      std::vector<float>& calibImageScore,
                      std::vector<std::size_t>& calibInputFrames,
                      std::vector<std::vector<cv::Point2f> >& calibImagePoints,
                      const std::size_t calibGridSize);

}//namespace bestImages
}//namespace openMVG

