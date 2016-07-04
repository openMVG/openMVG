#pragma once

#include <openMVG/matching/indMatch.hpp>
#include <openMVG/features/regions.hpp>
#include <openMVG/features/regions_factory.hpp>
#include <openMVG/image/image.hpp>
#include <openMVG/features/features.hpp>
#include <openMVG/matching_image_collection/Matcher_Regions_AllInMemory.hpp>

#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"

namespace openMVG {
namespace features {

/**
* @brief Compute the n best matches.
* @param[in] inputMatches Set of indices for (putative) matches.
* @param[in] regionsI Pointer to the regions of the left image.
* @param[in] regionsJ Pointer to the regions of the right image.
* @param[out] outputMatches Subset of inputMatches containing the best n matches, sorted.
*/
void sortMatches(
	const openMVG::matching::IndMatches& inputMatches,
	const openMVG::features::Feat_Regions<openMVG::features::SIOPointFeature> *regionsI,
	const openMVG::features::Feat_Regions<openMVG::features::SIOPointFeature> *regionsJ,
	openMVG::matching::IndMatches& outputMatches);


/**
* @brief Compare method used in the match sorting.
* @param[in] firstElem The first element to be compared.
* @param[in] secondElem The second element to be compared.
* @return True if firstElem is less than secondElem.
*/
bool matchCompare(const std::pair<float, size_t>& firstElem, const std::pair<float, size_t>& secondElem);

/**
* @brief Extracts by copy the first (and best) uNumMatchesToKeep.
* @param[out] outputMatches Set of image pairs and their respective sets of matches thresholded to the first uNumMatchesToKeep.
* @param[in] uNumMatchesToKeep The N best matches to keep.
*/
void thresholdMatches(openMVG::matching::IndMatches& outputMatches, const std::size_t uNumMatchesToKeep);

/**
 * @brief Perform the gris filtering on the matches
 * @param[in] lRegions The regions of the first picture
 * @param[in] rRegions The regions of the second picture
 * @param[in] indexImagePair The Pair of matched images
 * @param[in] sfm_data The sfm data file
 * @param[out] outMatches The remaining matches
 */
void matchesGridFiltering(const openMVG::features::Feat_Regions<openMVG::features::SIOPointFeature>* lRegions, 
        const openMVG::features::Feat_Regions<openMVG::features::SIOPointFeature>* rRegions, 
        const openMVG::Pair& indexImagePair,
        const openMVG::sfm::SfM_Data sfm_data, 
        openMVG::matching::IndMatches& outMatches);

}
}
