#include "selection.hpp"

namespace openMVG {
namespace features {

const size_t gridSize = 3;
  
/**
* @brief Sort the matches.
* @param[in] inputMatches Set of indices for (putative) matches.
* @param[in] regionsI Pointer to the regions of the left image.
* @param[in] regionsJ Pointer to the regions of the right image.
* @param[out] outputMatches Subset of inputMatches containing the best n matches, sorted.
*/
void sortMatches(
	const openMVG::matching::IndMatches& inputMatches,
	const openMVG::features::Feat_Regions<openMVG::features::SIOPointFeature> *regionsI,
	const openMVG::features::Feat_Regions<openMVG::features::SIOPointFeature> *regionsJ,
	openMVG::matching::IndMatches& outputMatches)
{
	if (!regionsI || !regionsJ)
  {
    throw std::runtime_error("Can't sort matches: One or both regions provided were NULL.");
	}
	const std::vector<openMVG::features::SIOPointFeature>& vecFeatureI = regionsI->Features();
	const std::vector<openMVG::features::SIOPointFeature>& vecFeatureJ = regionsJ->Features();

	//outputMatches will contain the sorted matches if inputMatches.
	outputMatches.reserve(inputMatches.size());

	//This vector is just a temporary container to link the index of the matches in the original vector inputMatches.
	//It will be used to retrieve the correct matches after the sort.
	std::vector<std::pair<float, size_t>> vecFeatureScale;

	for (size_t i = 0; i < inputMatches.size(); i++)
  {
		float scale1 = vecFeatureI[inputMatches[i]._i].scale();
		float scale2 = vecFeatureJ[inputMatches[i]._j].scale();
		vecFeatureScale.emplace_back((scale1 + scale2) / 2.0, i);
	}

	std::sort(vecFeatureScale.begin(), vecFeatureScale.end(), matchCompare);

	//The sorted match vector is filled according to the result of the sorting above.
	for (size_t i = 0; i < vecFeatureScale.size(); i++)
  {
		outputMatches.push_back(inputMatches[vecFeatureScale[i].second]);
	}
}

/**
* @brief Compare method used in the match sorting.
* @param[in] firstElem The first element to be compared.
* @param[in] secondElem The second element to be compared.
* @return True if firstElem is less than secondElem.
*/
bool matchCompare(const std::pair<float, size_t>& firstElem, const std::pair<float, size_t>& secondElem)
{
	return firstElem.first > secondElem.first;
}

/**
* @brief Extracts by copy the first (and best) uNumMatchesToKeep.
* @param[out] outputMatches Set of image pairs and their respective sets of matches thresholded to the first uNumMatchesToKeep.
* @param[in] uNumMatchesToKeep The N best matches to keep.
*/
void thresholdMatches(openMVG::matching::IndMatches& outputMatches, const std::size_t uNumMatchesToKeep)
{
	if (outputMatches.size() > uNumMatchesToKeep) {
		outputMatches.resize(uNumMatchesToKeep);
	}
}

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
        openMVG::matching::IndMatches& outMatches)
{
  const std::size_t lWidth = sfm_data.GetViews().at(indexImagePair.first)->ui_width;
  const std::size_t lHeight = sfm_data.GetViews().at(indexImagePair.first)->ui_height;
  const std::size_t rWidth = sfm_data.GetViews().at(indexImagePair.second)->ui_width;
  const std::size_t rHeight = sfm_data.GetViews().at(indexImagePair.second)->ui_height;
  
  const size_t leftCellHeight = std::ceil(lHeight / (float)gridSize);
  const size_t leftCellWidth = std::ceil(lWidth / (float)gridSize);
  const size_t rightCellHeight = std::ceil(rHeight / (float)gridSize);
  const size_t rightCellWidth = std::ceil(rWidth / (float)gridSize);

  std::vector< openMVG::matching::IndMatches > completeGrid(gridSize*gridSize*2);
  // Reserve all cells
  for(openMVG::matching::IndMatches& cell: completeGrid)
  {
    cell.reserve(outMatches.size()/completeGrid.size());
  }
  // Split matches in grid cells
  for(const auto& match: outMatches)
  {
    const openMVG::features::SIOPointFeature& leftPoint = lRegions->Features()[match._i];
    const openMVG::features::SIOPointFeature& rightPoint = rRegions->Features()[match._j];
    
    const size_t leftGridIndex = std::floor(leftPoint.x() / (float)leftCellWidth) + std::floor(leftPoint.y() / (float)leftCellHeight) * gridSize;
    const size_t rightGridIndex = std::floor(rightPoint.x() / (float)rightCellWidth) + std::floor(rightPoint.y() / (float)rightCellHeight) * gridSize;
    
    openMVG::matching::IndMatches& currentCaseL = completeGrid[leftGridIndex];
    openMVG::matching::IndMatches& currentCaseR = completeGrid[rightGridIndex + gridSize*gridSize];
    
    if(currentCaseL.size() > currentCaseR.size())
    {
      currentCaseR.push_back(match);
    }
    else
    {
      currentCaseL.push_back(match);
    }
  }
  
  // max Size of the cells:
  int maxSize = 0;
  for (const auto& cell: completeGrid)
  {
    if(cell.size() > maxSize)
    {
      maxSize = cell.size();
    }
  }
   
  openMVG::matching::IndMatches finalMatches;
  finalMatches.reserve(outMatches.size());
  
  // Combine all cells into a global ordered vector
  for (int cmpt = 0; cmpt < maxSize; ++cmpt)
  {
    for(const auto& cell: completeGrid)
    {
      if(cell.size() > cmpt)
      {
        finalMatches.push_back(cell[cmpt]);
      }
    }
  }
  
  outMatches.swap(finalMatches);
}

}
}
