#pragma once

#include "openMVG/matching/indMatch.hpp"

using namespace openMVG::matching;

#include <map>
#include <string>
#include <vector>

/// The structure used to store corresponding point index per images pairs
typedef std::map< std::pair<size_t, size_t>, std::vector<IndMatch> > IndexedMatchPerPair;

/// Implementation of an Image Collection Matcher
/// Compute putative matches between a collection of pictures
class ImageCollectionMatcher
{
  public:
  ImageCollectionMatcher() {};

  virtual ~ImageCollectionMatcher() {};

  /// Build point indexes correspondences lists between images ids
  virtual void Match(
    const std::vector<std::string> & vec_fileNames, // input filenames,
    IndexedMatchPerPair & map_putatives_matches // the parwise photometric corresponding points
    )const = 0;
};

