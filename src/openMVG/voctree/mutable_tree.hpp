#ifndef OPENMVG_VOCABULARY_TREE_MUTABLE_TREE_HPP
#define OPENMVG_VOCABULARY_TREE_MUTABLE_TREE_HPP

#include "vocabulary_tree.hpp"

namespace openMVG {
namespace voctree {

/**
 * @brief Vocabulary tree that exposes the hierarchical clustering centers. Mainly
 * intended for building a new tree.
 *
 * When loading and using an existing vocabulary tree, use VocabularyTree instead.
 */
template<class Feature, class Distance = L2<Feature>,
class FeatureAllocator = typename DefaultAllocator<Feature>::type>
class MutableVocabularyTree : public VocabularyTree<Feature, Distance, FeatureAllocator>
{
  typedef VocabularyTree<Feature, Distance, FeatureAllocator> BaseClass;

public:

  MutableVocabularyTree(Distance d = Distance())
  : BaseClass(d)
  {
  }

  void setSize(uint32_t levels, uint32_t splits)
  {
    this->levels_ = levels;
    this->k_ = splits;
    this->setNodeCounts();
  }

  uint32_t nodes() const
  {
    return this->word_start_ + this->num_words_;
  }

  std::vector<Feature, FeatureAllocator>& centers()
  {
    return this->centers_;
  }

  const std::vector<Feature, FeatureAllocator>& centers() const
  {
    return this->centers_;
  }

  std::vector<uint8_t>& validCenters()
  {
    return this->valid_centers_;
  }

  const std::vector<uint8_t>& validCenters() const
  {
    return this->valid_centers_;
  }
};

}
}

#endif //OPENMVG_VOCABULARY_TREE_MUTABLE_TREE_HPP
