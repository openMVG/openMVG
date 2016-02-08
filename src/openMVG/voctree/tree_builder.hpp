#ifndef OPENMVG_VOCABULARY_TREE_TREE_BUILDER_HPP
#define OPENMVG_VOCABULARY_TREE_TREE_BUILDER_HPP

#include "mutable_tree.hpp"
#include "simple_kmeans.hpp"
#include <deque>
//#include <cstdio> //DEBUG

namespace openMVG {
namespace voctree {

/**
 * @brief Class for building a new vocabulary by hierarchically clustering
 * a set of training features.
 */
template<class Feature, template<typename, typename> class DistanceT = L2,
class FeatureAllocator = typename DefaultAllocator<Feature>::type>
class TreeBuilder
{
public:
  typedef MutableVocabularyTree<Feature, DistanceT, FeatureAllocator> Tree;
  typedef DistanceT<Feature, Feature> Distance;
  typedef SimpleKmeans<Feature, Distance, FeatureAllocator> Kmeans;
  typedef std::vector<Feature, FeatureAllocator> FeatureVector;

  /**
   * @brief Constructor
   *
   * @param zero Object representing zero in the feature space
   * @param d    Functor for calculating squared distance
   */
  TreeBuilder(const Feature& zero = Feature(), Distance d = Distance(), unsigned char verbose = 0);

  /**
   * @brief Build a new vocabulary tree.
   *
   * The number of words in the resulting vocabulary is at most k ^ levels.
   *
   * @param training_features The set of training features to cluster.
   * @param k                 The branching factor, or max children of any node.
   * @param levels            The number of levels in the tree.
   */
  void build(const FeatureVector& training_features, uint32_t k, uint32_t levels);

  /// Get the built vocabulary tree.

  const Tree& tree() const
  {
    return tree_;
  }

  /// Get the k-means clusterer.

  Kmeans& kmeans()
  {
    return kmeans_;
  }
  /// Get the k-means clusterer.

  const Kmeans& kmeans() const
  {
    return kmeans_;
  }

  void setVerbose(unsigned char level)
  {
    verbose_ = level;
    kmeans_.setVerbose(level);
  }

  unsigned char getVerbose() const
  {
    return verbose_;
  }

protected:
  Tree tree_;
  Kmeans kmeans_;
  Feature zero_;
private:
  unsigned char verbose_;
};

template<class Feature, template<typename, typename> class DistanceT, class FeatureAllocator>
TreeBuilder<Feature, DistanceT, FeatureAllocator>::TreeBuilder(const Feature& zero, Distance d, unsigned char verbose)
: kmeans_(zero, d, verbose),
zero_(zero),
verbose_(verbose)
{
}

template<class Feature, template<typename, typename> class DistanceT, class FeatureAllocator>
void TreeBuilder<Feature, DistanceT, FeatureAllocator>::build(const FeatureVector& training_features,
                                                             uint32_t k, uint32_t levels)
{
  // Initial setup and memory allocation for the tree
  tree_.clear();
  tree_.setSize(levels, k);
  tree_.centers().reserve(tree_.nodes());
  tree_.validCenters().reserve(tree_.nodes());

  // We keep a queue of disjoint feature subsets to cluster.
  // Feature* is used to avoid copying features.
  std::deque< std::vector<Feature*> > subset_queue(1);

  // At first the queue contains one "subset" containing all the features.
  std::vector<Feature*> &feature_ptrs = subset_queue.front();
  feature_ptrs.reserve(training_features.size());
  BOOST_FOREACH(const Feature& f, training_features)
  feature_ptrs.push_back(const_cast<Feature*> (&f));

  FeatureVector centers; // always size k
  for(uint32_t level = 0; level < levels; ++level)
  {
    if(verbose_) printf("# Level %u\n", level);
    std::vector<unsigned int> membership;

    for(size_t i = 0, ie = subset_queue.size(); i < ie; ++i)
    {
      std::vector<Feature*> &subset = subset_queue.front();
      if(verbose_ > 1) printf("#\tClustering subset %lu/%lu of size %lu\n", i + 1, ie, subset.size());

      // If the subset already has k or fewer elements, just use those as the centers.
      if(subset.size() <= k)
      {
        if(verbose_ > 2) printf("#\tno need to cluster %lu elements\n", subset.size());
        for(size_t j = 0; j < subset.size(); ++j)
        {
          tree_.centers().push_back(*subset[j]);
          tree_.validCenters().push_back(1);
        }
        // Mark non-existent centers as invalid.
        tree_.centers().insert(tree_.centers().end(), k - subset.size(), zero_);
        tree_.validCenters().insert(tree_.validCenters().end(), k - subset.size(), 0);

        // Push k empty subsets into the queue so all children get marked invalid.
        subset_queue.pop_front();
        subset_queue.insert(subset_queue.end(), k, std::vector<Feature*>());
      }
      else
      {
        // Cluster the current subset into k centers.
        if(verbose_ > 2) printf("#\tclustering the current subset of %lu elements into %d centers\n", subset.size(), k);
        kmeans_.clusterPointers(subset, k, centers, membership);
        // Add the centers and mark them as valid.
        tree_.centers().insert(tree_.centers().end(), centers.begin(), centers.end());
        tree_.validCenters().insert(tree_.validCenters().end(), k, 1);
        //printf("LOL\n");
        // Partition the current subset into k new subsets based on the cluster assignments.
        std::vector< std::vector<Feature*> > new_subsets(k);
        for(size_t j = 0; j < subset.size(); ++j)
        {
          new_subsets[ membership[j] ].push_back(subset[j]);
        }
        //printf("LOL\n");
        // Update the queue
        subset_queue.pop_front();
        subset_queue.insert(subset_queue.end(), new_subsets.begin(), new_subsets.end());
      }
    }
    if(verbose_) printf("# centers so far = %lu\n", tree_.centers().size());
  }
}

}
}

#endif
