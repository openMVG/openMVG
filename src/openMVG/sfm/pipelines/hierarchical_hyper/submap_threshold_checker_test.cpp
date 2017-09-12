#include "openMVG/sfm/pipelines/hierarchical_hyper/submap_threshold_checker.hpp"

#include "openMVG/sfm/pipelines/pipelines_test.hpp"

#include "testing/testing.h"

HsfmSubmap generateBasicSubmap(int n_views, int n_points)
{
  HsfmSubmap smap;
  const nViewDatasetConfigurator config;
  const NViewDataSet dataset = NRealisticCamerasRing(n_views, n_points, config);
  smap.sfm_data = getInputScene(dataset, config, openMVG::cameras::PINHOLE_CAMERA_RADIAL3);
  const openMVG::sfm::Landmarks & landmarks = smap.sfm_data.structure;
  std::transform(landmarks.begin(), landmarks.end(),
                 std::inserter(smap.track_ids, smap.track_ids.begin()),
                 [](const std::pair<openMVG::IndexT, openMVG::sfm::Landmark> & landmark){return landmark.first;});
  // take only the first half of the points for the separator
  auto it_half = smap.track_ids.begin();
  std::advance(it_half, n_points / 2);
  smap.separator.insert(smap.track_ids.begin(), it_half);

  return smap;
}

TEST(tracksThresholdPredicate, submaps_returnExpectedPredicate)
{
  const int tracks_threshold = 100;

  openMVG::sfm::SubmapTracksThresholdChecker is_partitionable(tracks_threshold);

  // create submaps that don't meet the requirements of the predicate
  const HsfmSubmap smap_not_enough_views = generateBasicSubmap(1, tracks_threshold * 2);
  const HsfmSubmap smap_not_enough_tracks = generateBasicSubmap(10, tracks_threshold - 1);
  HsfmSubmap smap_parent_already = generateBasicSubmap(10, tracks_threshold * 2);
  smap_parent_already.is_parent = true;
  smap_parent_already.children_submaps = {1,2};

  EXPECT_FALSE(is_partitionable(smap_not_enough_views));
  EXPECT_FALSE(is_partitionable(smap_not_enough_tracks));
  EXPECT_FALSE(is_partitionable(smap_parent_already));

  const HsfmSubmap smap_good = generateBasicSubmap(10, tracks_threshold * 2);
  EXPECT_TRUE(is_partitionable(smap_good));
}

TEST(viewsThresholdPredicate, submaps_returnExpectedPredicate)
{
  const int views_threshold = 20;

  openMVG::sfm::SubmapViewThresholdChecker is_partitionable(views_threshold);

  // create submaps that don't meet the requirements of the predicate
  const HsfmSubmap smap_not_enough_views = generateBasicSubmap(views_threshold - 1, 100);
  HsfmSubmap smap_parent_already = generateBasicSubmap(views_threshold * 2, 100);
  smap_parent_already.is_parent = true;
  smap_parent_already.children_submaps = {1,2};

  EXPECT_FALSE(is_partitionable(smap_not_enough_views));
  EXPECT_FALSE(is_partitionable(smap_parent_already));

  const HsfmSubmap smap_good = generateBasicSubmap(views_threshold * 2, 100);
  EXPECT_TRUE(is_partitionable(smap_good));
}

TEST(SubmapPartitionablePredicate, works_asPredicateOnAlgorithms)
{
  const int tracks_threshold = 100;

  openMVG::sfm::SubmapTracksThresholdChecker is_partitionable(tracks_threshold);
  const HsfmSubmap smap1 = generateBasicSubmap(20, tracks_threshold * 2);
  const HsfmSubmap smap2 = generateBasicSubmap(20, tracks_threshold * 2);
  const HsfmSubmap smap3 = generateBasicSubmap(20, tracks_threshold * 2);
  const HsfmSubmap smap_not_enough_tracks = generateBasicSubmap(20, tracks_threshold - 1);
  const HsfmSubmaps submaps = {{1,smap1}, {2,smap2}, {3,smap3}, {4, smap_not_enough_tracks}};

  EXPECT_FALSE(std::all_of(submaps.begin(), submaps.end(), is_partitionable));
  EXPECT_TRUE(std::any_of(submaps.begin(), submaps.end(), is_partitionable));
}

/* ************************************************************************* */
int main() {
 TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
