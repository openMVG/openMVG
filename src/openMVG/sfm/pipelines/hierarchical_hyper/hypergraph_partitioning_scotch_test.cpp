#include "openMVG/sfm/pipelines/hierarchical_hyper/hypergraph_partitioning_scotch.hpp"

#include "testing/testing.h"

TEST(ScotchPartitionHyperGraph, emptyHyperGraph_returnsFalse)
{
  const std::map<std::set<openMVG::IndexT>, std::set<size_t>> hyper_graph;
  std::pair<std::set<openMVG::IndexT>, std::set<openMVG::IndexT>> view_id_partitions;

  const bool result = ScotchPartitionHyperGraph(hyper_graph, view_id_partitions);

  EXPECT_FALSE(result);
}

TEST(ScotchPartitionHyperGraph, twoNonOverlappingHyperEdges_goInSeparatePartitions)
{
  // create only two non-overlapping hyper edges
  const std::set<openMVG::IndexT> first_hyper_edge = {0,1,2};
  const std::set<openMVG::IndexT> second_hyper_edge = {3,4};

  // create hypergraph
  std::map<std::set<openMVG::IndexT>, std::set<size_t>> hyper_graph;
  hyper_graph[first_hyper_edge] = {0}; // content of "track" ids don't matter here
  hyper_graph[second_hyper_edge] = {1};

  // partition
  std::pair<std::set<openMVG::IndexT>, std::set<openMVG::IndexT>> view_id_partitions;
  const bool partitioning_succesful = ScotchPartitionHyperGraph(hyper_graph, view_id_partitions);
  EXPECT_TRUE(partitioning_succesful);

  bool first_hyper_edge_found_in_partition =
      view_id_partitions.first == first_hyper_edge || view_id_partitions.second == first_hyper_edge;
  bool second_hyper_edge_found_in_partition =
      view_id_partitions.first == second_hyper_edge || view_id_partitions.second == second_hyper_edge;
  EXPECT_TRUE(first_hyper_edge_found_in_partition && second_hyper_edge_found_in_partition);
}

/* ************************************************************************* */
int main() {
 TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
