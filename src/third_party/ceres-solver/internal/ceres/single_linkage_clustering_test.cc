// Ceres Solver - A fast non-linear least squares minimizer
// Copyright 2015 Google Inc. All rights reserved.
// http://ceres-solver.org/
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
// * Neither the name of Google Inc. nor the names of its contributors may be
//   used to endorse or promote products derived from this software without
//   specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Author: Sameer Agarwal (sameeragarwal@google.com)

// This include must come before any #ifndef check on Ceres compile options.
#include "ceres/internal/port.h"

#ifndef CERES_NO_SUITESPARSE

#include "ceres/single_linkage_clustering.h"

#include "ceres/collections_port.h"
#include "ceres/graph.h"
#include "gtest/gtest.h"

namespace ceres {
namespace internal {

TEST(SingleLinkageClustering, GraphHasTwoComponents) {
  WeightedGraph<int> graph;
  const int kNumVertices = 6;
  for (int i = 0; i < kNumVertices; ++i) {
    graph.AddVertex(i);
  }
  // Graph structure:
  //
  //  0-1-2-3 4-5
  graph.AddEdge(0, 1, 1.0);
  graph.AddEdge(1, 2, 1.0);
  graph.AddEdge(2, 3, 1.0);
  graph.AddEdge(4, 5, 1.0);

  SingleLinkageClusteringOptions options;
  HashMap<int, int> membership;
  ComputeSingleLinkageClustering(options, graph, &membership);
  EXPECT_EQ(membership.size(), kNumVertices);

  EXPECT_EQ(membership[1], membership[0]);
  EXPECT_EQ(membership[2], membership[0]);
  EXPECT_EQ(membership[3], membership[0]);
  EXPECT_NE(membership[4], membership[0]);
  EXPECT_NE(membership[5], membership[0]);
  EXPECT_EQ(membership[4], membership[5]);
}

TEST(SingleLinkageClustering, ComponentWithWeakLink) {
  WeightedGraph<int> graph;
  const int kNumVertices = 6;
  for (int i = 0; i < kNumVertices; ++i) {
    graph.AddVertex(i);
  }
  // Graph structure:
  //
  //  0-1-2-3 4-5
  graph.AddEdge(0, 1, 1.0);
  graph.AddEdge(1, 2, 1.0);
  graph.AddEdge(2, 3, 1.0);

  // This component should break up into two.
  graph.AddEdge(4, 5, 0.5);

  SingleLinkageClusteringOptions options;
  HashMap<int, int> membership;
  ComputeSingleLinkageClustering(options, graph, &membership);
  EXPECT_EQ(membership.size(), kNumVertices);

  EXPECT_EQ(membership[1], membership[0]);
  EXPECT_EQ(membership[2], membership[0]);
  EXPECT_EQ(membership[3], membership[0]);
  EXPECT_NE(membership[4], membership[0]);
  EXPECT_NE(membership[5], membership[0]);
  EXPECT_NE(membership[4], membership[5]);
}

TEST(SingleLinkageClustering, ComponentWithWeakLinkAndStrongLink) {
  WeightedGraph<int> graph;
  const int kNumVertices = 6;
  for (int i = 0; i < kNumVertices; ++i) {
    graph.AddVertex(i);
  }
  // Graph structure:
  //
  //  0-1-2-3 4-5
  graph.AddEdge(0, 1, 1.0);
  graph.AddEdge(1, 2, 1.0);
  graph.AddEdge(2, 3, 0.5);  // Weak link
  graph.AddEdge(0, 3, 1.0);

  // This component should break up into two.
  graph.AddEdge(4, 5, 1.0);

  SingleLinkageClusteringOptions options;
  HashMap<int, int> membership;
  ComputeSingleLinkageClustering(options, graph, &membership);
  EXPECT_EQ(membership.size(), kNumVertices);

  EXPECT_EQ(membership[1], membership[0]);
  EXPECT_EQ(membership[2], membership[0]);
  EXPECT_EQ(membership[3], membership[0]);
  EXPECT_EQ(membership[4], membership[5]);
}

}  // namespace internal
}  // namespace ceres

#endif  // CERES_NO_SUITESPARSE
