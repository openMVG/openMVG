// Ceres Solver - A fast non-linear least squares minimizer
// Copyright 2023 Google Inc. All rights reserved.
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
// Author: kushalav@google.com (Avanish Kushal)
//         sameeragarwal@google.com (Sameer Agarwal)

#include "ceres/visibility.h"

#include <memory>
#include <set>
#include <vector>

#include "ceres/block_structure.h"
#include "ceres/graph.h"
#include "glog/logging.h"
#include "gtest/gtest.h"

namespace ceres::internal {

class VisibilityTest : public ::testing::Test {};

TEST(VisibilityTest, SimpleMatrix) {
  //   A = [1 0 0 0 0 1
  //        1 0 0 1 0 0
  //        0 1 1 0 0 0
  //        0 1 0 0 1 0]

  int num_cols = 6;
  int num_eliminate_blocks = 2;
  CompressedRowBlockStructure bs;

  // Row 1
  {
    bs.rows.emplace_back();
    CompressedRow& row = bs.rows.back();
    row.block.size = 2;
    row.block.position = 0;
    row.cells.emplace_back(0, 0);
    row.cells.emplace_back(5, 0);
  }

  // Row 2
  {
    bs.rows.emplace_back();
    CompressedRow& row = bs.rows.back();
    row.block.size = 2;
    row.block.position = 2;
    row.cells.emplace_back(0, 1);
    row.cells.emplace_back(3, 1);
  }

  // Row 3
  {
    bs.rows.emplace_back();
    CompressedRow& row = bs.rows.back();
    row.block.size = 2;
    row.block.position = 4;
    row.cells.emplace_back(1, 2);
    row.cells.emplace_back(2, 2);
  }

  // Row 4
  {
    bs.rows.emplace_back();
    CompressedRow& row = bs.rows.back();
    row.block.size = 2;
    row.block.position = 6;
    row.cells.emplace_back(1, 3);
    row.cells.emplace_back(4, 3);
  }
  bs.cols.resize(num_cols);

  std::vector<std::set<int>> visibility;
  ComputeVisibility(bs, num_eliminate_blocks, &visibility);
  ASSERT_EQ(visibility.size(), num_cols - num_eliminate_blocks);
  for (const auto& visible : visibility) {
    ASSERT_EQ(visible.size(), 1);
  }

  std::unique_ptr<WeightedGraph<int>> graph(
      CreateSchurComplementGraph(visibility));
  EXPECT_EQ(graph->vertices().size(), visibility.size());
  for (int i = 0; i < visibility.size(); ++i) {
    EXPECT_EQ(graph->VertexWeight(i), 1.0);
  }

  for (int i = 0; i < visibility.size(); ++i) {
    for (int j = i; j < visibility.size(); ++j) {
      double edge_weight = 0.0;
      if ((i == 1 && j == 3) || (i == 0 && j == 2) || (i == j)) {
        edge_weight = 1.0;
      }

      EXPECT_EQ(graph->EdgeWeight(i, j), edge_weight)
          << "Edge: " << i << " " << j << " weight: " << graph->EdgeWeight(i, j)
          << " expected weight: " << edge_weight;
    }
  }
}

TEST(VisibilityTest, NoEBlocks) {
  //   A = [1 0 0 0 0 0
  //        1 0 0 0 0 0
  //        0 1 0 0 0 0
  //        0 1 0 0 0 0]

  int num_cols = 6;
  int num_eliminate_blocks = 2;
  CompressedRowBlockStructure bs;

  // Row 1
  {
    bs.rows.emplace_back();
    CompressedRow& row = bs.rows.back();
    row.block.size = 2;
    row.block.position = 0;
    row.cells.emplace_back(0, 0);
  }

  // Row 2
  {
    bs.rows.emplace_back();
    CompressedRow& row = bs.rows.back();
    row.block.size = 2;
    row.block.position = 2;
    row.cells.emplace_back(0, 1);
  }

  // Row 3
  {
    bs.rows.emplace_back();
    CompressedRow& row = bs.rows.back();
    row.block.size = 2;
    row.block.position = 4;
    row.cells.emplace_back(1, 2);
  }

  // Row 4
  {
    bs.rows.emplace_back();
    CompressedRow& row = bs.rows.back();
    row.block.size = 2;
    row.block.position = 6;
    row.cells.emplace_back(1, 3);
  }
  bs.cols.resize(num_cols);

  std::vector<std::set<int>> visibility;
  ComputeVisibility(bs, num_eliminate_blocks, &visibility);
  ASSERT_EQ(visibility.size(), num_cols - num_eliminate_blocks);
  for (const auto& visible : visibility) {
    ASSERT_EQ(visible.size(), 0);
  }

  std::unique_ptr<WeightedGraph<int>> graph(
      CreateSchurComplementGraph(visibility));
  EXPECT_EQ(graph->vertices().size(), visibility.size());
  for (int i = 0; i < visibility.size(); ++i) {
    EXPECT_EQ(graph->VertexWeight(i), 1.0);
  }

  for (int i = 0; i < visibility.size(); ++i) {
    for (int j = i; j < visibility.size(); ++j) {
      double edge_weight = 0.0;
      if (i == j) {
        edge_weight = 1.0;
      }
      EXPECT_EQ(graph->EdgeWeight(i, j), edge_weight)
          << "Edge: " << i << " " << j << " weight: " << graph->EdgeWeight(i, j)
          << " expected weight: " << edge_weight;
    }
  }
}

}  // namespace ceres::internal
