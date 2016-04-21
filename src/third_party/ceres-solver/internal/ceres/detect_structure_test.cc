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
// Author: sameeragarwal@google.com (Sameer Agarwal)

#include "Eigen/Core"
#include "glog/logging.h"
#include "gtest/gtest.h"
#include "ceres/block_structure.h"
#include "ceres/detect_structure.h"

namespace ceres {
namespace internal {

TEST(DetectStructure, EverythingStatic) {
  const int expected_row_block_size = 2;
  const int expected_e_block_size = 3;
  const int expected_f_block_size = 4;

  CompressedRowBlockStructure bs;

  bs.cols.push_back(Block());
  bs.cols.back().size = 3;
  bs.cols.back().position = 0;

  bs.cols.push_back(Block());
  bs.cols.back().size = 4;
  bs.cols.back().position = 3;

  bs.cols.push_back(Block());
  bs.cols.back().size = 4;
  bs.cols.back().position = 7;

  {
    bs.rows.push_back(CompressedRow());
    CompressedRow& row = bs.rows.back();
    row.block.size = 2;
    row.block.position = 0;
    row.cells.push_back(Cell(0, 0));
    row.cells.push_back(Cell(1, 0));
  }

  {
    bs.rows.push_back(CompressedRow());
    CompressedRow& row = bs.rows.back();
    row.block.size = 2;
    row.block.position = 2;
    row.cells.push_back(Cell(0, 0));
    row.cells.push_back(Cell(2, 0));
  }

  int row_block_size = 0;
  int e_block_size = 0;
  int f_block_size = 0;
  const int num_eliminate_blocks = 1;
  DetectStructure(bs,
                  num_eliminate_blocks,
                  &row_block_size,
                  &e_block_size,
                  &f_block_size);

  EXPECT_EQ(row_block_size, expected_row_block_size);
  EXPECT_EQ(e_block_size, expected_e_block_size);
  EXPECT_EQ(f_block_size, expected_f_block_size);
}

TEST(DetectStructure, DynamicRow) {
  const int expected_row_block_size = Eigen::Dynamic;
  const int expected_e_block_size = 3;
  const int expected_f_block_size = 4;

  CompressedRowBlockStructure bs;

  bs.cols.push_back(Block());
  bs.cols.back().size = 3;
  bs.cols.back().position = 0;

  bs.cols.push_back(Block());
  bs.cols.back().size = 4;
  bs.cols.back().position = 3;

  bs.cols.push_back(Block());
  bs.cols.back().size = 4;
  bs.cols.back().position = 7;

  {
    bs.rows.push_back(CompressedRow());
    CompressedRow& row = bs.rows.back();
    row.block.size = 2;
    row.block.position = 0;
    row.cells.push_back(Cell(0, 0));
    row.cells.push_back(Cell(1, 0));
  }

  {
    bs.rows.push_back(CompressedRow());
    CompressedRow& row = bs.rows.back();
    row.block.size = 1;
    row.block.position = 2;
    row.cells.push_back(Cell(0, 0));
    row.cells.push_back(Cell(2, 0));
  }

  int row_block_size = 0;
  int e_block_size = 0;
  int f_block_size = 0;
  const int num_eliminate_blocks = 1;
  DetectStructure(bs,
                  num_eliminate_blocks,
                  &row_block_size,
                  &e_block_size,
                  &f_block_size);

  EXPECT_EQ(row_block_size, expected_row_block_size);
  EXPECT_EQ(e_block_size, expected_e_block_size);
  EXPECT_EQ(f_block_size, expected_f_block_size);
}

TEST(DetectStructure, DynamicFBlockDifferentRows) {
  const int expected_row_block_size = 2;
  const int expected_e_block_size = 3;
  const int expected_f_block_size = Eigen::Dynamic;


  CompressedRowBlockStructure bs;

  bs.cols.push_back(Block());
  bs.cols.back().size = 3;
  bs.cols.back().position = 0;

  bs.cols.push_back(Block());
  bs.cols.back().size = 4;
  bs.cols.back().position = 3;

  bs.cols.push_back(Block());
  bs.cols.back().size = 3;
  bs.cols.back().position = 7;

  {
    bs.rows.push_back(CompressedRow());
    CompressedRow& row = bs.rows.back();
    row.block.size = 2;
    row.block.position = 0;
    row.cells.push_back(Cell(0, 0));
    row.cells.push_back(Cell(1, 0));
  }

  {
    bs.rows.push_back(CompressedRow());
    CompressedRow& row = bs.rows.back();
    row.block.size = 2;
    row.block.position = 2;
    row.cells.push_back(Cell(0, 0));
    row.cells.push_back(Cell(2, 0));
  }

  int row_block_size = 0;
  int e_block_size = 0;
  int f_block_size = 0;
  const int num_eliminate_blocks = 1;
  DetectStructure(bs,
                  num_eliminate_blocks,
                  &row_block_size,
                  &e_block_size,
                  &f_block_size);

  EXPECT_EQ(row_block_size, expected_row_block_size);
  EXPECT_EQ(e_block_size, expected_e_block_size);
  EXPECT_EQ(f_block_size, expected_f_block_size);
}

TEST(DetectStructure, DynamicEBlock) {
  const int expected_row_block_size = 2;
  const int expected_e_block_size = Eigen::Dynamic;
  const int expected_f_block_size = 3;

  CompressedRowBlockStructure bs;

  bs.cols.push_back(Block());
  bs.cols.back().size = 3;
  bs.cols.back().position = 0;

  bs.cols.push_back(Block());
  bs.cols.back().size = 4;
  bs.cols.back().position = 3;

  bs.cols.push_back(Block());
  bs.cols.back().size = 3;
  bs.cols.back().position = 7;

  {
    bs.rows.push_back(CompressedRow());
    CompressedRow& row = bs.rows.back();
    row.block.size = 2;
    row.block.position = 0;
    row.cells.push_back(Cell(0, 0));
    row.cells.push_back(Cell(2, 0));
  }

  {
    bs.rows.push_back(CompressedRow());
    CompressedRow& row = bs.rows.back();
    row.block.size = 2;
    row.block.position = 2;
    row.cells.push_back(Cell(1, 0));
    row.cells.push_back(Cell(2, 0));
  }

  int row_block_size = 0;
  int e_block_size = 0;
  int f_block_size = 0;
  const int num_eliminate_blocks = 2;
  DetectStructure(bs,
                  num_eliminate_blocks,
                  &row_block_size,
                  &e_block_size,
                  &f_block_size);

  EXPECT_EQ(row_block_size, expected_row_block_size);
  EXPECT_EQ(e_block_size, expected_e_block_size);
  EXPECT_EQ(f_block_size, expected_f_block_size);
}

TEST(DetectStructure, DynamicFBlockSameRow) {
  const int expected_row_block_size = 2;
  const int expected_e_block_size = 3;
  const int expected_f_block_size = Eigen::Dynamic;

  CompressedRowBlockStructure bs;

  bs.cols.push_back(Block());
  bs.cols.back().size = 3;
  bs.cols.back().position = 0;

  bs.cols.push_back(Block());
  bs.cols.back().size = 4;
  bs.cols.back().position = 3;

  bs.cols.push_back(Block());
  bs.cols.back().size = 3;
  bs.cols.back().position = 7;

  {
    bs.rows.push_back(CompressedRow());
    CompressedRow& row = bs.rows.back();
    row.block.size = 2;
    row.block.position = 0;
    row.cells.push_back(Cell(0, 0));
    row.cells.push_back(Cell(1, 0));
    row.cells.push_back(Cell(2, 0));
  }

  int row_block_size = 0;
  int e_block_size = 0;
  int f_block_size = 0;
  const int num_eliminate_blocks = 1;
  DetectStructure(bs,
                  num_eliminate_blocks,
                  &row_block_size,
                  &e_block_size,
                  &f_block_size);

  EXPECT_EQ(row_block_size, expected_row_block_size);
  EXPECT_EQ(e_block_size, expected_e_block_size);
  EXPECT_EQ(f_block_size, expected_f_block_size);
}

}  // namespace internal
}  // namespace ceres
