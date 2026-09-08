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
// Author: sameeragarwal@google.com (Sameer Agarwal)

#include "ceres/inner_product_computer.h"

#include <memory>
#include <numeric>
#include <random>

#include "Eigen/SparseCore"
#include "ceres/block_sparse_matrix.h"
#include "ceres/internal/eigen.h"
#include "ceres/triplet_sparse_matrix.h"
#include "glog/logging.h"
#include "gtest/gtest.h"

namespace ceres {
namespace internal {

#define COMPUTE_AND_COMPARE                                                   \
  {                                                                           \
    inner_product_computer->Compute();                                        \
    CompressedRowSparseMatrix* actual_product_crsm =                          \
        inner_product_computer->mutable_result();                             \
    Matrix actual_inner_product =                                             \
        Eigen::Map<Eigen::SparseMatrix<double, Eigen::ColMajor>>(             \
            actual_product_crsm->num_rows(),                                  \
            actual_product_crsm->num_rows(),                                  \
            actual_product_crsm->num_nonzeros(),                              \
            actual_product_crsm->mutable_rows(),                              \
            actual_product_crsm->mutable_cols(),                              \
            actual_product_crsm->mutable_values());                           \
    EXPECT_EQ(actual_inner_product.rows(), actual_inner_product.cols());      \
    EXPECT_EQ(expected_inner_product.rows(), expected_inner_product.cols());  \
    EXPECT_EQ(actual_inner_product.rows(), expected_inner_product.rows());    \
    Matrix expected_t, actual_t;                                              \
    if (actual_product_crsm->storage_type() ==                                \
        CompressedRowSparseMatrix::StorageType::LOWER_TRIANGULAR) {           \
      expected_t = expected_inner_product.triangularView<Eigen::Upper>();     \
      actual_t = actual_inner_product.triangularView<Eigen::Upper>();         \
    } else {                                                                  \
      expected_t = expected_inner_product.triangularView<Eigen::Lower>();     \
      actual_t = actual_inner_product.triangularView<Eigen::Lower>();         \
    }                                                                         \
    EXPECT_LE((expected_t - actual_t).norm(),                                 \
              100 * std::numeric_limits<double>::epsilon() * actual_t.norm()) \
        << "expected: \n"                                                     \
        << expected_t << "\nactual: \n"                                       \
        << actual_t;                                                          \
  }

TEST(InnerProductComputer, NormalOperation) {
  const int kMaxNumRowBlocks = 10;
  const int kMaxNumColBlocks = 10;
  const int kNumTrials = 10;
  std::mt19937 prng;
  std::uniform_real_distribution<double> distribution(0.01, 1.0);

  // Create a random matrix, compute its outer product using Eigen and
  // ComputeOuterProduct. Convert both matrices to dense matrices and
  // compare their upper triangular parts.
  for (int num_row_blocks = 1; num_row_blocks < kMaxNumRowBlocks;
       ++num_row_blocks) {
    for (int num_col_blocks = 1; num_col_blocks < kMaxNumColBlocks;
         ++num_col_blocks) {
      for (int trial = 0; trial < kNumTrials; ++trial) {
        BlockSparseMatrix::RandomMatrixOptions options;
        options.num_row_blocks = num_row_blocks;
        options.num_col_blocks = num_col_blocks;
        options.min_row_block_size = 1;
        options.max_row_block_size = 5;
        options.min_col_block_size = 1;
        options.max_col_block_size = 10;
        options.block_density = distribution(prng);

        VLOG(2) << "num row blocks: " << options.num_row_blocks;
        VLOG(2) << "num col blocks: " << options.num_col_blocks;
        VLOG(2) << "min row block size: " << options.min_row_block_size;
        VLOG(2) << "max row block size: " << options.max_row_block_size;
        VLOG(2) << "min col block size: " << options.min_col_block_size;
        VLOG(2) << "max col block size: " << options.max_col_block_size;
        VLOG(2) << "block density: " << options.block_density;

        std::unique_ptr<BlockSparseMatrix> random_matrix(
            BlockSparseMatrix::CreateRandomMatrix(options, prng));

        TripletSparseMatrix tsm(random_matrix->num_rows(),
                                random_matrix->num_cols(),
                                random_matrix->num_nonzeros());
        random_matrix->ToTripletSparseMatrix(&tsm);
        std::vector<Eigen::Triplet<double>> triplets;
        for (int i = 0; i < tsm.num_nonzeros(); ++i) {
          triplets.emplace_back(tsm.rows()[i], tsm.cols()[i], tsm.values()[i]);
        }
        Eigen::SparseMatrix<double> eigen_random_matrix(
            random_matrix->num_rows(), random_matrix->num_cols());
        eigen_random_matrix.setFromTriplets(triplets.begin(), triplets.end());
        Matrix expected_inner_product =
            eigen_random_matrix.transpose() * eigen_random_matrix;

        std::unique_ptr<InnerProductComputer> inner_product_computer;

        inner_product_computer = InnerProductComputer::Create(
            *random_matrix,
            CompressedRowSparseMatrix::StorageType::LOWER_TRIANGULAR);
        COMPUTE_AND_COMPARE;
        inner_product_computer = InnerProductComputer::Create(
            *random_matrix,
            CompressedRowSparseMatrix::StorageType::UPPER_TRIANGULAR);
        COMPUTE_AND_COMPARE;
      }
    }
  }
}

TEST(InnerProductComputer, SubMatrix) {
  const int kNumRowBlocks = 10;
  const int kNumColBlocks = 20;
  const int kNumTrials = 5;
  std::mt19937 prng;
  std::uniform_real_distribution<double> distribution(0.01, 1.0);

  // Create a random matrix, compute its outer product using Eigen and
  // ComputeInnerProductComputer. Convert both matrices to dense matrices and
  // compare their upper triangular parts.
  for (int trial = 0; trial < kNumTrials; ++trial) {
    BlockSparseMatrix::RandomMatrixOptions options;
    options.num_row_blocks = kNumRowBlocks;
    options.num_col_blocks = kNumColBlocks;
    options.min_row_block_size = 1;
    options.max_row_block_size = 5;
    options.min_col_block_size = 1;
    options.max_col_block_size = 10;
    options.block_density = distribution(prng);

    VLOG(2) << "num row blocks: " << options.num_row_blocks;
    VLOG(2) << "num col blocks: " << options.num_col_blocks;
    VLOG(2) << "min row block size: " << options.min_row_block_size;
    VLOG(2) << "max row block size: " << options.max_row_block_size;
    VLOG(2) << "min col block size: " << options.min_col_block_size;
    VLOG(2) << "max col block size: " << options.max_col_block_size;
    VLOG(2) << "block density: " << options.block_density;

    std::unique_ptr<BlockSparseMatrix> random_matrix(
        BlockSparseMatrix::CreateRandomMatrix(options, prng));

    const std::vector<CompressedRow>& row_blocks =
        random_matrix->block_structure()->rows;
    const int num_row_blocks = row_blocks.size();

    for (int start_row_block = 0; start_row_block < num_row_blocks - 1;
         ++start_row_block) {
      for (int end_row_block = start_row_block + 1;
           end_row_block < num_row_blocks;
           ++end_row_block) {
        const int start_row = row_blocks[start_row_block].block.position;
        const int end_row = row_blocks[end_row_block].block.position;

        TripletSparseMatrix tsm(random_matrix->num_rows(),
                                random_matrix->num_cols(),
                                random_matrix->num_nonzeros());
        random_matrix->ToTripletSparseMatrix(&tsm);
        std::vector<Eigen::Triplet<double>> triplets;
        for (int i = 0; i < tsm.num_nonzeros(); ++i) {
          if (tsm.rows()[i] >= start_row && tsm.rows()[i] < end_row) {
            triplets.emplace_back(
                tsm.rows()[i], tsm.cols()[i], tsm.values()[i]);
          }
        }

        Eigen::SparseMatrix<double> eigen_random_matrix(
            random_matrix->num_rows(), random_matrix->num_cols());
        eigen_random_matrix.setFromTriplets(triplets.begin(), triplets.end());

        Matrix expected_inner_product =
            eigen_random_matrix.transpose() * eigen_random_matrix;

        std::unique_ptr<InnerProductComputer> inner_product_computer;
        inner_product_computer = InnerProductComputer::Create(
            *random_matrix,
            start_row_block,
            end_row_block,
            CompressedRowSparseMatrix::StorageType::LOWER_TRIANGULAR);
        COMPUTE_AND_COMPARE;
        inner_product_computer = InnerProductComputer::Create(
            *random_matrix,
            start_row_block,
            end_row_block,
            CompressedRowSparseMatrix::StorageType::UPPER_TRIANGULAR);
        COMPUTE_AND_COMPARE;
      }
    }
  }
}

#undef COMPUTE_AND_COMPARE
}  // namespace internal
}  // namespace ceres
