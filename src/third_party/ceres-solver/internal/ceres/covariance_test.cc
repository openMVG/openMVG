// Ceres Solver - A fast non-linear least squares minimizer
// Copyright 2013 Google Inc. All rights reserved.
// http://code.google.com/p/ceres-solver/
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

#include "ceres/covariance.h"

#include <algorithm>
#include <cmath>
#include "ceres/compressed_row_sparse_matrix.h"
#include "ceres/cost_function.h"
#include "ceres/covariance_impl.h"
#include "ceres/local_parameterization.h"
#include "ceres/map_util.h"
#include "ceres/problem_impl.h"
#include "gtest/gtest.h"

namespace ceres {
namespace internal {

TEST(CovarianceImpl, ComputeCovarianceSparsity) {
  double parameters[10];

  double* block1 = parameters;
  double* block2 = block1 + 1;
  double* block3 = block2 + 2;
  double* block4 = block3 + 3;

  ProblemImpl problem;

  // Add in random order
  problem.AddParameterBlock(block1, 1);
  problem.AddParameterBlock(block4, 4);
  problem.AddParameterBlock(block3, 3);
  problem.AddParameterBlock(block2, 2);

  // Sparsity pattern
  //
  //  x 0 0 0 0 0 x x x x
  //  0 x x x x x 0 0 0 0
  //  0 x x x x x 0 0 0 0
  //  0 0 0 x x x 0 0 0 0
  //  0 0 0 x x x 0 0 0 0
  //  0 0 0 x x x 0 0 0 0
  //  0 0 0 0 0 0 x x x x
  //  0 0 0 0 0 0 x x x x
  //  0 0 0 0 0 0 x x x x
  //  0 0 0 0 0 0 x x x x

  int expected_rows[] = {0, 5, 10, 15, 18, 21, 24, 28, 32, 36, 40};
  int expected_cols[] = {0, 6, 7, 8, 9,
                         1, 2, 3, 4, 5,
                         1, 2, 3, 4, 5,
                         3, 4, 5,
                         3, 4, 5,
                         3, 4, 5,
                         6, 7, 8, 9,
                         6, 7, 8, 9,
                         6, 7, 8, 9,
                         6, 7, 8, 9};


  vector<pair<const double*, const double*> > covariance_blocks;
  covariance_blocks.push_back(make_pair(block1, block1));
  covariance_blocks.push_back(make_pair(block4, block4));
  covariance_blocks.push_back(make_pair(block2, block2));
  covariance_blocks.push_back(make_pair(block3, block3));
  covariance_blocks.push_back(make_pair(block2, block3));
  covariance_blocks.push_back(make_pair(block4, block1));  // reversed

  Covariance::Options options;
  CovarianceImpl covariance_impl(options);
  EXPECT_TRUE(covariance_impl
              .ComputeCovarianceSparsity(covariance_blocks, &problem));

  const CompressedRowSparseMatrix* crsm = covariance_impl.covariance_matrix();

  EXPECT_EQ(crsm->num_rows(), 10);
  EXPECT_EQ(crsm->num_cols(), 10);
  EXPECT_EQ(crsm->num_nonzeros(), 40);

  const int* rows = crsm->rows();
  for (int r = 0; r < crsm->num_rows() + 1; ++r) {
    EXPECT_EQ(rows[r], expected_rows[r])
        << r << " "
        << rows[r] << " "
        << expected_rows[r];
  }

  const int* cols = crsm->cols();
  for (int c = 0; c < crsm->num_nonzeros(); ++c) {
    EXPECT_EQ(cols[c], expected_cols[c])
        << c << " "
        << cols[c] << " "
        << expected_cols[c];
  }
}


class UnaryCostFunction: public CostFunction {
 public:
  UnaryCostFunction(const int num_residuals,
                    const int16 parameter_block_size,
                    const double* jacobian)
      : jacobian_(jacobian, jacobian + num_residuals * parameter_block_size) {
    set_num_residuals(num_residuals);
    mutable_parameter_block_sizes()->push_back(parameter_block_size);
  }

  virtual bool Evaluate(double const* const* parameters,
                        double* residuals,
                        double** jacobians) const {
    for (int i = 0; i < num_residuals(); ++i) {
      residuals[i] = 1;
    }

    if (jacobians == NULL) {
      return true;
    }

    if (jacobians[0] != NULL) {
      copy(jacobian_.begin(), jacobian_.end(), jacobians[0]);
    }

    return true;
  }

 private:
  vector<double> jacobian_;
};


class BinaryCostFunction: public CostFunction {
 public:
  BinaryCostFunction(const int num_residuals,
                     const int16 parameter_block1_size,
                     const int16 parameter_block2_size,
                     const double* jacobian1,
                     const double* jacobian2)
      : jacobian1_(jacobian1,
                   jacobian1 + num_residuals * parameter_block1_size),
        jacobian2_(jacobian2,
                   jacobian2 + num_residuals * parameter_block2_size) {
    set_num_residuals(num_residuals);
    mutable_parameter_block_sizes()->push_back(parameter_block1_size);
    mutable_parameter_block_sizes()->push_back(parameter_block2_size);
  }

  virtual bool Evaluate(double const* const* parameters,
                        double* residuals,
                        double** jacobians) const {
    for (int i = 0; i < num_residuals(); ++i) {
      residuals[i] = 2;
    }

    if (jacobians == NULL) {
      return true;
    }

    if (jacobians[0] != NULL) {
      copy(jacobian1_.begin(), jacobian1_.end(), jacobians[0]);
    }

    if (jacobians[1] != NULL) {
      copy(jacobian2_.begin(), jacobian2_.end(), jacobians[1]);
    }

    return true;
  }

 private:
  vector<double> jacobian1_;
  vector<double> jacobian2_;
};

// x_plus_delta = delta * x;
class PolynomialParameterization : public LocalParameterization {
 public:
  virtual ~PolynomialParameterization() {}

  virtual bool Plus(const double* x,
                    const double* delta,
                    double* x_plus_delta) const {
    x_plus_delta[0] = delta[0] * x[0];
    x_plus_delta[1] = delta[0] * x[1];
    return true;
  }

  virtual bool ComputeJacobian(const double* x, double* jacobian) const {
    jacobian[0] = x[0];
    jacobian[1] = x[1];
    return true;
  }

  virtual int GlobalSize() const { return 2; }
  virtual int LocalSize() const { return 1; }
};

class CovarianceTest : public ::testing::Test {
 protected:
  virtual void SetUp() {
    double* x = parameters_;
    double* y = x + 2;
    double* z = y + 3;

    x[0] = 1;
    x[1] = 1;
    y[0] = 2;
    y[1] = 2;
    y[2] = 2;
    z[0] = 3;

    {
      double jacobian[] = { 1.0, 0.0, 0.0, 1.0};
      problem_.AddResidualBlock(new UnaryCostFunction(2, 2, jacobian), NULL, x);
    }

    {
      double jacobian[] = { 2.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 2.0 };
      problem_.AddResidualBlock(new UnaryCostFunction(3, 3, jacobian), NULL, y);
    }

    {
      double jacobian = 5.0;
      problem_.AddResidualBlock(new UnaryCostFunction(1, 1, &jacobian), NULL, z);
    }

    {
      double jacobian1[] = { 1.0, 2.0, 3.0 };
      double jacobian2[] = { -5.0, -6.0 };
      problem_.AddResidualBlock(
          new BinaryCostFunction(1, 3, 2, jacobian1, jacobian2),
          NULL,
          y,
          x);
    }

    {
      double jacobian1[] = {2.0 };
      double jacobian2[] = { 3.0, -2.0 };
      problem_.AddResidualBlock(
          new BinaryCostFunction(1, 1, 2, jacobian1, jacobian2),
          NULL,
          z,
          x);
    }

    all_covariance_blocks_.push_back(make_pair(x, x));
    all_covariance_blocks_.push_back(make_pair(y, y));
    all_covariance_blocks_.push_back(make_pair(z, z));
    all_covariance_blocks_.push_back(make_pair(x, y));
    all_covariance_blocks_.push_back(make_pair(x, z));
    all_covariance_blocks_.push_back(make_pair(y, z));

    column_bounds_[x] = make_pair(0, 2);
    column_bounds_[y] = make_pair(2, 5);
    column_bounds_[z] = make_pair(5, 6);
  }

  void ComputeAndCompareCovarianceBlocks(const Covariance::Options& options,
                                         const double* expected_covariance) {
    // Generate all possible combination of block pairs and check if the
    // covariance computation is correct.
    for (int i = 1; i <= 64; ++i) {
      vector<pair<const double*, const double*> > covariance_blocks;
      if (i & 1) {
        covariance_blocks.push_back(all_covariance_blocks_[0]);
      }

      if (i & 2) {
        covariance_blocks.push_back(all_covariance_blocks_[1]);
      }

      if (i & 4) {
        covariance_blocks.push_back(all_covariance_blocks_[2]);
      }

      if (i & 8) {
        covariance_blocks.push_back(all_covariance_blocks_[3]);
      }

      if (i & 16) {
        covariance_blocks.push_back(all_covariance_blocks_[4]);
      }

      if (i & 32) {
        covariance_blocks.push_back(all_covariance_blocks_[5]);
      }

      Covariance covariance(options);
      EXPECT_TRUE(covariance.Compute(covariance_blocks, &problem_));

      for (int i = 0; i < covariance_blocks.size(); ++i) {
        const double* block1 = covariance_blocks[i].first;
        const double* block2 = covariance_blocks[i].second;
        // block1, block2
        GetCovarianceBlockAndCompare(block1, block2, covariance, expected_covariance);
        // block2, block1
        GetCovarianceBlockAndCompare(block2, block1, covariance, expected_covariance);
      }
    }
  }

  void GetCovarianceBlockAndCompare(const double* block1,
                                    const double* block2,
                                    const Covariance& covariance,
                                    const double* expected_covariance) {
    const int row_begin = FindOrDie(column_bounds_, block1).first;
    const int row_end = FindOrDie(column_bounds_, block1).second;
    const int col_begin = FindOrDie(column_bounds_, block2).first;
    const int col_end = FindOrDie(column_bounds_, block2).second;

    Matrix actual(row_end - row_begin, col_end - col_begin);
    EXPECT_TRUE(covariance.GetCovarianceBlock(block1,
                                              block2,
                                              actual.data()));

    ConstMatrixRef expected(expected_covariance, 6, 6);
    double diff_norm = (expected.block(row_begin,
                                       col_begin,
                                       row_end - row_begin,
                                       col_end - col_begin) - actual).norm();
    diff_norm /= (row_end - row_begin) * (col_end - col_begin);

    const double kTolerance = 1e-5;
    EXPECT_NEAR(diff_norm, 0.0, kTolerance)
        << "rows: " << row_begin << " " << row_end << "  "
        << "cols: " << col_begin << " " << col_end << "  "
        << "\n\n expected: \n " << expected.block(row_begin,
                                                  col_begin,
                                                  row_end - row_begin,
                                                  col_end - col_begin)
        << "\n\n actual: \n " << actual
        << "\n\n full expected: \n" << expected;
  }

  double parameters_[10];
  Problem problem_;
  vector<pair<const double*, const double*> > all_covariance_blocks_;
  map<const double*, pair<int, int> > column_bounds_;
};


TEST_F(CovarianceTest, NormalBehavior) {
  // J
  //
  //   1  0  0  0  0  0
  //   0  1  0  0  0  0
  //   0  0  2  0  0  0
  //   0  0  0  2  0  0
  //   0  0  0  0  2  0
  //   0  0  0  0  0  5
  //  -5 -6  1  2  3  0
  //   3 -2  0  0  0  2

  // J'J
  //
  //   35  24 -5 -10 -15  6
  //   24  41 -6 -12 -18 -4
  //   -5  -6  5   2   3  0
  //  -10 -12  2   8   6  0
  //  -15 -18  3   6  13  0
  //    6  -4  0   0   0 29

  // inv(J'J) computed using octave.
  double expected_covariance[] = {
     7.0747e-02,  -8.4923e-03,   1.6821e-02,   3.3643e-02,   5.0464e-02,  -1.5809e-02,  // NOLINT
    -8.4923e-03,   8.1352e-02,   2.4758e-02,   4.9517e-02,   7.4275e-02,   1.2978e-02,  // NOLINT
     1.6821e-02,   2.4758e-02,   2.4904e-01,  -1.9271e-03,  -2.8906e-03,  -6.5325e-05,  // NOLINT
     3.3643e-02,   4.9517e-02,  -1.9271e-03,   2.4615e-01,  -5.7813e-03,  -1.3065e-04,  // NOLINT
     5.0464e-02,   7.4275e-02,  -2.8906e-03,  -5.7813e-03,   2.4133e-01,  -1.9598e-04,  // NOLINT
    -1.5809e-02,   1.2978e-02,  -6.5325e-05,  -1.3065e-04,  -1.9598e-04,   3.9544e-02,  // NOLINT
  };

  Covariance::Options options;

#ifndef CERES_NO_SUITESPARSE
  options.algorithm_type = SPARSE_CHOLESKY;
  ComputeAndCompareCovarianceBlocks(options, expected_covariance);

  options.algorithm_type = SPARSE_QR;
  ComputeAndCompareCovarianceBlocks(options, expected_covariance);
#endif

  options.algorithm_type = DENSE_SVD;
  ComputeAndCompareCovarianceBlocks(options, expected_covariance);
}

#ifdef CERES_USE_OPENMP

TEST_F(CovarianceTest, ThreadedNormalBehavior) {
  // J
  //
  //   1  0  0  0  0  0
  //   0  1  0  0  0  0
  //   0  0  2  0  0  0
  //   0  0  0  2  0  0
  //   0  0  0  0  2  0
  //   0  0  0  0  0  5
  //  -5 -6  1  2  3  0
  //   3 -2  0  0  0  2

  // J'J
  //
  //   35  24 -5 -10 -15  6
  //   24  41 -6 -12 -18 -4
  //   -5  -6  5   2   3  0
  //  -10 -12  2   8   6  0
  //  -15 -18  3   6  13  0
  //    6  -4  0   0   0 29

  // inv(J'J) computed using octave.
  double expected_covariance[] = {
     7.0747e-02,  -8.4923e-03,   1.6821e-02,   3.3643e-02,   5.0464e-02,  -1.5809e-02,  // NOLINT
    -8.4923e-03,   8.1352e-02,   2.4758e-02,   4.9517e-02,   7.4275e-02,   1.2978e-02,  // NOLINT
     1.6821e-02,   2.4758e-02,   2.4904e-01,  -1.9271e-03,  -2.8906e-03,  -6.5325e-05,  // NOLINT
     3.3643e-02,   4.9517e-02,  -1.9271e-03,   2.4615e-01,  -5.7813e-03,  -1.3065e-04,  // NOLINT
     5.0464e-02,   7.4275e-02,  -2.8906e-03,  -5.7813e-03,   2.4133e-01,  -1.9598e-04,  // NOLINT
    -1.5809e-02,   1.2978e-02,  -6.5325e-05,  -1.3065e-04,  -1.9598e-04,   3.9544e-02,  // NOLINT
  };

  Covariance::Options options;
  options.num_threads = 4;

#ifndef CERES_NO_SUITESPARSE
  options.algorithm_type = SPARSE_CHOLESKY;
  ComputeAndCompareCovarianceBlocks(options, expected_covariance);

  options.algorithm_type = SPARSE_QR;
  ComputeAndCompareCovarianceBlocks(options, expected_covariance);
#endif

  options.algorithm_type = DENSE_SVD;
  ComputeAndCompareCovarianceBlocks(options, expected_covariance);
}

#endif  // CERES_USE_OPENMP

TEST_F(CovarianceTest, ConstantParameterBlock) {
  problem_.SetParameterBlockConstant(parameters_);

  // J
  //
  //  0  0  0  0  0  0
  //  0  0  0  0  0  0
  //  0  0  2  0  0  0
  //  0  0  0  2  0  0
  //  0  0  0  0  2  0
  //  0  0  0  0  0  5
  //  0  0  1  2  3  0
  //  0  0  0  0  0  2

  // J'J
  //
  //  0  0  0  0  0  0
  //  0  0  0  0  0  0
  //  0  0  5  2  3  0
  //  0  0  2  8  6  0
  //  0  0  3  6 13  0
  //  0  0  0  0  0 29

  // pinv(J'J) computed using octave.
  double expected_covariance[] = {
              0,            0,            0,            0,            0,            0,  // NOLINT
              0,            0,            0,            0,            0,            0,  // NOLINT
              0,            0,      0.23611,     -0.02778,     -0.04167,     -0.00000,  // NOLINT
              0,            0,     -0.02778,      0.19444,     -0.08333,     -0.00000,  // NOLINT
              0,            0,     -0.04167,     -0.08333,      0.12500,     -0.00000,  // NOLINT
              0,            0,     -0.00000,     -0.00000,     -0.00000,      0.03448   // NOLINT
  };

  Covariance::Options options;

#ifndef CERES_NO_SUITESPARSE
  options.algorithm_type = SPARSE_CHOLESKY;
  ComputeAndCompareCovarianceBlocks(options, expected_covariance);

  options.algorithm_type = SPARSE_QR;
  ComputeAndCompareCovarianceBlocks(options, expected_covariance);
#endif

  options.algorithm_type = DENSE_SVD;
  ComputeAndCompareCovarianceBlocks(options, expected_covariance);
}

TEST_F(CovarianceTest, LocalParameterization) {
  double* x = parameters_;
  double* y = x + 2;

  problem_.SetParameterization(x, new PolynomialParameterization);

  vector<int> subset;
  subset.push_back(2);
  problem_.SetParameterization(y, new SubsetParameterization(3, subset));

  // Raw Jacobian: J
  //
  //   1   0  0  0  0  0
  //   0   1  0  0  0  0
  //   0   0  2  0  0  0
  //   0   0  0  2  0  0
  //   0   0  0  0  0  0
  //   0   0  0  0  0  5
  //  -5  -6  1  2  0  0
  //   3  -2  0  0  0  2

  // Global to local jacobian: A
  //
  //
  //  1   0   0   0   0
  //  1   0   0   0   0
  //  0   1   0   0   0
  //  0   0   1   0   0
  //  0   0   0   1   0
  //  0   0   0   0   1

  // A * pinv((J*A)'*(J*A)) * A'
  // Computed using octave.
  double expected_covariance[] = {
    0.01766,   0.01766,   0.02158,   0.04316,   0.00000,  -0.00122,
    0.01766,   0.01766,   0.02158,   0.04316,   0.00000,  -0.00122,
    0.02158,   0.02158,   0.24860,  -0.00281,   0.00000,  -0.00149,
    0.04316,   0.04316,  -0.00281,   0.24439,   0.00000,  -0.00298,
    0.00000,   0.00000,   0.00000,   0.00000,   0.00000,   0.00000,
   -0.00122,  -0.00122,  -0.00149,  -0.00298,   0.00000,   0.03457
  };

  Covariance::Options options;

#ifndef CERES_NO_SUITESPARSE
  options.algorithm_type = SPARSE_CHOLESKY;
  ComputeAndCompareCovarianceBlocks(options, expected_covariance);

  options.algorithm_type = SPARSE_QR;
  ComputeAndCompareCovarianceBlocks(options, expected_covariance);
#endif

  options.algorithm_type = DENSE_SVD;
  ComputeAndCompareCovarianceBlocks(options, expected_covariance);
}


TEST_F(CovarianceTest, TruncatedRank) {
  // J
  //
  //   1  0  0  0  0  0
  //   0  1  0  0  0  0
  //   0  0  2  0  0  0
  //   0  0  0  2  0  0
  //   0  0  0  0  2  0
  //   0  0  0  0  0  5
  //  -5 -6  1  2  3  0
  //   3 -2  0  0  0  2

  // J'J
  //
  //   35  24 -5 -10 -15  6
  //   24  41 -6 -12 -18 -4
  //   -5  -6  5   2   3  0
  //  -10 -12  2   8   6  0
  //  -15 -18  3   6  13  0
  //    6  -4  0   0   0 29

  // 3.4142 is the smallest eigen value of J'J. The following matrix
  // was obtained by dropping the eigenvector corresponding to this
  // eigenvalue.
  double expected_covariance[] = {
     5.4135e-02,  -3.5121e-02,   1.7257e-04,   3.4514e-04,   5.1771e-04,  -1.6076e-02,
    -3.5121e-02,   3.8667e-02,  -1.9288e-03,  -3.8576e-03,  -5.7864e-03,   1.2549e-02,
     1.7257e-04,  -1.9288e-03,   2.3235e-01,  -3.5297e-02,  -5.2946e-02,  -3.3329e-04,
     3.4514e-04,  -3.8576e-03,  -3.5297e-02,   1.7941e-01,  -1.0589e-01,  -6.6659e-04,
     5.1771e-04,  -5.7864e-03,  -5.2946e-02,  -1.0589e-01,   9.1162e-02,  -9.9988e-04,
    -1.6076e-02,   1.2549e-02,  -3.3329e-04,  -6.6659e-04,  -9.9988e-04,   3.9539e-02
  };


  {
    Covariance::Options options;
    options.algorithm_type = DENSE_SVD;
    // Force dropping of the smallest eigenvector.
    options.null_space_rank = 1;
    ComputeAndCompareCovarianceBlocks(options, expected_covariance);
  }

  {
    Covariance::Options options;
    options.algorithm_type = DENSE_SVD;
    // Force dropping of the smallest eigenvector via the ratio but
    // automatic truncation.
    options.min_reciprocal_condition_number = 0.044494;
    options.null_space_rank = -1;
    ComputeAndCompareCovarianceBlocks(options, expected_covariance);
  }
}

class RankDeficientCovarianceTest : public CovarianceTest {
 protected:
  virtual void SetUp() {
    double* x = parameters_;
    double* y = x + 2;
    double* z = y + 3;

    {
      double jacobian[] = { 1.0, 0.0, 0.0, 1.0};
      problem_.AddResidualBlock(new UnaryCostFunction(2, 2, jacobian), NULL, x);
    }

    {
      double jacobian[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
      problem_.AddResidualBlock(new UnaryCostFunction(3, 3, jacobian), NULL, y);
    }

    {
      double jacobian = 5.0;
      problem_.AddResidualBlock(new UnaryCostFunction(1, 1, &jacobian), NULL, z);
    }

    {
      double jacobian1[] = { 0.0, 0.0, 0.0 };
      double jacobian2[] = { -5.0, -6.0 };
      problem_.AddResidualBlock(
          new BinaryCostFunction(1, 3, 2, jacobian1, jacobian2),
          NULL,
          y,
          x);
    }

    {
      double jacobian1[] = {2.0 };
      double jacobian2[] = { 3.0, -2.0 };
      problem_.AddResidualBlock(
          new BinaryCostFunction(1, 1, 2, jacobian1, jacobian2),
          NULL,
          z,
          x);
    }

    all_covariance_blocks_.push_back(make_pair(x, x));
    all_covariance_blocks_.push_back(make_pair(y, y));
    all_covariance_blocks_.push_back(make_pair(z, z));
    all_covariance_blocks_.push_back(make_pair(x, y));
    all_covariance_blocks_.push_back(make_pair(x, z));
    all_covariance_blocks_.push_back(make_pair(y, z));

    column_bounds_[x] = make_pair(0, 2);
    column_bounds_[y] = make_pair(2, 5);
    column_bounds_[z] = make_pair(5, 6);
  }
};

TEST_F(RankDeficientCovarianceTest, AutomaticTruncation) {
  // J
  //
  //   1  0  0  0  0  0
  //   0  1  0  0  0  0
  //   0  0  0  0  0  0
  //   0  0  0  0  0  0
  //   0  0  0  0  0  0
  //   0  0  0  0  0  5
  //  -5 -6  0  0  0  0
  //   3 -2  0  0  0  2

  // J'J
  //
  //  35 24  0  0  0  6
  //  24 41  0  0  0 -4
  //   0  0  0  0  0  0
  //   0  0  0  0  0  0
  //   0  0  0  0  0  0
  //   6 -4  0  0  0 29

  // pinv(J'J) computed using octave.
  double expected_covariance[] = {
     0.053998,  -0.033145,   0.000000,   0.000000,   0.000000,  -0.015744,
    -0.033145,   0.045067,   0.000000,   0.000000,   0.000000,   0.013074,
     0.000000,   0.000000,   0.000000,   0.000000,   0.000000,   0.000000,
     0.000000,   0.000000,   0.000000,   0.000000,   0.000000,   0.000000,
     0.000000,   0.000000,   0.000000,   0.000000,   0.000000,   0.000000,
    -0.015744,   0.013074,   0.000000,   0.000000,   0.000000,   0.039543
  };

  Covariance::Options options;
  options.algorithm_type = DENSE_SVD;
  options.null_space_rank = -1;
  ComputeAndCompareCovarianceBlocks(options, expected_covariance);
}

class LargeScaleCovarianceTest : public ::testing::Test {
 protected:
  virtual void SetUp() {
    num_parameter_blocks_ = 2000;
    parameter_block_size_ = 5;
    parameters_.reset(new double[parameter_block_size_ * num_parameter_blocks_]);

    Matrix jacobian(parameter_block_size_, parameter_block_size_);
    for (int i = 0; i < num_parameter_blocks_; ++i) {
      jacobian.setIdentity();
      jacobian *= (i + 1);

      double* block_i = parameters_.get() + i * parameter_block_size_;
      problem_.AddResidualBlock(new UnaryCostFunction(parameter_block_size_,
                                                      parameter_block_size_,
                                                      jacobian.data()),
                                NULL,
                                block_i );
      for (int j = i; j < num_parameter_blocks_; ++j) {
        double* block_j = parameters_.get() + j * parameter_block_size_;
        all_covariance_blocks_.push_back(make_pair(block_i, block_j));
      }
    }
  }

  void ComputeAndCompare(CovarianceAlgorithmType algorithm_type,
                         int num_threads) {
    Covariance::Options options;
    options.algorithm_type = algorithm_type;
    options.num_threads = num_threads;
    Covariance covariance(options);
    EXPECT_TRUE(covariance.Compute(all_covariance_blocks_, &problem_));

    Matrix expected(parameter_block_size_, parameter_block_size_);
    Matrix actual(parameter_block_size_, parameter_block_size_);
    const double kTolerance = 1e-16;

    for (int i = 0; i < num_parameter_blocks_; ++i) {
      expected.setIdentity();
      expected /= (i + 1.0) * (i + 1.0);

      double* block_i = parameters_.get() + i * parameter_block_size_;
      covariance.GetCovarianceBlock(block_i, block_i, actual.data());
      EXPECT_NEAR((expected - actual).norm(), 0.0, kTolerance)
          << "block: " << i << ", " << i << "\n"
          << "expected: \n" << expected << "\n"
          << "actual: \n" << actual;

      expected.setZero();
      for (int j = i + 1; j < num_parameter_blocks_; ++j) {
        double* block_j = parameters_.get() + j * parameter_block_size_;
        covariance.GetCovarianceBlock(block_i, block_j, actual.data());
        EXPECT_NEAR((expected - actual).norm(), 0.0, kTolerance)
            << "block: " << i << ", " << j << "\n"
            << "expected: \n" << expected << "\n"
            << "actual: \n" << actual;
      }
    }
  }

  scoped_array<double> parameters_;
  int parameter_block_size_;
  int num_parameter_blocks_;

  Problem problem_;
  vector<pair<const double*, const double*> > all_covariance_blocks_;
};

#if !defined(CERES_NO_SUITESPARSE) && defined(CERES_USE_OPENMP)

TEST_F(LargeScaleCovarianceTest, Parallel) {
  ComputeAndCompare(SPARSE_CHOLESKY, 4);
  ComputeAndCompare(SPARSE_QR, 4);
}

#endif  // !defined(CERES_NO_SUITESPARSE) && defined(CERES_USE_OPENMP)

}  // namespace internal
}  // namespace ceres
