// Copyright (c) 2014 Matthieu Tourne

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "testing/testing.h"
#include "openMVG/numeric/numeric.h"
#include "openMVG/multiview/projection.hpp"
#include "openMVG/robust_estimation/score_evaluator.hpp"

#include "solver_affine_fund_kernel.hpp"

using namespace openMVG;

struct TestData {
  //-- Dataset that encapsulate :
  // 3D points and their projection given P1 and P2
  // Link between the two camera [R|t]
  Mat3X X;
  Mat3 R;
  Vec3 t;
  Mat34 P1, P2;
  Mat2X x1, x2;
};

TestData SomeTestData() {
  TestData d;

  // --
  // Modeling random 3D points,
  // Consider first camera as [R=I|t=0],
  // Second camera as [R=Rx*Ry*Rz|t=random],
  // Compute projection of the 3D points onto image plane.
  // --
  d.X = Mat3X::Random(3,4);

  //-- Make point in front to the cameras.
  d.X.row(0).array() -= .5;
  d.X.row(1).array() -= .5;
  d.X.row(2).array() += 1.0;

  d.R = RotationAroundZ(0.8) * RotationAroundX(0.5) * RotationAroundY(0.7);
  d.t.setRandom();
  //d.t(1)  = 10;

  P_From_KRt(Mat3::Identity(), Mat3::Identity(), Vec3::Zero(), &d.P1);
  P_From_KRt(Mat3::Identity(), d.R, d.t, &d.P2);

  Project(d.P1, d.X, &d.x1);
  Project(d.P2, d.X, &d.x2);

  return d;
}

TEST(gs_vs_minimal, easy) {
    TestData d = SomeTestData();

    std::vector<Mat3> Fa;
    Vec3 e1;
    Vec3 e2;

    // use gold standard
    fund_affine::kernel::GSAffineSolver gs_solver;
    fund_affine::kernel::AffineSampsonError cost_function;

    gs_solver.Solve(d.x1, d.x2, &Fa);
    double cost = cost_function.Error(Fa[0], d.x1, d.x2);

    // Fa's are only defined up to a scale, normalize by e
    Fa[0] /= Fa[0](2,2);
    std::cout << "GS Fundamental Affine: " << std::endl << Fa[0] << std::endl;
    std::cout << "cost : " << std::endl << cost << std::endl;

    // use the minimal configuration
    fund_affine::kernel::MinimalAffineSolver min_solver;
    min_solver.Solve(d.x1, d.x2, &Fa);

    cost = cost_function.Error(Fa[1], d.x1, d.x2);

    Fa[1] /= Fa[1](2,2);
    std::cout << "Minimal Configuration Fundamental Affine: "
              << std::endl << Fa[1] << std::endl;
    std::cout << "cost : " << std::endl << cost << std::endl;

    EXPECT_NEAR(0, (Fa[0] - Fa[1]).squaredNorm(), 1e-20);
}

template<typename TMat>
bool ExpectFundamentalAffineProperties(const TMat &F,
                                       const Mat &ptsA,
                                       const Mat &ptsB,
                                       double precision) {
  bool bOk = true;
  // test the shape
  bOk &= (F(0,0) == 0 && F(0,1) == 0 && F(0,2) != 0
          && F(1,0) == 0 && F(1,1) == 0 && F(1,2) != 0
          && F(2,0) != 0 && F(2,1) != 0 && F(2,2) != 0);
  Mat hptsA, hptsB;
  EuclideanToHomogeneous(ptsA, &hptsA);
  EuclideanToHomogeneous(ptsB, &hptsB);
  for (int i = 0; i < ptsA.cols(); ++i) {
      double residual = hptsB.col(i).dot(F * hptsA.col(i));
      bOk &= residual < precision;
  }
  return bOk;
}

template <class Kernel>
bool ExpectKernelProperties(const Mat &x1,
                            const Mat &x2) {
    bool bOk = true;
    Kernel kernel(x1, x2);
    std::vector<size_t> samples;
    for (size_t i = 0; i < x1.cols(); ++i) {
        samples.push_back(i);
    }

    std::vector<Mat3> Fs;
    kernel.Fit(samples, &Fs);

    robust::ScorerEvaluator<Kernel> scorer(0.3);
    std::vector<size_t> inliers;

    bOk &= (!Fs.empty());
    for (int i = 0; i < Fs.size(); ++i) {
        bOk &= ExpectFundamentalAffineProperties(Fs[i], x1, x2, 1e-8);
        double cost = scorer.Score(kernel, Fs[i], samples, &inliers);
        std::cout << "cost: " << cost;
    }
    return bOk;
}

TEST(solving_minimal, easy) {
    TestData d = SomeTestData();

    typedef fund_affine::kernel::MinimalAffineKernel Kernel;
    EXPECT_TRUE(ExpectKernelProperties<Kernel>(d.x1, d.x2));
}


/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
