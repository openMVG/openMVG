// Copyright (c) 2014 Matthieu Tourne

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "solver_affine_fund_kernel.hpp"
#include "solver_affine.hpp"
#include <iostream>

namespace openMVG {
namespace fund_affine {
namespace kernel {

// Gold standard algorithm (H.Z algo 14.1 p.351)
// minimize the sum of the distances (cost) of
// x = (x, y)^t         measured point
// x_^ = (x, y, z)^t    ideal point, estimated
// sum distance(x_i, x_^_i)^2 + distance(x_i', x_^_i')^2
void GSAffineSolver::Solve(const Mat2X &x1, const Mat2X &x2, std::vector<Mat3> *Fa) {
    assert(x1.rows() == 2);
    assert(x1.cols() >= 4);
    assert(x1.cols() == x2.cols);

    // point correspondances
    Mat4X X(4, x1.cols());

    // stack vertically
    X << x2,
         x1;

    // centroid
    Vec4 centroid = X.rowwise().mean();
    // DeltaX = X - centroid
    X.colwise() -= centroid;

    Eigen::JacobiSVD<Mat> svd(X.transpose(), Eigen::ComputeFullU|Eigen::ComputeFullV);
    Mat Vt = svd.matrixV();
    Vec4 v = Vt.col(3);

    Mat3 model_Fa(3,3);
    model_Fa <<  0,      0,         v[0],
                 0,      0,         v[1],
                v[2],  v[3], -v.transpose() * centroid;

    Fa->push_back(model_Fa);
}

// minimal configuration with 4 points (H.Z algo 14.2 p.352)
void MinimalAffineSolver::Solve(const Mat &x1, const Mat &x2, std::vector<Mat3> *Fa) {
    assert(x1.rows() == 2);
    assert(x1.cols() >= 4);
    assert(x1.rows() == x2.rows());
    assert(x1.cols() == x2.cols());

    // x2 = affine_H * x1
    Mat x1_minimal(x1.leftCols<3>());
    Mat x2_minimal(x2.leftCols<3>());

    // fourth point
    Vec3 x1_4 ( x1(0,3), x1(1,3), 1 );
    Vec3 x2_4 ( x2(0,3), x2(1,3), 1 );

    Mat3 affine_H;
    Mat3 model_Fa;

    if (Affine2DFromCorrespondencesLinear(x1_minimal, x2_minimal, &affine_H, 1e-8)) {

        // use 4th point to compute the epipolar line
        Vec3 l2 = CrossProductMatrix(affine_H * x1_4) * x2_4;

        // epipole for image2
        Vec3 e2(-l2[1], l2[0], 0);

        model_Fa = CrossProductMatrix(e2) * affine_H;

        // XX (mtourne) enforce shape constraint on the matrix here ?
        Fa->push_back(model_Fa);
    }
}

} // namespace kernel
} // namespace fund_affine
} // namespace openMVG
