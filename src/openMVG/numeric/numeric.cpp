
// Copyright (c) 2007, 2008 libmv authors.
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to
// deal in the Software without restriction, including without limitation the
// rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
// sell copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
// IN THE SOFTWARE.

// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/numeric/numeric.h"
#include <iostream>
#include <fstream>
#include <string>

namespace openMVG {

Mat3 CrossProductMatrix(const Vec3 &x) {
  Mat3 X;
  X << 0, -x(2),  x(1),
    x(2),     0, -x(0),
   -x(1),  x(0),     0;
  return X;
}

Mat3 RotationAroundX(double angle) {
  return Eigen::AngleAxisd(angle, Vec3::UnitX()).toRotationMatrix();
}

Mat3 RotationAroundY(double angle) {
  return Eigen::AngleAxisd(angle, Vec3::UnitY()).toRotationMatrix();
}

Mat3 RotationAroundZ(double angle) {
  return Eigen::AngleAxisd(angle, Vec3::UnitZ()).toRotationMatrix();
}

double getRotationMagnitude(const Mat3 & R2) {
  const Mat3 R1 = Mat3::Identity();
  double cos_theta = (R1.array() * R2.array()).sum() / 3.0;
  cos_theta = clamp(cos_theta, -1.0, 1.0);
  return std::acos(cos_theta);
}

Mat3 LookAt(const Vec3 &center, const Vec3 & up) {
  Vec3 zc = center.normalized();
  Vec3 xc = up.cross(zc).normalized();
  Vec3 yc = zc.cross(xc);
  Mat3 R;
  R.row(0) = xc;
  R.row(1) = yc;
  R.row(2) = zc;
  return R;
}

//eyePosition3D is a XYZ position. This is where you are (your eye is).
//center3D is the XYZ position where you want to look at.
//upVector3D is a XYZ normalized vector. Quite often 0.0, 1.0, 0.0

Mat3 LookAt2(const Vec3 &eyePosition3D,
  const Vec3 &center3D,
  const Vec3 &upVector3D )
{
  Vec3 forward, side, up;
  Mat3 matrix2, resultMatrix;
  //------------------
  forward = center3D - eyePosition3D;
  forward.normalize();
  //------------------
  //Side = forward x up
  //ComputeNormalOfPlane(side, forward, upVector3D);
  side[0]=(forward[1]*upVector3D[2])-(forward[2]*upVector3D[1]);
  side[1]=(forward[2]*upVector3D[0])-(forward[0]*upVector3D[2]);
  side[2]=(forward[0]*upVector3D[1])-(forward[1]*upVector3D[0]);
  side.normalize();
  //------------------
  //Recompute up as: up = side x forward
  //ComputeNormalOfPlane(up, side, forward);
  up[0]=(side[1]*forward[2])-(side[2]*forward[1]);
  up[1]=(side[2]*forward[0])-(side[0]*forward[2]);
  up[2]=(side[0]*forward[1])-(side[1]*forward[0]);

  //------------------
  matrix2(0) = side[0];
  matrix2(1) = side[1];
  matrix2(2) = side[2];
  //------------------
  matrix2(3) = up[0];
  matrix2(4) = up[1];
  matrix2(5) = up[2];
  //------------------
  matrix2(6) = -forward[0];
  matrix2(7) = -forward[1];
  matrix2(8) = -forward[2];

  return matrix2;
}

void MeanAndVarianceAlongRows(const Mat &A,
  Vec *mean_pointer,
  Vec *variance_pointer) {
    Vec &mean = *mean_pointer;
    Vec &variance = *variance_pointer;
    Mat::Index n = A.rows();
    double m = static_cast<double>(A.cols());
    mean = variance = Vec::Zero(n);

    for (Mat::Index i = 0; i < n; ++i) {
      mean(i) += A.row(i).array().sum();
      variance(i) += (A.row(i).array() * A.row(i).array()).array().sum();
    }

    mean /= m;
    for (Mat::Index i = 0; i < n; ++i) {
      variance(i) = variance(i) / m - Square(mean(i));
    }
}

bool exportMatToTextFile(const Mat & mat, const std::string & filename,
  const std::string & sPrefix)
{
  bool bOk = false;
  std::ofstream outfile;
  outfile.open(filename.c_str(), std::ios_base::out);
  if (outfile.is_open()) {
    outfile << sPrefix << "=[" << std::endl;
    for (int j=0; j < mat.rows(); ++j)  {
      for (int i=0; i < mat.cols(); ++i)  {
        outfile << mat(j,i) << " ";
      }
      outfile << ";\n";
    }
    outfile << "];";
    bOk = true;
  }
  outfile.close();
  return bOk;
}

}  // namespace openMVG

