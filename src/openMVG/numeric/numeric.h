
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

#ifndef OPENMVG_NUMERIC_NUMERIC_H
#define OPENMVG_NUMERIC_NUMERIC_H

//--
// Eigen
// http://eigen.tuxfamily.org/dox-devel/QuickRefPage.html
//--
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <Eigen/Geometry>
#include <Eigen/LU>
#include <Eigen/QR>
#include <Eigen/SparseCore>
#include <Eigen/SVD>

#include <cmath>
#include <numeric>
#include <string>
#include <iostream>
#include <vector>

namespace openMVG {


  using Eigen::Map;

  typedef Eigen::NumTraits<double> EigenDoubleTraits;

  typedef Eigen::Vector3d Vec3;
  typedef Eigen::Vector2i Vec2i;
  typedef Eigen::Vector2f Vec2f;
  typedef Eigen::Vector3f Vec3f;
  typedef Eigen::Matrix<double, 9, 1> Vec9;

  typedef Eigen::Quaternion<double> Quaternion;

  typedef Eigen::Matrix<double, 3, 3> Mat3;

#if defined(_WIN32) || defined(WIN32)
  // Handle alignment issue with Mat34, Vec2, Vec4, Vec6 on win32 with old compiler
  enum { NeedsToAlignMat34 = (sizeof(Eigen::Matrix<double, 3, 4>)%16)==0 };
  typedef Eigen::Matrix<double, 3, 4, ((NeedsToAlignMat34)==0 ? Eigen::Aligned : Eigen::DontAlign)> Mat34;

  enum { NeedsToAlignVec2= (sizeof(Eigen::Vector2d)%16)==0 };
  typedef Eigen::Matrix<double, 2, 1, ((NeedsToAlignVec2)==0 ? Eigen::Aligned : Eigen::DontAlign)> Vec2;

  enum { NeedsToAlignVec4= (sizeof(Eigen::Vector4d)%16)==0 };
  typedef Eigen::Matrix<double, 4, 1, ((NeedsToAlignVec4)==0 ? Eigen::Aligned : Eigen::DontAlign)> Vec4;

  enum { NeedsToAlignVec6= (sizeof(Eigen::Matrix<double, 6, 1>)%16)==0 };
  typedef Eigen::Matrix<double, 6, 1, ((NeedsToAlignVec6)==0 ? Eigen::Aligned : Eigen::DontAlign)> Vec6;
#else // defined(_WIN32) || defined(WIN32)
  typedef Eigen::Matrix<double, 3, 4> Mat34;
  typedef Eigen::Vector2d Vec2;
  typedef Eigen::Vector4d Vec4;
  typedef Eigen::Matrix<double, 6, 1> Vec6;
#endif // defined(_WIN32) || defined(WIN32)

  typedef Eigen::Matrix<double, 4, 4> Mat4;
  typedef Eigen::Matrix<unsigned int, Eigen::Dynamic, Eigen::Dynamic> Matu;

  typedef Eigen::Matrix<double, 3, 3, Eigen::RowMajor> RMat3;

  //-- General purpose Matrix and Vector
  typedef Eigen::MatrixXd Mat;
  typedef Eigen::VectorXd Vec;
  typedef Eigen::Matrix<unsigned int, Eigen::Dynamic, 1> Vecu;
  typedef Eigen::MatrixXf Matf;
  typedef Eigen::VectorXf Vecf;

  typedef Eigen::Matrix<double, 2, Eigen::Dynamic> Mat2X;
  typedef Eigen::Matrix<double, 3, Eigen::Dynamic> Mat3X;
  typedef Eigen::Matrix<double, 4, Eigen::Dynamic> Mat4X;

  typedef Eigen::Matrix<double, Eigen::Dynamic, 9> MatX9;

  //-- Sparse Matrix (Column major, and row major)
  typedef Eigen::SparseMatrix<double> sMat;
  typedef Eigen::SparseMatrix<double,Eigen::RowMajor> sRMat;

  //--------------
  //-- Function --
  //--------------

  /// Return the square of a number.
  template<typename T>
  inline T Square(T x) {
    return x * x;
  }

  /// Clamp return the number if inside range, else min or max range.
  template<typename T>
  inline T clamp(const T & val, const T& min, const T & max)  {
    return std::max(min, std::min(val, max));
    //(val < min) ? val : ((val>max) ? val : max);
  }

  Mat3 CrossProductMatrix(const Vec3 &x);

  Mat3 RotationAroundX(double angle);

  Mat3 RotationAroundY(double angle);

  Mat3 RotationAroundZ(double angle);

  // Degree to Radian (suppose input in [0;360])
  inline double D2R(double degree)  {
    return degree* M_PI / 180.0;
  }

  // Radian to degree
  inline double R2D(double radian)  {
    return radian / M_PI * 180.0;
  }

  /// Return in radian the mean rotation amplitude of the given rotation matrix
  /// Computed as the mean of matrix column dot products to an Identity matrix
  double  getRotationMagnitude(const Mat3 & R2);

  inline double SIGN(double x) {
    return x < 0.0 ? -1.0 : 1.0;
  }

  // L1 norm = Sum (|x0| + |x1| + |xn|)
  template<typename TVec>
  inline double NormL1(const TVec &x) {
    return x.array().abs().sum();
  }

  // L2 norm = Sqrt (Sum (x0^2 + x1^2 + xn^2))
  template<typename TVec>
  inline double NormL2(const TVec &x) {
    return x.norm();
  }

  // LInfinity norm = max (|x0|, |x1|, ..., |xn|)
  template<typename TVec>
  inline double NormLInfinity(const TVec &x) {
    return x.array().abs().maxCoeff();
  }

  template<typename TVec>
  inline double DistanceL1(const TVec &x, const TVec &y) {
    return (x - y).array().abs().sum();
  }

  template<typename TVec>
  inline double DistanceL2(const TVec &x, const TVec &y) {
    return (x - y).norm();
  }

  template<typename TVec>
  inline double DistanceLInfinity(const TVec &x, const TVec &y) {
    return NormLInfinity(x - y);
  }

  // Solve the linear system Ax = 0 via SVD. Store the solution in x, such that
  // ||x|| = 1.0. Return the singular value corresponding to the solution.
  // Destroys A and resizes x if necessary.
  template <typename TMat, typename TVec>
  double Nullspace(TMat *A, TVec *nullspace) {
    if (A->rows() >= A->cols()) {
      Eigen::JacobiSVD<TMat> svd(*A, Eigen::ComputeFullV);
      (*nullspace) = svd.matrixV().col(A->cols()-1);
      return svd.singularValues()(A->cols()-1);
    }
    // Extend A with rows of zeros to make it square. It's a hack, but is
    // necessary until Eigen supports SVD with more columns than rows.
    TMat A_extended(A->cols(), A->cols());
    A_extended.block(A->rows(), 0, A->cols() - A->rows(), A->cols()).setZero();
    A_extended.block(0,0, A->rows(), A->cols()) = (*A);
    return Nullspace(&A_extended, nullspace);
  }

  /// Solve the linear system Ax = 0 via SVD. Finds two solutions, x1 and x2, such
  /// that x1 is the best solution and x2 is the next best solution (in the L2
  /// norm sense). Store the solution in x1 and x2, such that ||x|| = 1.0. Return
  /// the singular value corresponding to the solution x1. Destroys A and resizes
  /// x if necessary.
  template <typename TMat, typename TVec1, typename TVec2>
  inline double Nullspace2(TMat *A, TVec1 *x1, TVec2 *x2) {
    if (A->rows() >= A->cols()) {
      Eigen::JacobiSVD<TMat> svd(*A,Eigen::ComputeFullV);
      TMat V = svd.matrixV();
      *x1 = V.col(A->cols() - 1);
      *x2 = V.col(A->cols() - 2);
      return svd.singularValues()(A->cols()-1);
    }
    // Extend A with rows of zeros to make it square. It's a hack, but is
    // necessary until Eigen supports SVD with more columns than rows.
    TMat A_extended(A->cols(), A->cols());
    A_extended.block(A->rows(), 0, A->cols() - A->rows(), A->cols()).setZero();
    A_extended.block(0,0, A->rows(), A->cols()) = (*A);
    return Nullspace2(&A_extended, x1, x2);
  }

  // Make a rotation matrix such that center becomes the direction of the
  // positive z-axis, and y is oriented close to up by default.
  Mat3 LookAt(const Vec3 &center, const Vec3 & up = Vec3::UnitY());

  Mat3 LookAt2(const Vec3 &eyePosition3D,
    const Vec3 &center3D = Vec3::Zero(),
    const Vec3 &upVector3D = Vec3::UnitY() );

#define SUM_OR_DYNAMIC(x,y) (x==Eigen::Dynamic||y==Eigen::Dynamic)?Eigen::Dynamic:(x+y)

  template<typename Derived1, typename Derived2>
  struct hstack_return {
    typedef typename Derived1::Scalar Scalar;
    enum {
      RowsAtCompileTime = Derived1::RowsAtCompileTime,
      ColsAtCompileTime = SUM_OR_DYNAMIC(Derived1::ColsAtCompileTime, Derived2::ColsAtCompileTime),
      Options = Derived1::Flags&Eigen::RowMajorBit ? Eigen::RowMajor : 0,
      MaxRowsAtCompileTime = Derived1::MaxRowsAtCompileTime,
      MaxColsAtCompileTime = SUM_OR_DYNAMIC(Derived1::MaxColsAtCompileTime, Derived2::MaxColsAtCompileTime)
    };
    typedef Eigen::Matrix<Scalar,
      RowsAtCompileTime,
      ColsAtCompileTime,
      Options,
      MaxRowsAtCompileTime,
      MaxColsAtCompileTime> type;
  };

  template<typename Derived1, typename Derived2>
  typename hstack_return<Derived1,Derived2>::type
    HStack (const Eigen::MatrixBase<Derived1>& lhs, const Eigen::MatrixBase<Derived2>& rhs) {
      typename hstack_return<Derived1,Derived2>::type res;
      res.resize(lhs.rows(), lhs.cols()+rhs.cols());
      res << lhs, rhs;
      return res;
  };


  template<typename Derived1, typename Derived2>
  struct vstack_return {
    typedef typename Derived1::Scalar Scalar;
    enum {
      RowsAtCompileTime = SUM_OR_DYNAMIC(Derived1::RowsAtCompileTime, Derived2::RowsAtCompileTime),
      ColsAtCompileTime = Derived1::ColsAtCompileTime,
      Options = Derived1::Flags&Eigen::RowMajorBit ? Eigen::RowMajor : 0,
      MaxRowsAtCompileTime = SUM_OR_DYNAMIC(Derived1::MaxRowsAtCompileTime, Derived2::MaxRowsAtCompileTime),
      MaxColsAtCompileTime = Derived1::MaxColsAtCompileTime
    };
    typedef Eigen::Matrix<Scalar,
      RowsAtCompileTime,
      ColsAtCompileTime,
      Options,
      MaxRowsAtCompileTime,
      MaxColsAtCompileTime> type;
  };

  template<typename Derived1, typename Derived2>
  typename vstack_return<Derived1,Derived2>::type
    VStack (const Eigen::MatrixBase<Derived1>& lhs, const Eigen::MatrixBase<Derived2>& rhs) {
      typename vstack_return<Derived1,Derived2>::type res;
      res.resize(lhs.rows()+rhs.rows(), lhs.cols());
      res << lhs, rhs;
      return res;
  };
#undef SUM_OR_DYNAMIC

  template<typename TMat>
  inline double FrobeniusNorm(const TMat &A) {
    return sqrt(A.array().abs2().sum());
  }

  template<typename TMat>
  inline double FrobeniusDistance(const TMat &A, const TMat &B) {
    return FrobeniusNorm(A - B);
  }

  template<class TMat>
  double CosinusBetweenMatrices(const TMat &a, const TMat &b) {
    return (a.array() * b.array()).sum() /
      FrobeniusNorm(a) / FrobeniusNorm(b);
  }

  template <typename TMat, typename TCols>
  TMat ExtractColumns(const TMat &A, const TCols &columns) {
    TMat compressed(A.rows(), columns.size());
    for (size_t i = 0; i < static_cast<size_t>(columns.size()); ++i) {
      compressed.col(i) = A.col(columns[i]);
    }
    return compressed;
  }

  void MeanAndVarianceAlongRows(const Mat &A,
    Vec *mean_pointer,
    Vec *variance_pointer);

  bool exportMatToTextFile(const Mat & mat, const std::string & filename,
    const std::string & sPrefix = "A");

  inline int is_finite(const double val)
  {
#ifdef _WIN32
    return _finite(val);
#else
    return std::isfinite(val);
#endif
  }

  /// Get back the min, mean, median and the max
  ///  values of an iterable sequence.
  template <typename Type, typename DataInputIterator>
  void minMaxMeanMedian(DataInputIterator begin, DataInputIterator end,
        Type & min, Type & max, Type & mean, Type & median)
  {
    if(std::distance(begin,end)<1)
      return;

    std::vector<Type> vec_val(begin, end);
    std::sort(vec_val.begin(), vec_val.end());
    min = vec_val[0];
    max = vec_val[vec_val.size()-1];
    mean = accumulate(vec_val.begin(), vec_val.end(), 0.0)
        / static_cast<double>(vec_val.size());
    median = vec_val[vec_val.size()/2];
  }

  /// Display to the console the min, mean, median and the max
  ///  values of an iterable sequence.
  template <typename Type, typename DataInputIterator>
  void minMaxMeanMedian(DataInputIterator begin, DataInputIterator end)
  {
    Type min, max, mean, median;
    minMaxMeanMedian(begin,end, min, max, mean, median);
    std::cout << "\n"
      << "\t min: " << min << "\n"
      << "\t mean: " << mean << "\n"
      << "\t median: " << median << std::endl
      << "\t max: " << max << std::endl;
  }

} // namespace openMVG


#endif  // OPENMVG_NUMERIC_NUMERIC_H
