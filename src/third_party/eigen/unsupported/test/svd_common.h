// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2008 Gael Guennebaud <gael.guennebaud@inria.fr>
// Copyright (C) 2009 Benoit Jacob <jacob.benoit.1@gmail.com>
//
// Copyright (C) 2013 Gauthier Brun <brun.gauthier@gmail.com>
// Copyright (C) 2013 Nicolas Carre <nicolas.carre@ensimag.fr>
// Copyright (C) 2013 Jean Ceccato <jean.ceccato@ensimag.fr>
// Copyright (C) 2013 Pierre Zoppitelli <pierre.zoppitelli@ensimag.fr>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

// discard stack allocation as that too bypasses malloc
#define EIGEN_STACK_ALLOCATION_LIMIT 0
#define EIGEN_RUNTIME_NO_MALLOC

#include "main.h"
#include <unsupported/Eigen/SVD>
#include <Eigen/LU>


// check if "svd" is the good image of "m"  
template<typename MatrixType, typename SVD>
void svd_check_full(const MatrixType& m, const SVD& svd)
{
  typedef typename MatrixType::Index Index;
  Index rows = m.rows();
  Index cols = m.cols();
  enum {
    RowsAtCompileTime = MatrixType::RowsAtCompileTime,
    ColsAtCompileTime = MatrixType::ColsAtCompileTime
  };

  typedef typename MatrixType::Scalar Scalar;
  typedef Matrix<Scalar, RowsAtCompileTime, RowsAtCompileTime> MatrixUType;
  typedef Matrix<Scalar, ColsAtCompileTime, ColsAtCompileTime> MatrixVType;

  
  MatrixType sigma = MatrixType::Zero(rows, cols);
  sigma.diagonal() = svd.singularValues().template cast<Scalar>();
  MatrixUType u = svd.matrixU();
  MatrixVType v = svd.matrixV();
  VERIFY_IS_APPROX(m, u * sigma * v.adjoint());
  VERIFY_IS_UNITARY(u);
  VERIFY_IS_UNITARY(v);
} // end svd_check_full



// Compare to a reference value
template<typename MatrixType, typename SVD>
void svd_compare_to_full(const MatrixType& m,
			 unsigned int computationOptions,
			 const SVD& referenceSvd)
{
  typedef typename MatrixType::Index Index;
  Index rows = m.rows();
  Index cols = m.cols();
  Index diagSize = (std::min)(rows, cols);

  SVD svd(m, computationOptions);

  VERIFY_IS_APPROX(svd.singularValues(), referenceSvd.singularValues());
  if(computationOptions & ComputeFullU)
    VERIFY_IS_APPROX(svd.matrixU(), referenceSvd.matrixU());
  if(computationOptions & ComputeThinU)
    VERIFY_IS_APPROX(svd.matrixU(), referenceSvd.matrixU().leftCols(diagSize));
  if(computationOptions & ComputeFullV)
    VERIFY_IS_APPROX(svd.matrixV(), referenceSvd.matrixV());
  if(computationOptions & ComputeThinV)
    VERIFY_IS_APPROX(svd.matrixV(), referenceSvd.matrixV().leftCols(diagSize));
} // end svd_compare_to_full



template<typename MatrixType, typename SVD>
void svd_solve(const MatrixType& m, unsigned int computationOptions)
{
  typedef typename MatrixType::Scalar Scalar;
  typedef typename MatrixType::Index Index;
  Index rows = m.rows();
  Index cols = m.cols();

  enum {
    RowsAtCompileTime = MatrixType::RowsAtCompileTime,
    ColsAtCompileTime = MatrixType::ColsAtCompileTime
  };

  typedef Matrix<Scalar, RowsAtCompileTime, Dynamic> RhsType;
  typedef Matrix<Scalar, ColsAtCompileTime, Dynamic> SolutionType;

  RhsType rhs = RhsType::Random(rows, internal::random<Index>(1, cols));
  SVD svd(m, computationOptions);
  SolutionType x = svd.solve(rhs);
  // evaluate normal equation which works also for least-squares solutions
  VERIFY_IS_APPROX(m.adjoint()*m*x,m.adjoint()*rhs);
} // end svd_solve


// test computations options
// 2 functions because Jacobisvd can return before the second function
template<typename MatrixType, typename SVD>
void svd_test_computation_options_1(const MatrixType& m, const SVD& fullSvd)
{
  svd_check_full< MatrixType, SVD >(m, fullSvd);
  svd_solve< MatrixType, SVD >(m, ComputeFullU | ComputeFullV);
}


template<typename MatrixType, typename SVD>
void svd_test_computation_options_2(const MatrixType& m, const SVD& fullSvd)
{
  svd_compare_to_full< MatrixType, SVD >(m, ComputeFullU, fullSvd);
  svd_compare_to_full< MatrixType, SVD >(m, ComputeFullV, fullSvd);
  svd_compare_to_full< MatrixType, SVD >(m, 0, fullSvd);

  if (MatrixType::ColsAtCompileTime == Dynamic) {
    // thin U/V are only available with dynamic number of columns
 
    svd_compare_to_full< MatrixType, SVD >(m, ComputeFullU|ComputeThinV, fullSvd);
    svd_compare_to_full< MatrixType, SVD >(m,              ComputeThinV, fullSvd);
    svd_compare_to_full< MatrixType, SVD >(m, ComputeThinU|ComputeFullV, fullSvd);
    svd_compare_to_full< MatrixType, SVD >(m, ComputeThinU             , fullSvd);
    svd_compare_to_full< MatrixType, SVD >(m, ComputeThinU|ComputeThinV, fullSvd);
    svd_solve<MatrixType, SVD>(m, ComputeFullU | ComputeThinV);
    svd_solve<MatrixType, SVD>(m, ComputeThinU | ComputeFullV);
    svd_solve<MatrixType, SVD>(m, ComputeThinU | ComputeThinV);
    
    typedef typename MatrixType::Index Index;
    Index diagSize = (std::min)(m.rows(), m.cols());
    SVD svd(m, ComputeThinU | ComputeThinV);
    VERIFY_IS_APPROX(m, svd.matrixU().leftCols(diagSize) * svd.singularValues().asDiagonal() * svd.matrixV().leftCols(diagSize).adjoint());
  }
}

template<typename MatrixType, typename SVD> 
void svd_verify_assert(const MatrixType& m)
{
  typedef typename MatrixType::Scalar Scalar;
  typedef typename MatrixType::Index Index;
  Index rows = m.rows();
  Index cols = m.cols();

  enum {
    RowsAtCompileTime = MatrixType::RowsAtCompileTime,
    ColsAtCompileTime = MatrixType::ColsAtCompileTime
  };

  typedef Matrix<Scalar, RowsAtCompileTime, 1> RhsType;
  RhsType rhs(rows);
  SVD svd;
  VERIFY_RAISES_ASSERT(svd.matrixU())
  VERIFY_RAISES_ASSERT(svd.singularValues())
  VERIFY_RAISES_ASSERT(svd.matrixV())
  VERIFY_RAISES_ASSERT(svd.solve(rhs))
  MatrixType a = MatrixType::Zero(rows, cols);
  a.setZero();
  svd.compute(a, 0);
  VERIFY_RAISES_ASSERT(svd.matrixU())
  VERIFY_RAISES_ASSERT(svd.matrixV())
  svd.singularValues();
  VERIFY_RAISES_ASSERT(svd.solve(rhs))
    
  if (ColsAtCompileTime == Dynamic)
  {
    svd.compute(a, ComputeThinU);
    svd.matrixU();
    VERIFY_RAISES_ASSERT(svd.matrixV())
    VERIFY_RAISES_ASSERT(svd.solve(rhs))
    svd.compute(a, ComputeThinV);
    svd.matrixV();
    VERIFY_RAISES_ASSERT(svd.matrixU())
    VERIFY_RAISES_ASSERT(svd.solve(rhs))
  }
  else
  {
    VERIFY_RAISES_ASSERT(svd.compute(a, ComputeThinU))
    VERIFY_RAISES_ASSERT(svd.compute(a, ComputeThinV))
  }
}

// work around stupid msvc error when constructing at compile time an expression that involves
// a division by zero, even if the numeric type has floating point
template<typename Scalar>
EIGEN_DONT_INLINE Scalar zero() { return Scalar(0); }

// workaround aggressive optimization in ICC
template<typename T> EIGEN_DONT_INLINE  T sub(T a, T b) { return a - b; }


template<typename MatrixType, typename SVD>
void svd_inf_nan()
{
  // all this function does is verify we don't iterate infinitely on nan/inf values

  SVD svd;
  typedef typename MatrixType::Scalar Scalar;
  Scalar some_inf = Scalar(1) / zero<Scalar>();
  VERIFY(sub(some_inf, some_inf) != sub(some_inf, some_inf));
  svd.compute(MatrixType::Constant(10,10,some_inf), ComputeFullU | ComputeFullV);

  Scalar some_nan = zero<Scalar> () / zero<Scalar> ();
  VERIFY(some_nan != some_nan);
  svd.compute(MatrixType::Constant(10,10,some_nan), ComputeFullU | ComputeFullV);

  MatrixType m = MatrixType::Zero(10,10);
  m(internal::random<int>(0,9), internal::random<int>(0,9)) = some_inf;
  svd.compute(m, ComputeFullU | ComputeFullV);

  m = MatrixType::Zero(10,10);
  m(internal::random<int>(0,9), internal::random<int>(0,9)) = some_nan;
  svd.compute(m, ComputeFullU | ComputeFullV);
}


template<typename SVD>
void svd_preallocate()
{
  Vector3f v(3.f, 2.f, 1.f);
  MatrixXf m = v.asDiagonal();

  internal::set_is_malloc_allowed(false);
  VERIFY_RAISES_ASSERT(VectorXf v(10);)
    SVD svd;
  internal::set_is_malloc_allowed(true);
  svd.compute(m);
  VERIFY_IS_APPROX(svd.singularValues(), v);

  SVD svd2(3,3);
  internal::set_is_malloc_allowed(false);
  svd2.compute(m);
  internal::set_is_malloc_allowed(true);
  VERIFY_IS_APPROX(svd2.singularValues(), v);
  VERIFY_RAISES_ASSERT(svd2.matrixU());
  VERIFY_RAISES_ASSERT(svd2.matrixV());
  svd2.compute(m, ComputeFullU | ComputeFullV);
  VERIFY_IS_APPROX(svd2.matrixU(), Matrix3f::Identity());
  VERIFY_IS_APPROX(svd2.matrixV(), Matrix3f::Identity());
  internal::set_is_malloc_allowed(false);
  svd2.compute(m);
  internal::set_is_malloc_allowed(true);

  SVD svd3(3,3,ComputeFullU|ComputeFullV);
  internal::set_is_malloc_allowed(false);
  svd2.compute(m);
  internal::set_is_malloc_allowed(true);
  VERIFY_IS_APPROX(svd2.singularValues(), v);
  VERIFY_IS_APPROX(svd2.matrixU(), Matrix3f::Identity());
  VERIFY_IS_APPROX(svd2.matrixV(), Matrix3f::Identity());
  internal::set_is_malloc_allowed(false);
  svd2.compute(m, ComputeFullU|ComputeFullV);
  internal::set_is_malloc_allowed(true);
}





