// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2008 Gael Guennebaud <gael.guennebaud@inria.fr>
// Copyright (C) 2009 Benoit Jacob <jacob.benoit.1@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "svd_common.h"

template<typename MatrixType, int QRPreconditioner>
void jacobisvd_check_full(const MatrixType& m, const JacobiSVD<MatrixType, QRPreconditioner>& svd)
{
  svd_check_full<MatrixType, JacobiSVD<MatrixType, QRPreconditioner > >(m, svd);
}

template<typename MatrixType, int QRPreconditioner>
void jacobisvd_compare_to_full(const MatrixType& m,
                               unsigned int computationOptions,
                               const JacobiSVD<MatrixType, QRPreconditioner>& referenceSvd)
{
  svd_compare_to_full<MatrixType, JacobiSVD<MatrixType, QRPreconditioner> >(m, computationOptions, referenceSvd);
}


template<typename MatrixType, int QRPreconditioner>
void jacobisvd_solve(const MatrixType& m, unsigned int computationOptions)
{
  svd_solve< MatrixType, JacobiSVD< MatrixType, QRPreconditioner > >(m, computationOptions);
}



template<typename MatrixType, int QRPreconditioner>
void jacobisvd_test_all_computation_options(const MatrixType& m)
{
  
  if (QRPreconditioner == NoQRPreconditioner && m.rows() != m.cols())
    return;

  JacobiSVD< MatrixType, QRPreconditioner > fullSvd(m, ComputeFullU|ComputeFullV);
  svd_test_computation_options_1< MatrixType, JacobiSVD< MatrixType, QRPreconditioner > >(m, fullSvd);

  if(QRPreconditioner == FullPivHouseholderQRPreconditioner)
    return;
  svd_test_computation_options_2< MatrixType, JacobiSVD< MatrixType, QRPreconditioner > >(m, fullSvd);

}

template<typename MatrixType>
void jacobisvd(const MatrixType& a = MatrixType(), bool pickrandom = true)
{
  MatrixType m = pickrandom ? MatrixType::Random(a.rows(), a.cols()) : a;

  jacobisvd_test_all_computation_options<MatrixType, FullPivHouseholderQRPreconditioner>(m);
  jacobisvd_test_all_computation_options<MatrixType, ColPivHouseholderQRPreconditioner>(m);
  jacobisvd_test_all_computation_options<MatrixType, HouseholderQRPreconditioner>(m);
  jacobisvd_test_all_computation_options<MatrixType, NoQRPreconditioner>(m);
}


template<typename MatrixType> 
void jacobisvd_verify_assert(const MatrixType& m)
{
  
  svd_verify_assert<MatrixType, JacobiSVD< MatrixType > >(m);

  typedef typename MatrixType::Index Index;
  Index rows = m.rows();
  Index cols = m.cols();

  enum {
    RowsAtCompileTime = MatrixType::RowsAtCompileTime,
    ColsAtCompileTime = MatrixType::ColsAtCompileTime
  };

  MatrixType a = MatrixType::Zero(rows, cols);
  a.setZero();

  if (ColsAtCompileTime == Dynamic)
  {
    JacobiSVD<MatrixType, FullPivHouseholderQRPreconditioner> svd_fullqr;
    VERIFY_RAISES_ASSERT(svd_fullqr.compute(a, ComputeFullU|ComputeThinV))
    VERIFY_RAISES_ASSERT(svd_fullqr.compute(a, ComputeThinU|ComputeThinV))
    VERIFY_RAISES_ASSERT(svd_fullqr.compute(a, ComputeThinU|ComputeFullV))
  }
}

template<typename MatrixType>
void jacobisvd_method()
{
  enum { Size = MatrixType::RowsAtCompileTime };
  typedef typename MatrixType::RealScalar RealScalar;
  typedef Matrix<RealScalar, Size, 1> RealVecType;
  MatrixType m = MatrixType::Identity();
  VERIFY_IS_APPROX(m.jacobiSvd().singularValues(), RealVecType::Ones());
  VERIFY_RAISES_ASSERT(m.jacobiSvd().matrixU());
  VERIFY_RAISES_ASSERT(m.jacobiSvd().matrixV());
  VERIFY_IS_APPROX(m.jacobiSvd(ComputeFullU|ComputeFullV).solve(m), m);
}



template<typename MatrixType>
void jacobisvd_inf_nan()
{
  svd_inf_nan<MatrixType, JacobiSVD< MatrixType > >();
}


// Regression test for bug 286: JacobiSVD loops indefinitely with some
// matrices containing denormal numbers.
void jacobisvd_bug286()
{
#if defined __INTEL_COMPILER
// shut up warning #239: floating point underflow
#pragma warning push
#pragma warning disable 239
#endif
  Matrix2d M;
  M << -7.90884e-313, -4.94e-324,
                 0, 5.60844e-313;
#if defined __INTEL_COMPILER
#pragma warning pop
#endif
  JacobiSVD<Matrix2d> svd;
  svd.compute(M); // just check we don't loop indefinitely
}


void jacobisvd_preallocate()
{
  svd_preallocate< JacobiSVD <MatrixXf> >();
}

void test_jacobisvd()
{
  CALL_SUBTEST_11(( jacobisvd<Matrix<double,Dynamic,Dynamic> >
		    (Matrix<double,Dynamic,Dynamic>(16, 6)) ));

  CALL_SUBTEST_3(( jacobisvd_verify_assert(Matrix3f()) ));
  CALL_SUBTEST_4(( jacobisvd_verify_assert(Matrix4d()) ));
  CALL_SUBTEST_7(( jacobisvd_verify_assert(MatrixXf(10,12)) ));
  CALL_SUBTEST_8(( jacobisvd_verify_assert(MatrixXcd(7,5)) ));

  for(int i = 0; i < g_repeat; i++) {
    Matrix2cd m;
    m << 0, 1,
         0, 1;
    CALL_SUBTEST_1(( jacobisvd(m, false) ));
    m << 1, 0,
         1, 0;
    CALL_SUBTEST_1(( jacobisvd(m, false) ));

    Matrix2d n;
    n << 0, 0,
         0, 0;
    CALL_SUBTEST_2(( jacobisvd(n, false) ));
    n << 0, 0,
         0, 1;
    CALL_SUBTEST_2(( jacobisvd(n, false) ));
    
    CALL_SUBTEST_3(( jacobisvd<Matrix3f>() ));
    CALL_SUBTEST_4(( jacobisvd<Matrix4d>() ));
    CALL_SUBTEST_5(( jacobisvd<Matrix<float,3,5> >() ));
    CALL_SUBTEST_6(( jacobisvd<Matrix<double,Dynamic,2> >(Matrix<double,Dynamic,2>(10,2)) ));

    int r = internal::random<int>(1, 30),
        c = internal::random<int>(1, 30);
    CALL_SUBTEST_7(( jacobisvd<MatrixXf>(MatrixXf(r,c)) ));
    CALL_SUBTEST_8(( jacobisvd<MatrixXcd>(MatrixXcd(r,c)) ));
    (void) r;
    (void) c;

    // Test on inf/nan matrix
    CALL_SUBTEST_7( jacobisvd_inf_nan<MatrixXf>() );
  }

  CALL_SUBTEST_7(( jacobisvd<MatrixXf>(MatrixXf(internal::random<int>(EIGEN_TEST_MAX_SIZE/4, EIGEN_TEST_MAX_SIZE/2), internal::random<int>(EIGEN_TEST_MAX_SIZE/4, EIGEN_TEST_MAX_SIZE/2))) ));
  CALL_SUBTEST_8(( jacobisvd<MatrixXcd>(MatrixXcd(internal::random<int>(EIGEN_TEST_MAX_SIZE/4, EIGEN_TEST_MAX_SIZE/3), internal::random<int>(EIGEN_TEST_MAX_SIZE/4, EIGEN_TEST_MAX_SIZE/3))) ));


  // test matrixbase method
  CALL_SUBTEST_1(( jacobisvd_method<Matrix2cd>() ));
  CALL_SUBTEST_3(( jacobisvd_method<Matrix3f>() ));


  // Test problem size constructors
  CALL_SUBTEST_7( JacobiSVD<MatrixXf>(10,10) );

  // Check that preallocation avoids subsequent mallocs
  CALL_SUBTEST_9( jacobisvd_preallocate() );

  // Regression check for bug 286
  CALL_SUBTEST_2( jacobisvd_bug286() );
}
