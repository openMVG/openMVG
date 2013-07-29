// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2013 Gauthier Brun <brun.gauthier@gmail.com>
// Copyright (C) 2013 Nicolas Carre <nicolas.carre@ensimag.fr>
// Copyright (C) 2013 Jean Ceccato <jean.ceccato@ensimag.fr>
// Copyright (C) 2013 Pierre Zoppitelli <pierre.zoppitelli@ensimag.fr>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/

#include "svd_common.h"
#include <iostream>
#include <Eigen/LU>

// check if "svd" is the good image of "m"  
template<typename MatrixType>
void bdcsvd_check_full(const MatrixType& m, const BDCSVD<MatrixType>& svd)
{
  svd_check_full< MatrixType, BDCSVD< MatrixType > >(m, svd);
}

// Compare to a reference value
template<typename MatrixType>
void bdcsvd_compare_to_full(const MatrixType& m,
			    unsigned int computationOptions,
			    const BDCSVD<MatrixType>& referenceSvd)
{
  svd_compare_to_full< MatrixType, BDCSVD< MatrixType > >(m, computationOptions, referenceSvd);
} // end bdcsvd_compare_to_full


template<typename MatrixType>
void bdcsvd_solve(const MatrixType& m, unsigned int computationOptions)
{
  svd_solve< MatrixType, BDCSVD< MatrixType > >(m, computationOptions);
} //  end template bdcsvd_solve


// test the computations options
template<typename MatrixType>
void bdcsvd_test_all_computation_options(const MatrixType& m)
{
  BDCSVD<MatrixType> fullSvd(m, ComputeFullU|ComputeFullV);
  svd_test_computation_options_1< MatrixType, BDCSVD< MatrixType > >(m, fullSvd); 
  svd_test_computation_options_2< MatrixType, BDCSVD< MatrixType > >(m, fullSvd); 
} // end bdcsvd_test_all_computation_options


// Call a test with all the computations options
template<typename MatrixType>
void bdcsvd(const MatrixType& a = MatrixType(), bool pickrandom = true)
{
  MatrixType m = pickrandom ? MatrixType::Random(a.rows(), a.cols()) : a;
  bdcsvd_test_all_computation_options<MatrixType>(m);
} // end template bdcsvd


// verify assert
template<typename MatrixType> 
void bdcsvd_verify_assert(const MatrixType& m)
{
  svd_verify_assert< MatrixType, BDCSVD< MatrixType > >(m);
}// end template bdcsvd_verify_assert


// test weird values
template<typename MatrixType>
void bdcsvd_inf_nan()
{
  svd_inf_nan< MatrixType, BDCSVD< MatrixType > >();
}// end template bdcsvd_inf_nan



void bdcsvd_preallocate()
{
  svd_preallocate< BDCSVD< MatrixXf > >();
} // end bdcsvd_preallocate


// compare the Singular values returned with Jacobi and Bdc
template<typename MatrixType> 
void compare_bdc_jacobi(const MatrixType& a = MatrixType(), unsigned int computationOptions = 0)
{
  std::cout << "debut compare" << std::endl;
  MatrixType m = MatrixType::Random(a.rows(), a.cols());
  BDCSVD<MatrixType> bdc_svd(m);
  JacobiSVD<MatrixType> jacobi_svd(m);
  VERIFY_IS_APPROX(bdc_svd.singularValues(), jacobi_svd.singularValues());
  if(computationOptions & ComputeFullU)
    VERIFY_IS_APPROX(bdc_svd.matrixU(), jacobi_svd.matrixU());
  if(computationOptions & ComputeThinU)
    VERIFY_IS_APPROX(bdc_svd.matrixU(), jacobi_svd.matrixU());
  if(computationOptions & ComputeFullV)
    VERIFY_IS_APPROX(bdc_svd.matrixV(), jacobi_svd.matrixV());
  if(computationOptions & ComputeThinV)
    VERIFY_IS_APPROX(bdc_svd.matrixV(), jacobi_svd.matrixV());
  std::cout << "fin compare" << std::endl;
} // end template compare_bdc_jacobi


// call the tests
void test_bdcsvd()
{
  // test of Dynamic defined Matrix (42, 42) of float 
  CALL_SUBTEST_11(( bdcsvd_verify_assert<Matrix<float,Dynamic,Dynamic> >
		    (Matrix<float,Dynamic,Dynamic>(42,42)) ));
  CALL_SUBTEST_11(( compare_bdc_jacobi<Matrix<float,Dynamic,Dynamic> >
		    (Matrix<float,Dynamic,Dynamic>(42,42), 0) ));
  CALL_SUBTEST_11(( bdcsvd<Matrix<float,Dynamic,Dynamic> >
		    (Matrix<float,Dynamic,Dynamic>(42,42)) ));

  // test of Dynamic defined Matrix (50, 50) of double 
  CALL_SUBTEST_13(( bdcsvd_verify_assert<Matrix<double,Dynamic,Dynamic> >
		    (Matrix<double,Dynamic,Dynamic>(50,50)) ));
  CALL_SUBTEST_13(( compare_bdc_jacobi<Matrix<double,Dynamic,Dynamic> >
		    (Matrix<double,Dynamic,Dynamic>(50,50), 0) ));
  CALL_SUBTEST_13(( bdcsvd<Matrix<double,Dynamic,Dynamic> >
		    (Matrix<double,Dynamic,Dynamic>(50, 50)) )); 

  // test of Dynamic defined Matrix (22, 22) of complex double
  CALL_SUBTEST_14(( bdcsvd_verify_assert<Matrix<std::complex<double>,Dynamic,Dynamic> >
  		    (Matrix<std::complex<double>,Dynamic,Dynamic>(22,22)) ));
  CALL_SUBTEST_14(( compare_bdc_jacobi<Matrix<std::complex<double>,Dynamic,Dynamic> >
  		    (Matrix<std::complex<double>, Dynamic, Dynamic> (22,22), 0) ));
  CALL_SUBTEST_14(( bdcsvd<Matrix<std::complex<double>,Dynamic,Dynamic> >
  		    (Matrix<std::complex<double>,Dynamic,Dynamic>(22, 22)) )); 

  // test of Dynamic defined Matrix (10, 10) of int
  //CALL_SUBTEST_15(( bdcsvd_verify_assert<Matrix<int,Dynamic,Dynamic> >
  //		    (Matrix<int,Dynamic,Dynamic>(10,10)) ));		    
  //CALL_SUBTEST_15(( compare_bdc_jacobi<Matrix<int,Dynamic,Dynamic> >
  //		    (Matrix<int,Dynamic,Dynamic>(10,10), 0) ));
  //CALL_SUBTEST_15(( bdcsvd<Matrix<int,Dynamic,Dynamic> >
  //		    (Matrix<int,Dynamic,Dynamic>(10, 10)) )); 
  

  // test of Dynamic defined Matrix (8, 6) of double 
 
  CALL_SUBTEST_16(( bdcsvd_verify_assert<Matrix<double,Dynamic,Dynamic> >
		    (Matrix<double,Dynamic,Dynamic>(8,6)) ));
  CALL_SUBTEST_16(( compare_bdc_jacobi<Matrix<double,Dynamic,Dynamic> >
		    (Matrix<double,Dynamic,Dynamic>(8, 6), 0) )); 
  CALL_SUBTEST_16(( bdcsvd<Matrix<double,Dynamic,Dynamic> >
		    (Matrix<double,Dynamic,Dynamic>(8, 6)) ));


  
  // test of Dynamic defined Matrix (36, 12) of float
  CALL_SUBTEST_17(( compare_bdc_jacobi<Matrix<float,Dynamic,Dynamic> >
		    (Matrix<float,Dynamic,Dynamic>(36, 12), 0) )); 
  CALL_SUBTEST_17(( bdcsvd<Matrix<float,Dynamic,Dynamic> >
		    (Matrix<float,Dynamic,Dynamic>(36, 12)) )); 

  // test of Dynamic defined Matrix (5, 8) of double 
  CALL_SUBTEST_18(( compare_bdc_jacobi<Matrix<double,Dynamic,Dynamic> >
		    (Matrix<double,Dynamic,Dynamic>(5, 8), 0) )); 
  CALL_SUBTEST_18(( bdcsvd<Matrix<double,Dynamic,Dynamic> >
		    (Matrix<double,Dynamic,Dynamic>(5, 8)) )); 


  // non regression tests
  CALL_SUBTEST_3(( bdcsvd_verify_assert(Matrix3f()) ));
  CALL_SUBTEST_4(( bdcsvd_verify_assert(Matrix4d()) ));
  CALL_SUBTEST_7(( bdcsvd_verify_assert(MatrixXf(10,12)) ));
  CALL_SUBTEST_8(( bdcsvd_verify_assert(MatrixXcd(7,5)) ));

  // SUBTESTS 1 and 2 on specifics matrix
  for(int i = 0; i < g_repeat; i++) {
    Matrix2cd m;
    m << 0, 1,
      0, 1;
    CALL_SUBTEST_1(( bdcsvd(m, false) ));
    m << 1, 0,
      1, 0;
    CALL_SUBTEST_1(( bdcsvd(m, false) ));

    Matrix2d n;
    n << 0, 0,
      0, 0;
    CALL_SUBTEST_2(( bdcsvd(n, false) ));
    n << 0, 0,
      0, 1;
    CALL_SUBTEST_2(( bdcsvd(n, false) ));
    
    // Statics matrix don't work with BDSVD yet
    // bdc algo on a random 3x3 float matrix
    // CALL_SUBTEST_3(( bdcsvd<Matrix3f>() ));
    // bdc algo on a random 4x4 double matrix
    // CALL_SUBTEST_4(( bdcsvd<Matrix4d>() ));
    // bdc algo on a random 3x5 float matrix
    // CALL_SUBTEST_5(( bdcsvd<Matrix<float,3,5> >() ));

    int r = internal::random<int>(1, 30),
      c = internal::random<int>(1, 30);
    CALL_SUBTEST_7(( bdcsvd<MatrixXf>(MatrixXf(r,c)) ));
    CALL_SUBTEST_8(( bdcsvd<MatrixXcd>(MatrixXcd(r,c)) ));
    (void) r;
    (void) c;

    // Test on inf/nan matrix
    CALL_SUBTEST_7( bdcsvd_inf_nan<MatrixXf>() );
  }

  CALL_SUBTEST_7(( bdcsvd<MatrixXf>(MatrixXf(internal::random<int>(EIGEN_TEST_MAX_SIZE/4, EIGEN_TEST_MAX_SIZE/2), internal::random<int>(EIGEN_TEST_MAX_SIZE/4, EIGEN_TEST_MAX_SIZE/2))) ));
  CALL_SUBTEST_8(( bdcsvd<MatrixXcd>(MatrixXcd(internal::random<int>(EIGEN_TEST_MAX_SIZE/4, EIGEN_TEST_MAX_SIZE/3), internal::random<int>(EIGEN_TEST_MAX_SIZE/4, EIGEN_TEST_MAX_SIZE/3))) ));

  // Test problem size constructors
  CALL_SUBTEST_7( BDCSVD<MatrixXf>(10,10) );

} // end test_bdcsvd
