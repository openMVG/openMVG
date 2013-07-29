// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2009-2011 Jitse Niesen <jitse@maths.leeds.ac.uk>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "main.h"
#include <unsupported/Eigen/MatrixFunctions>

template <typename MatrixType, int IsComplex = NumTraits<typename internal::traits<MatrixType>::Scalar>::IsComplex>
struct generateTestMatrix;

// for real matrices, make sure none of the eigenvalues are negative
template <typename MatrixType>
struct generateTestMatrix<MatrixType,0>
{
  static void run(MatrixType& result, typename MatrixType::Index size)
  {
    MatrixType mat = MatrixType::Random(size, size);
    EigenSolver<MatrixType> es(mat);
    typename EigenSolver<MatrixType>::EigenvalueType eivals = es.eigenvalues();
    for (typename MatrixType::Index i = 0; i < size; ++i) {
      if (eivals(i).imag() == 0 && eivals(i).real() < 0)
	eivals(i) = -eivals(i);
    }
    result = (es.eigenvectors() * eivals.asDiagonal() * es.eigenvectors().inverse()).real();
  }
};

// for complex matrices, any matrix is fine
template <typename MatrixType>
struct generateTestMatrix<MatrixType,1>
{
  static void run(MatrixType& result, typename MatrixType::Index size)
  {
    result = MatrixType::Random(size, size);
  }
};

template <typename Derived, typename OtherDerived>
double relerr(const MatrixBase<Derived>& A, const MatrixBase<OtherDerived>& B)
{
  return std::sqrt((A - B).cwiseAbs2().sum() / (std::min)(A.cwiseAbs2().sum(), B.cwiseAbs2().sum()));
}
