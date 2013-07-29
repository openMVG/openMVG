// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2009-2010 Benoit Jacob <jacob.benoit.1@gmail.com>
//
// Copyright (C) 2013 Gauthier Brun <brun.gauthier@gmail.com>
// Copyright (C) 2013 Nicolas Carre <nicolas.carre@ensimag.fr>
// Copyright (C) 2013 Jean Ceccato <jean.ceccato@ensimag.fr>
// Copyright (C) 2013 Pierre Zoppitelli <pierre.zoppitelli@ensimag.fr>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef EIGEN_SVD_H
#define EIGEN_SVD_H

namespace Eigen {
/** \ingroup SVD_Module
 *
 *
 * \class SVDBase
 *
 * \brief Mother class of SVD classes algorithms
 *
 * \param MatrixType the type of the matrix of which we are computing the SVD decomposition
 * SVD decomposition consists in decomposing any n-by-p matrix \a A as a product
 *   \f[ A = U S V^* \f]
 * where \a U is a n-by-n unitary, \a V is a p-by-p unitary, and \a S is a n-by-p real positive matrix which is zero outside of its main diagonal;
 * the diagonal entries of S are known as the \em singular \em values of \a A and the columns of \a U and \a V are known as the left
 * and right \em singular \em vectors of \a A respectively.
 *
 * Singular values are always sorted in decreasing order.
 *
 * 
 * You can ask for only \em thin \a U or \a V to be computed, meaning the following. In case of a rectangular n-by-p matrix, letting \a m be the
 * smaller value among \a n and \a p, there are only \a m singular vectors; the remaining columns of \a U and \a V do not correspond to actual
 * singular vectors. Asking for \em thin \a U or \a V means asking for only their \a m first columns to be formed. So \a U is then a n-by-m matrix,
 * and \a V is then a p-by-m matrix. Notice that thin \a U and \a V are all you need for (least squares) solving.
 *  
 * If the input matrix has inf or nan coefficients, the result of the computation is undefined, but the computation is guaranteed to
 * terminate in finite (and reasonable) time.
 * \sa MatrixBase::genericSvd()
 */
template<typename _MatrixType> 
class SVDBase
{

public:
  typedef _MatrixType MatrixType;
  typedef typename MatrixType::Scalar Scalar;
  typedef typename NumTraits<typename MatrixType::Scalar>::Real RealScalar;
  typedef typename MatrixType::Index Index;
  enum {
    RowsAtCompileTime = MatrixType::RowsAtCompileTime,
    ColsAtCompileTime = MatrixType::ColsAtCompileTime,
    DiagSizeAtCompileTime = EIGEN_SIZE_MIN_PREFER_DYNAMIC(RowsAtCompileTime,ColsAtCompileTime),
    MaxRowsAtCompileTime = MatrixType::MaxRowsAtCompileTime,
    MaxColsAtCompileTime = MatrixType::MaxColsAtCompileTime,
    MaxDiagSizeAtCompileTime = EIGEN_SIZE_MIN_PREFER_FIXED(MaxRowsAtCompileTime,MaxColsAtCompileTime),
    MatrixOptions = MatrixType::Options
  };

  typedef Matrix<Scalar, RowsAtCompileTime, RowsAtCompileTime,
		 MatrixOptions, MaxRowsAtCompileTime, MaxRowsAtCompileTime>
  MatrixUType;
  typedef Matrix<Scalar, ColsAtCompileTime, ColsAtCompileTime,
		 MatrixOptions, MaxColsAtCompileTime, MaxColsAtCompileTime>
  MatrixVType;
  typedef typename internal::plain_diag_type<MatrixType, RealScalar>::type SingularValuesType;
  typedef typename internal::plain_row_type<MatrixType>::type RowType;
  typedef typename internal::plain_col_type<MatrixType>::type ColType;
  typedef Matrix<Scalar, DiagSizeAtCompileTime, DiagSizeAtCompileTime,
		 MatrixOptions, MaxDiagSizeAtCompileTime, MaxDiagSizeAtCompileTime>
  WorkMatrixType;
	



  /** \brief Method performing the decomposition of given matrix using custom options.
   *
   * \param matrix the matrix to decompose
   * \param computationOptions optional parameter allowing to specify if you want full or thin U or V unitaries to be computed.
   *                           By default, none is computed. This is a bit-field, the possible bits are #ComputeFullU, #ComputeThinU,
   *                           #ComputeFullV, #ComputeThinV.
   *
   * Thin unitaries are only available if your matrix type has a Dynamic number of columns (for example MatrixXf). They also are not
   * available with the (non-default) FullPivHouseholderQR preconditioner.
   */
  SVDBase& compute(const MatrixType& matrix, unsigned int computationOptions);

  /** \brief Method performing the decomposition of given matrix using current options.
   *
   * \param matrix the matrix to decompose
   *
   * This method uses the current \a computationOptions, as already passed to the constructor or to compute(const MatrixType&, unsigned int).
   */
  //virtual SVDBase& compute(const MatrixType& matrix) = 0;
  SVDBase& compute(const MatrixType& matrix);

  /** \returns the \a U matrix.
   *
   * For the SVDBase decomposition of a n-by-p matrix, letting \a m be the minimum of \a n and \a p,
   * the U matrix is n-by-n if you asked for #ComputeFullU, and is n-by-m if you asked for #ComputeThinU.
   *
   * The \a m first columns of \a U are the left singular vectors of the matrix being decomposed.
   *
   * This method asserts that you asked for \a U to be computed.
   */
  const MatrixUType& matrixU() const
  {
    eigen_assert(m_isInitialized && "SVD is not initialized.");
    eigen_assert(computeU() && "This SVD decomposition didn't compute U. Did you ask for it?");
    return m_matrixU;
  }

  /** \returns the \a V matrix.
   *
   * For the SVD decomposition of a n-by-p matrix, letting \a m be the minimum of \a n and \a p,
   * the V matrix is p-by-p if you asked for #ComputeFullV, and is p-by-m if you asked for ComputeThinV.
   *
   * The \a m first columns of \a V are the right singular vectors of the matrix being decomposed.
   *
   * This method asserts that you asked for \a V to be computed.
   */
  const MatrixVType& matrixV() const
  {
    eigen_assert(m_isInitialized && "SVD is not initialized.");
    eigen_assert(computeV() && "This SVD decomposition didn't compute V. Did you ask for it?");
    return m_matrixV;
  }

  /** \returns the vector of singular values.
   *
   * For the SVD decomposition of a n-by-p matrix, letting \a m be the minimum of \a n and \a p, the
   * returned vector has size \a m.  Singular values are always sorted in decreasing order.
   */
  const SingularValuesType& singularValues() const
  {
    eigen_assert(m_isInitialized && "SVD is not initialized.");
    return m_singularValues;
  }

  

  /** \returns the number of singular values that are not exactly 0 */
  Index nonzeroSingularValues() const
  {
    eigen_assert(m_isInitialized && "SVD is not initialized.");
    return m_nonzeroSingularValues;
  }


  /** \returns true if \a U (full or thin) is asked for in this SVD decomposition */
  inline bool computeU() const { return m_computeFullU || m_computeThinU; }
  /** \returns true if \a V (full or thin) is asked for in this SVD decomposition */
  inline bool computeV() const { return m_computeFullV || m_computeThinV; }


  inline Index rows() const { return m_rows; }
  inline Index cols() const { return m_cols; }


protected:
  // return true if already allocated
  bool allocate(Index rows, Index cols, unsigned int computationOptions) ;

  MatrixUType m_matrixU;
  MatrixVType m_matrixV;
  SingularValuesType m_singularValues;
  bool m_isInitialized, m_isAllocated;
  bool m_computeFullU, m_computeThinU;
  bool m_computeFullV, m_computeThinV;
  unsigned int m_computationOptions;
  Index m_nonzeroSingularValues, m_rows, m_cols, m_diagSize;


  /** \brief Default Constructor.
   *
   * Default constructor of SVDBase
   */
  SVDBase()
    : m_isInitialized(false),
      m_isAllocated(false),
      m_computationOptions(0),
      m_rows(-1), m_cols(-1)
  {}


};


template<typename MatrixType>
bool SVDBase<MatrixType>::allocate(Index rows, Index cols, unsigned int computationOptions)
{
  eigen_assert(rows >= 0 && cols >= 0);

  if (m_isAllocated &&
      rows == m_rows &&
      cols == m_cols &&
      computationOptions == m_computationOptions)
  {
    return true;
  }

  m_rows = rows;
  m_cols = cols;
  m_isInitialized = false;
  m_isAllocated = true;
  m_computationOptions = computationOptions;
  m_computeFullU = (computationOptions & ComputeFullU) != 0;
  m_computeThinU = (computationOptions & ComputeThinU) != 0;
  m_computeFullV = (computationOptions & ComputeFullV) != 0;
  m_computeThinV = (computationOptions & ComputeThinV) != 0;
  eigen_assert(!(m_computeFullU && m_computeThinU) && "SVDBase: you can't ask for both full and thin U");
  eigen_assert(!(m_computeFullV && m_computeThinV) && "SVDBase: you can't ask for both full and thin V");
  eigen_assert(EIGEN_IMPLIES(m_computeThinU || m_computeThinV, MatrixType::ColsAtCompileTime==Dynamic) &&
	       "SVDBase: thin U and V are only available when your matrix has a dynamic number of columns.");

  m_diagSize = (std::min)(m_rows, m_cols);
  m_singularValues.resize(m_diagSize);
  if(RowsAtCompileTime==Dynamic)
    m_matrixU.resize(m_rows, m_computeFullU ? m_rows
		     : m_computeThinU ? m_diagSize
		     : 0);
  if(ColsAtCompileTime==Dynamic)
    m_matrixV.resize(m_cols, m_computeFullV ? m_cols
		     : m_computeThinV ? m_diagSize
		     : 0);

  return false;
}

}// end namespace

#endif // EIGEN_SVD_H
