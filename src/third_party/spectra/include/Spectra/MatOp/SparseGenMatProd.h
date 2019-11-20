// Copyright (C) 2016-2019 Yixuan Qiu <yixuan.qiu@cos.name>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SPARSE_GEN_MAT_PROD_H
#define SPARSE_GEN_MAT_PROD_H

#include <Eigen/Core>
#include <Eigen/SparseCore>

namespace Spectra {


///
/// \ingroup MatOp
///
/// This class defines the matrix-vector multiplication operation on a
/// sparse real matrix \f$A\f$, i.e., calculating \f$y=Ax\f$ for any vector
/// \f$x\f$. It is mainly used in the GenEigsSolver and SymEigsSolver
/// eigen solvers.
///
template <typename Scalar, int Flags = 0, typename StorageIndex = int>
class SparseGenMatProd
{
private:
    typedef Eigen::Index Index;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;
    typedef Eigen::Map<const Vector> MapConstVec;
    typedef Eigen::Map<Vector> MapVec;
    typedef Eigen::SparseMatrix<Scalar, Flags, StorageIndex> SparseMatrix;
    typedef const Eigen::Ref<const SparseMatrix> ConstGenericSparseMatrix;

    ConstGenericSparseMatrix m_mat;

public:
    ///
    /// Constructor to create the matrix operation object.
    ///
    /// \param mat An **Eigen** sparse matrix object, whose type can be
    /// `Eigen::SparseMatrix<Scalar, ...>` or its mapped version
    /// `Eigen::Map<Eigen::SparseMatrix<Scalar, ...> >`.
    ///
    SparseGenMatProd(ConstGenericSparseMatrix& mat) :
        m_mat(mat)
    {}

    ///
    /// Return the number of rows of the underlying matrix.
    ///
    Index rows() const { return m_mat.rows(); }
    ///
    /// Return the number of columns of the underlying matrix.
    ///
    Index cols() const { return m_mat.cols(); }

    ///
    /// Perform the matrix-vector multiplication operation \f$y=Ax\f$.
    ///
    /// \param x_in  Pointer to the \f$x\f$ vector.
    /// \param y_out Pointer to the \f$y\f$ vector.
    ///
    // y_out = A * x_in
    void perform_op(const Scalar* x_in, Scalar* y_out) const
    {
        MapConstVec x(x_in,  m_mat.cols());
        MapVec      y(y_out, m_mat.rows());
        y.noalias() = m_mat * x;
    }
};


} // namespace Spectra

#endif // SPARSE_GEN_MAT_PROD_H
