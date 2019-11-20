// Copyright (C) 2017-2019 Yixuan Qiu <yixuan.qiu@cos.name>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SPARSE_REGULAR_INVERSE_H
#define SPARSE_REGULAR_INVERSE_H

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>
#include <stdexcept>

namespace Spectra {


///
/// \ingroup MatOp
///
/// This class defines matrix operations required by the generalized eigen solver
/// in the regular inverse mode. For a sparse and positive definite matrix \f$B\f$,
/// it implements the matrix-vector product \f$y=Bx\f$ and the linear equation
/// solving operation \f$y=B^{-1}x\f$.
///
/// This class is intended to be used with the SymGEigsSolver generalized eigen solver
/// in the regular inverse mode.
///
template <typename Scalar, int Uplo = Eigen::Lower, int Flags = 0, typename StorageIndex = int>
class SparseRegularInverse
{
private:
    typedef Eigen::Index Index;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;
    typedef Eigen::Map<const Vector> MapConstVec;
    typedef Eigen::Map<Vector> MapVec;
    typedef Eigen::SparseMatrix<Scalar, Flags, StorageIndex> SparseMatrix;
    typedef const Eigen::Ref<const SparseMatrix> ConstGenericSparseMatrix;

    ConstGenericSparseMatrix m_mat;
    const int m_n;
    Eigen::ConjugateGradient<SparseMatrix> m_cg;

public:
    ///
    /// Constructor to create the matrix operation object.
    ///
    /// \param mat An **Eigen** sparse matrix object, whose type can be
    /// `Eigen::SparseMatrix<Scalar, ...>` or its mapped version
    /// `Eigen::Map<Eigen::SparseMatrix<Scalar, ...> >`.
    ///
    SparseRegularInverse(ConstGenericSparseMatrix& mat) :
        m_mat(mat), m_n(mat.rows())
    {
        if(mat.rows() != mat.cols())
            throw std::invalid_argument("SparseRegularInverse: matrix must be square");

        m_cg.compute(mat);
    }

    ///
    /// Return the number of rows of the underlying matrix.
    ///
    Index rows() const { return m_n; }
    ///
    /// Return the number of columns of the underlying matrix.
    ///
    Index cols() const { return m_n; }

    ///
    /// Perform the solving operation \f$y=B^{-1}x\f$.
    ///
    /// \param x_in  Pointer to the \f$x\f$ vector.
    /// \param y_out Pointer to the \f$y\f$ vector.
    ///
    // y_out = inv(B) * x_in
    void solve(const Scalar* x_in, Scalar* y_out) const
    {
        MapConstVec x(x_in,  m_n);
        MapVec      y(y_out, m_n);
        y.noalias() = m_cg.solve(x);
    }

    ///
    /// Perform the matrix-vector multiplication operation \f$y=Bx\f$.
    ///
    /// \param x_in  Pointer to the \f$x\f$ vector.
    /// \param y_out Pointer to the \f$y\f$ vector.
    ///
    // y_out = B * x_in
    void mat_prod(const Scalar* x_in, Scalar* y_out) const
    {
        MapConstVec x(x_in,  m_n);
        MapVec      y(y_out, m_n);
        y.noalias() = m_mat.template selfadjointView<Uplo>() * x;
    }
};


} // namespace Spectra

#endif // SPARSE_REGULAR_INVERSE_H
