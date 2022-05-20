// Copyright (C) 2016-2022 Yixuan Qiu <yixuan.qiu@cos.name>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SPECTRA_SPARSE_CHOLESKY_H
#define SPECTRA_SPARSE_CHOLESKY_H

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/SparseCholesky>
#include <stdexcept>

#include "../Util/CompInfo.h"

namespace Spectra {

///
/// \ingroup MatOp
///
/// This class defines the operations related to Cholesky decomposition on a
/// sparse positive definite matrix, \f$B=LL'\f$, where \f$L\f$ is a lower triangular
/// matrix. It is mainly used in the SymGEigsSolver generalized eigen solver
/// in the Cholesky decomposition mode.
///
/// \tparam Scalar_      The element type of the matrix, for example,
///                      `float`, `double`, and `long double`.
/// \tparam Uplo         Either `Eigen::Lower` or `Eigen::Upper`, indicating which
///                      triangular part of the matrix is used.
/// \tparam Flags        Either `Eigen::ColMajor` or `Eigen::RowMajor`, indicating
///                      the storage format of the input matrix.
/// \tparam StorageIndex The type of the indices for the sparse matrix.
///
template <typename Scalar_, int Uplo = Eigen::Lower, int Flags = Eigen::ColMajor, typename StorageIndex = int>
class SparseCholesky
{
public:
    ///
    /// Element type of the matrix.
    ///
    using Scalar = Scalar_;

private:
    using Index = Eigen::Index;
    using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
    using MapConstVec = Eigen::Map<const Vector>;
    using MapVec = Eigen::Map<Vector>;
    using SparseMatrix = Eigen::SparseMatrix<Scalar, Flags, StorageIndex>;

    const Index m_n;
    Eigen::SimplicialLLT<SparseMatrix, Uplo> m_decomp;
    CompInfo m_info;  // status of the decomposition

public:
    ///
    /// Constructor to create the matrix operation object.
    ///
    /// \param mat An **Eigen** sparse matrix object, whose type can be
    /// `Eigen::SparseMatrix<Scalar, ...>` or its mapped version
    /// `Eigen::Map<Eigen::SparseMatrix<Scalar, ...> >`.
    ///
    template <typename Derived>
    SparseCholesky(const Eigen::SparseMatrixBase<Derived>& mat) :
        m_n(mat.rows())
    {
        static_assert(
            static_cast<int>(Derived::PlainObject::IsRowMajor) == static_cast<int>(SparseMatrix::IsRowMajor),
            "SparseCholesky: the \"Flags\" template parameter does not match the input matrix (Eigen::ColMajor/Eigen::RowMajor)");

        if (mat.rows() != mat.cols())
            throw std::invalid_argument("SparseCholesky: matrix must be square");

        m_decomp.compute(mat);
        m_info = (m_decomp.info() == Eigen::Success) ?
            CompInfo::Successful :
            CompInfo::NumericalIssue;
    }

    ///
    /// Returns the number of rows of the underlying matrix.
    ///
    Index rows() const { return m_n; }
    ///
    /// Returns the number of columns of the underlying matrix.
    ///
    Index cols() const { return m_n; }

    ///
    /// Returns the status of the computation.
    /// The full list of enumeration values can be found in \ref Enumerations.
    ///
    CompInfo info() const { return m_info; }

    ///
    /// Performs the lower triangular solving operation \f$y=L^{-1}x\f$.
    ///
    /// \param x_in  Pointer to the \f$x\f$ vector.
    /// \param y_out Pointer to the \f$y\f$ vector.
    ///
    // y_out = inv(L) * x_in
    void lower_triangular_solve(const Scalar* x_in, Scalar* y_out) const
    {
        MapConstVec x(x_in, m_n);
        MapVec y(y_out, m_n);
        y.noalias() = m_decomp.permutationP() * x;
        m_decomp.matrixL().solveInPlace(y);
    }

    ///
    /// Performs the upper triangular solving operation \f$y=(L')^{-1}x\f$.
    ///
    /// \param x_in  Pointer to the \f$x\f$ vector.
    /// \param y_out Pointer to the \f$y\f$ vector.
    ///
    // y_out = inv(L') * x_in
    void upper_triangular_solve(const Scalar* x_in, Scalar* y_out) const
    {
        MapConstVec x(x_in, m_n);
        MapVec y(y_out, m_n);
        y.noalias() = m_decomp.matrixU().solve(x);
        y = m_decomp.permutationPinv() * y;
    }
};

}  // namespace Spectra

#endif  // SPECTRA_SPARSE_CHOLESKY_H
