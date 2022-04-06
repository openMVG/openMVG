// Copyright (C) 2017-2022 Yixuan Qiu <yixuan.qiu@cos.name>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SPECTRA_SPARSE_REGULAR_INVERSE_H
#define SPECTRA_SPARSE_REGULAR_INVERSE_H

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
/// \tparam Scalar_      The element type of the matrix, for example,
///                      `float`, `double`, and `long double`.
/// \tparam Uplo         Either `Eigen::Lower` or `Eigen::Upper`, indicating which
///                      triangular part of the matrix is used.
/// \tparam Flags        Either `Eigen::ColMajor` or `Eigen::RowMajor`, indicating
///                      the storage format of the input matrix.
/// \tparam StorageIndex The type of the indices for the sparse matrix.
///
template <typename Scalar_, int Uplo = Eigen::Lower, int Flags = Eigen::ColMajor, typename StorageIndex = int>
class SparseRegularInverse
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
    using ConstGenericSparseMatrix = const Eigen::Ref<const SparseMatrix>;

    ConstGenericSparseMatrix m_mat;
    const Index m_n;
    Eigen::ConjugateGradient<SparseMatrix> m_cg;
    mutable CompInfo m_info;

public:
    ///
    /// Constructor to create the matrix operation object.
    ///
    /// \param mat An **Eigen** sparse matrix object, whose type can be
    /// `Eigen::SparseMatrix<Scalar, ...>` or its mapped version
    /// `Eigen::Map<Eigen::SparseMatrix<Scalar, ...> >`.
    ///
    template <typename Derived>
    SparseRegularInverse(const Eigen::SparseMatrixBase<Derived>& mat) :
        m_mat(mat), m_n(mat.rows())
    {
        static_assert(
            static_cast<int>(Derived::PlainObject::IsRowMajor) == static_cast<int>(SparseMatrix::IsRowMajor),
            "SparseRegularInverse: the \"Flags\" template parameter does not match the input matrix (Eigen::ColMajor/Eigen::RowMajor)");

        if (mat.rows() != mat.cols())
            throw std::invalid_argument("SparseRegularInverse: matrix must be square");

        m_cg.compute(mat);
        m_info = (m_cg.info() == Eigen::Success) ?
            CompInfo::Successful :
            CompInfo::NumericalIssue;
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
    /// Returns the status of the computation.
    /// The full list of enumeration values can be found in \ref Enumerations.
    ///
    CompInfo info() const { return m_info; }

    ///
    /// Perform the solving operation \f$y=B^{-1}x\f$.
    ///
    /// \param x_in  Pointer to the \f$x\f$ vector.
    /// \param y_out Pointer to the \f$y\f$ vector.
    ///
    // y_out = inv(B) * x_in
    void solve(const Scalar* x_in, Scalar* y_out) const
    {
        MapConstVec x(x_in, m_n);
        MapVec y(y_out, m_n);
        y.noalias() = m_cg.solve(x);

        m_info = (m_cg.info() == Eigen::Success) ?
            CompInfo::Successful :
            CompInfo::NotConverging;
        if (m_info != CompInfo::Successful)
            throw std::runtime_error("SparseRegularInverse: CG solver does not converge");
    }

    ///
    /// Perform the matrix-vector multiplication operation \f$y=Bx\f$.
    ///
    /// \param x_in  Pointer to the \f$x\f$ vector.
    /// \param y_out Pointer to the \f$y\f$ vector.
    ///
    // y_out = B * x_in
    void perform_op(const Scalar* x_in, Scalar* y_out) const
    {
        MapConstVec x(x_in, m_n);
        MapVec y(y_out, m_n);
        y.noalias() = m_mat.template selfadjointView<Uplo>() * x;
    }
};

}  // namespace Spectra

#endif  // SPECTRA_SPARSE_REGULAR_INVERSE_H
