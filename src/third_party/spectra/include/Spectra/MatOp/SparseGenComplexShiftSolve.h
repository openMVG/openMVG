// Copyright (C) 2020-2022 Yixuan Qiu <yixuan.qiu@cos.name>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SPECTRA_SPARSE_GEN_COMPLEX_SHIFT_SOLVE_H
#define SPECTRA_SPARSE_GEN_COMPLEX_SHIFT_SOLVE_H

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <stdexcept>

namespace Spectra {

///
/// \ingroup MatOp
///
/// This class defines the complex shift-solve operation on a sparse real matrix \f$A\f$,
/// i.e., calculating \f$y=\mathrm{Re}\{(A-\sigma I)^{-1}x\}\f$ for any complex-valued
/// \f$\sigma\f$ and real-valued vector \f$x\f$. It is mainly used in the
/// GenEigsComplexShiftSolver eigen solver.
///
/// \tparam Scalar_      The element type of the matrix, for example,
///                      `float`, `double`, and `long double`.
/// \tparam Flags        Either `Eigen::ColMajor` or `Eigen::RowMajor`, indicating
///                      the storage format of the input matrix.
/// \tparam StorageIndex The type of the indices for the sparse matrix.
///
template <typename Scalar_, int Flags = Eigen::ColMajor, typename StorageIndex = int>
class SparseGenComplexShiftSolve
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

    using Complex = std::complex<Scalar>;
    using ComplexVector = Eigen::Matrix<Complex, Eigen::Dynamic, 1>;
    using SparseComplexMatrix = Eigen::SparseMatrix<Complex, Flags, StorageIndex>;

    using ComplexSolver = Eigen::SparseLU<SparseComplexMatrix>;

    ConstGenericSparseMatrix m_mat;
    const Index m_n;
    ComplexSolver m_solver;
    mutable ComplexVector m_x_cache;

public:
    ///
    /// Constructor to create the matrix operation object.
    ///
    /// \param mat An **Eigen** sparse matrix object, whose type can be
    /// `Eigen::SparseMatrix<Scalar, ...>` or its mapped version
    /// `Eigen::Map<Eigen::SparseMatrix<Scalar, ...> >`.
    ///
    template <typename Derived>
    SparseGenComplexShiftSolve(const Eigen::SparseMatrixBase<Derived>& mat) :
        m_mat(mat), m_n(mat.rows())
    {
        static_assert(
            static_cast<int>(Derived::PlainObject::IsRowMajor) == static_cast<int>(SparseMatrix::IsRowMajor),
            "SparseGenComplexShiftSolve: the \"Flags\" template parameter does not match the input matrix (Eigen::ColMajor/Eigen::RowMajor)");

        if (mat.rows() != mat.cols())
            throw std::invalid_argument("SparseGenComplexShiftSolve: matrix must be square");
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
    /// Set the complex shift \f$\sigma\f$.
    ///
    /// \param sigmar Real part of \f$\sigma\f$.
    /// \param sigmai Imaginary part of \f$\sigma\f$.
    ///
    void set_shift(const Scalar& sigmar, const Scalar& sigmai)
    {
        // Create a sparse idendity matrix (1 + 0i on diagonal)
        SparseComplexMatrix I(m_n, m_n);
        I.setIdentity();
        // Sparse LU decomposition
        m_solver.compute(m_mat.template cast<Complex>() - Complex(sigmar, sigmai) * I);
        // Set cache to zero
        m_x_cache.resize(m_n);
        m_x_cache.setZero();
    }

    ///
    /// Perform the complex shift-solve operation
    /// \f$y=\mathrm{Re}\{(A-\sigma I)^{-1}x\}\f$.
    ///
    /// \param x_in  Pointer to the \f$x\f$ vector.
    /// \param y_out Pointer to the \f$y\f$ vector.
    ///
    // y_out = Re( inv(A - sigma * I) * x_in )
    void perform_op(const Scalar* x_in, Scalar* y_out) const
    {
        m_x_cache.real() = MapConstVec(x_in, m_n);
        MapVec y(y_out, m_n);
        y.noalias() = m_solver.solve(m_x_cache).real();
    }
};

}  // namespace Spectra

#endif  // SPECTRA_SPARSE_GEN_COMPLEX_SHIFT_SOLVE_H
