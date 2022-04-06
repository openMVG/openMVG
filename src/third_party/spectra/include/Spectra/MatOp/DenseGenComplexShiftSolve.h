// Copyright (C) 2016-2022 Yixuan Qiu <yixuan.qiu@cos.name>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SPECTRA_DENSE_GEN_COMPLEX_SHIFT_SOLVE_H
#define SPECTRA_DENSE_GEN_COMPLEX_SHIFT_SOLVE_H

#include <Eigen/Core>
#include <Eigen/LU>
#include <stdexcept>

namespace Spectra {

///
/// \ingroup MatOp
///
/// This class defines the complex shift-solve operation on a general real matrix \f$A\f$,
/// i.e., calculating \f$y=\mathrm{Re}\{(A-\sigma I)^{-1}x\}\f$ for any complex-valued
/// \f$\sigma\f$ and real-valued vector \f$x\f$. It is mainly used in the
/// GenEigsComplexShiftSolver eigen solver.
///
/// \tparam Scalar_ The element type of the matrix, for example,
///                 `float`, `double`, and `long double`.
/// \tparam Flags   Either `Eigen::ColMajor` or `Eigen::RowMajor`, indicating
///                 the storage format of the input matrix.
///
template <typename Scalar_, int Flags = Eigen::ColMajor>
class DenseGenComplexShiftSolve
{
public:
    ///
    /// Element type of the matrix.
    ///
    using Scalar = Scalar_;

private:
    using Index = Eigen::Index;
    using Matrix = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Flags>;
    using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
    using MapConstVec = Eigen::Map<const Vector>;
    using MapVec = Eigen::Map<Vector>;
    using ConstGenericMatrix = const Eigen::Ref<const Matrix>;

    using Complex = std::complex<Scalar>;
    using ComplexMatrix = Eigen::Matrix<Complex, Eigen::Dynamic, Eigen::Dynamic, Flags>;
    using ComplexVector = Eigen::Matrix<Complex, Eigen::Dynamic, 1>;

    using ComplexSolver = Eigen::PartialPivLU<ComplexMatrix>;

    ConstGenericMatrix m_mat;
    const Index m_n;
    ComplexSolver m_solver;
    mutable ComplexVector m_x_cache;

public:
    ///
    /// Constructor to create the matrix operation object.
    ///
    /// \param mat An **Eigen** matrix object, whose type can be
    /// `Eigen::Matrix<Scalar, ...>` (e.g. `Eigen::MatrixXd` and
    /// `Eigen::MatrixXf`), or its mapped version
    /// (e.g. `Eigen::Map<Eigen::MatrixXd>`).
    ///
    template <typename Derived>
    DenseGenComplexShiftSolve(const Eigen::MatrixBase<Derived>& mat) :
        m_mat(mat), m_n(mat.rows())
    {
        static_assert(
            static_cast<int>(Derived::PlainObject::IsRowMajor) == static_cast<int>(Matrix::IsRowMajor),
            "DenseGenComplexShiftSolve: the \"Flags\" template parameter does not match the input matrix (Eigen::ColMajor/Eigen::RowMajor)");

        if (mat.rows() != mat.cols())
            throw std::invalid_argument("DenseGenComplexShiftSolve: matrix must be square");
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
        m_solver.compute(m_mat.template cast<Complex>() - Complex(sigmar, sigmai) * ComplexMatrix::Identity(m_n, m_n));
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

#endif  // SPECTRA_DENSE_GEN_COMPLEX_SHIFT_SOLVE_H
