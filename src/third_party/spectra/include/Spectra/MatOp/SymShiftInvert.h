// Copyright (C) 2020-2022 Yixuan Qiu <yixuan.qiu@cos.name>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SPECTRA_SYM_SHIFT_INVERT_H
#define SPECTRA_SYM_SHIFT_INVERT_H

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <stdexcept>
#include <type_traits>  // std::conditional, std::is_same

#include "../LinAlg/BKLDLT.h"
#include "../Util/CompInfo.h"

namespace Spectra {

/// \cond

// Compute and factorize A-sigma*B without unnecessary copying
// Default case: A is sparse, B is sparse
template <bool AIsSparse, bool BIsSparse, int UploA, int UploB>
class SymShiftInvertHelper
{
public:
    template <typename Scalar, typename Fac, typename ArgA, typename ArgB>
    static bool factorize(Fac& fac, const ArgA& A, const ArgB& B, const Scalar& sigma)
    {
        using SpMat = typename ArgA::PlainObject;
        SpMat matA = A.template selfadjointView<UploA>();
        SpMat matB = B.template selfadjointView<UploB>();
        SpMat mat = matA - sigma * matB;
        // SparseLU solver
        fac.isSymmetric(true);
        fac.compute(mat);
        // Return true if successful
        return fac.info() == Eigen::Success;
    }
};

// A is dense, B is dense or sparse
template <bool BIsSparse, int UploA, int UploB>
class SymShiftInvertHelper<false, BIsSparse, UploA, UploB>
{
public:
    template <typename Scalar, typename Fac, typename ArgA, typename ArgB>
    static bool factorize(Fac& fac, const ArgA& A, const ArgB& B, const Scalar& sigma)
    {
        using Matrix = typename ArgA::PlainObject;
        // Make a copy of the <UploA> triangular part of A
        Matrix mat(A.rows(), A.cols());
        mat.template triangularView<UploA>() = A;
        // Update <UploA> triangular part of mat
        if (UploA == UploB)
            mat -= (B * sigma).template triangularView<UploA>();
        else
            mat -= (B * sigma).template triangularView<UploB>().transpose();
        // BKLDLT solver
        fac.compute(mat, UploA);
        // Return true if successful
        return fac.info() == CompInfo::Successful;
    }
};

// A is sparse, B is dense
template <int UploA, int UploB>
class SymShiftInvertHelper<true, false, UploA, UploB>
{
public:
    template <typename Scalar, typename Fac, typename ArgA, typename ArgB>
    static bool factorize(Fac& fac, const ArgA& A, const ArgB& B, const Scalar& sigma)
    {
        using Matrix = typename ArgB::PlainObject;
        // Construct the <UploB> triangular part of -sigma*B
        Matrix mat(B.rows(), B.cols());
        mat.template triangularView<UploB>() = -sigma * B;
        // Update <UploB> triangular part of mat
        if (UploA == UploB)
            mat += A.template triangularView<UploB>();
        else
            mat += A.template triangularView<UploA>().transpose();
        // BKLDLT solver
        fac.compute(mat, UploB);
        // Return true if successful
        return fac.info() == CompInfo::Successful;
    }
};

/// \endcond

///
/// \ingroup MatOp
///
/// This class defines matrix operations required by the generalized eigen solver
/// in the shift-and-invert mode. Given two symmetric matrices \f$A\f$ and \f$B\f$,
/// it solves the linear equation \f$y=(A-\sigma B)^{-1}x\f$, where \f$\sigma\f$ is a real shift.
/// Each of \f$A\f$ and \f$B\f$ can be a dense or sparse matrix.
///
/// This class is intended to be used with the SymGEigsShiftSolver generalized eigen solver.
///
/// \tparam Scalar_        The element type of the matrices.
///                        Currently supported types are `float`, `double`, and `long double`.
/// \tparam TypeA          The type of the \f$A\f$ matrix, indicating whether \f$A\f$ is
///                        dense or sparse. Possible values are `Eigen::Dense` and `Eigen::Sparse`.
/// \tparam TypeB          The type of the \f$B\f$ matrix, indicating whether \f$B\f$ is
///                        dense or sparse. Possible values are `Eigen::Dense` and `Eigen::Sparse`.
/// \tparam UploA          Whether the lower or upper triangular part of \f$A\f$ should be used.
///                        Possible values are `Eigen::Lower` and `Eigen::Upper`.
/// \tparam UploB          Whether the lower or upper triangular part of \f$B\f$ should be used.
///                        Possible values are `Eigen::Lower` and `Eigen::Upper`.
/// \tparam FlagsA         Additional flags for the matrix class of \f$A\f$.
///                        Possible values are `Eigen::ColMajor` and `Eigen::RowMajor`.
/// \tparam FlagsB         Additional flags for the matrix class of \f$B\f$.
///                        Possible values are `Eigen::ColMajor` and `Eigen::RowMajor`.
/// \tparam StorageIndexA  The storage index type of the \f$A\f$ matrix, only used when \f$A\f$
///                        is a sparse matrix.
/// \tparam StorageIndexB  The storage index type of the \f$B\f$ matrix, only used when \f$B\f$
///                        is a sparse matrix.
///
template <typename Scalar_, typename TypeA = Eigen::Sparse, typename TypeB = Eigen::Sparse,
          int UploA = Eigen::Lower, int UploB = Eigen::Lower,
          int FlagsA = Eigen::ColMajor, int FlagsB = Eigen::ColMajor,
          typename StorageIndexA = int, typename StorageIndexB = int>
class SymShiftInvert
{
public:
    ///
    /// Element type of the matrix.
    ///
    using Scalar = Scalar_;

private:
    using Index = Eigen::Index;

    // Hypothetical type of the A matrix, either dense or sparse
    using DenseTypeA = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, FlagsA>;
    using SparseTypeA = Eigen::SparseMatrix<Scalar, FlagsA, StorageIndexA>;
    // Whether A is sparse
    using ASparse = std::is_same<TypeA, Eigen::Sparse>;
    // Actual type of the A matrix
    using MatrixA = typename std::conditional<ASparse::value, SparseTypeA, DenseTypeA>::type;

    // Hypothetical type of the B matrix, either dense or sparse
    using DenseTypeB = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, FlagsB>;
    using SparseTypeB = Eigen::SparseMatrix<Scalar, FlagsB, StorageIndexB>;
    // Whether B is sparse
    using BSparse = std::is_same<TypeB, Eigen::Sparse>;
    // Actual type of the B matrix
    using MatrixB = typename std::conditional<BSparse::value, SparseTypeB, DenseTypeB>::type;

    using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
    using MapConstVec = Eigen::Map<const Vector>;
    using MapVec = Eigen::Map<Vector>;

    // The type of A-sigma*B if one of A and B is dense
    // DenseType = if (A is dense) MatrixA else MatrixB
    using DenseType = typename std::conditional<ASparse::value, MatrixB, MatrixA>::type;
    // The type of A-sigma*B
    // If both A and B are sparse, the result is MatrixA; otherwise the result is DenseType
    using ResType = typename std::conditional<ASparse::value && BSparse::value, MatrixA, DenseType>::type;

    // If both A and B are sparse, then the result A-sigma*B is sparse, so we use
    // sparseLU for factorization; otherwise A-sigma*B is dense, and we use BKLDLT
    using FacType = typename std::conditional<
        ASparse::value && BSparse::value,
        Eigen::SparseLU<ResType>,
        BKLDLT<Scalar>>::type;

    using ConstGenericMatrixA = const Eigen::Ref<const MatrixA>;
    using ConstGenericMatrixB = const Eigen::Ref<const MatrixB>;

    ConstGenericMatrixA m_matA;
    ConstGenericMatrixB m_matB;
    const Index m_n;
    FacType m_solver;

public:
    ///
    /// Constructor to create the matrix operation object.
    ///
    /// \param A A dense or sparse matrix object, whose type can be `Eigen::Matrix<...>`,
    ///          `Eigen::SparseMatrix<...>`, `Eigen::Map<Eigen::Matrix<...>>`,
    ///          `Eigen::Map<Eigen::SparseMatrix<...>>`, `Eigen::Ref<Eigen::Matrix<...>>`,
    ///          `Eigen::Ref<Eigen::SparseMatrix<...>>`, etc.
    /// \param B A dense or sparse matrix object.
    ///
    template <typename DerivedA, typename DerivedB>
    SymShiftInvert(const Eigen::EigenBase<DerivedA>& A, const Eigen::EigenBase<DerivedB>& B) :
        m_matA(A.derived()), m_matB(B.derived()), m_n(A.rows())
    {
        static_assert(
            static_cast<int>(DerivedA::PlainObject::IsRowMajor) == static_cast<int>(MatrixA::IsRowMajor),
            "SymShiftInvert: the \"FlagsA\" template parameter does not match the input matrix (Eigen::ColMajor/Eigen::RowMajor)");

        static_assert(
            static_cast<int>(DerivedB::PlainObject::IsRowMajor) == static_cast<int>(MatrixB::IsRowMajor),
            "SymShiftInvert: the \"FlagsB\" template parameter does not match the input matrix (Eigen::ColMajor/Eigen::RowMajor)");

        if (m_n != A.cols() || m_n != B.rows() || m_n != B.cols())
            throw std::invalid_argument("SymShiftInvert: A and B must be square matrices of the same size");
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
    /// Set the real shift \f$\sigma\f$.
    ///
    void set_shift(const Scalar& sigma)
    {
        constexpr bool AIsSparse = ASparse::value;
        constexpr bool BIsSparse = BSparse::value;
        using Helper = SymShiftInvertHelper<AIsSparse, BIsSparse, UploA, UploB>;
        const bool success = Helper::factorize(m_solver, m_matA, m_matB, sigma);
        if (!success)
            throw std::invalid_argument("SymShiftInvert: factorization failed with the given shift");
    }

    ///
    /// Perform the shift-invert operation \f$y=(A-\sigma B)^{-1}x\f$.
    ///
    /// \param x_in  Pointer to the \f$x\f$ vector.
    /// \param y_out Pointer to the \f$y\f$ vector.
    ///
    // y_out = inv(A - sigma * B) * x_in
    void perform_op(const Scalar* x_in, Scalar* y_out) const
    {
        MapConstVec x(x_in, m_n);
        MapVec y(y_out, m_n);
        y.noalias() = m_solver.solve(x);
    }
};

}  // namespace Spectra

#endif  // SPECTRA_SYM_SHIFT_INVERT_H
