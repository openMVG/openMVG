// Copyright (C) 2016-2022 Yixuan Qiu <yixuan.qiu@cos.name>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SPECTRA_SYM_EIGS_SHIFT_SOLVER_H
#define SPECTRA_SYM_EIGS_SHIFT_SOLVER_H

#include <Eigen/Core>

#include "SymEigsBase.h"
#include "Util/SelectionRule.h"
#include "MatOp/DenseSymShiftSolve.h"

namespace Spectra {

///
/// \ingroup EigenSolver
///
/// This class implements the eigen solver for real symmetric matrices using
/// the **shift-and-invert mode**. The background information of the symmetric
/// eigen solver is documented in the SymEigsSolver class. Here we focus on
/// explaining the shift-and-invert mode.
///
/// The shift-and-invert mode is based on the following fact:
/// If \f$\lambda\f$ and \f$x\f$ are a pair of eigenvalue and eigenvector of
/// matrix \f$A\f$, such that \f$Ax=\lambda x\f$, then for any \f$\sigma\f$,
/// we have
/// \f[(A-\sigma I)^{-1}x=\nu x\f]
/// where
/// \f[\nu=\frac{1}{\lambda-\sigma}\f]
/// which indicates that \f$(\nu, x)\f$ is an eigenpair of the matrix
/// \f$(A-\sigma I)^{-1}\f$.
///
/// Therefore, if we pass the matrix operation \f$(A-\sigma I)^{-1}y\f$
/// (rather than \f$Ay\f$) to the eigen solver, then we would get the desired
/// values of \f$\nu\f$, and \f$\lambda\f$ can also be easily obtained by noting
/// that \f$\lambda=\sigma+\nu^{-1}\f$.
///
/// The reason why we need this type of manipulation is that
/// the algorithm of **Spectra** (and also **ARPACK**)
/// is good at finding eigenvalues with large magnitude, but may fail in looking
/// for eigenvalues that are close to zero. However, if we really need them, we
/// can set \f$\sigma=0\f$, find the largest eigenvalues of \f$A^{-1}\f$, and then
/// transform back to \f$\lambda\f$, since in this case largest values of \f$\nu\f$
/// implies smallest values of \f$\lambda\f$.
///
/// To summarize, in the shift-and-invert mode, the selection rule will apply to
/// \f$\nu=1/(\lambda-\sigma)\f$ rather than \f$\lambda\f$. So a selection rule
/// of `LARGEST_MAGN` combined with shift \f$\sigma\f$ will find eigenvalues of
/// \f$A\f$ that are closest to \f$\sigma\f$. But note that the eigenvalues()
/// method will always return the eigenvalues in the original problem (i.e.,
/// returning \f$\lambda\f$ rather than \f$\nu\f$), and eigenvectors are the
/// same for both the original problem and the shifted-and-inverted problem.
///
/// \tparam OpType  The name of the matrix operation class. Users could either
///                 use the wrapper classes such as DenseSymShiftSolve and
///                 SparseSymShiftSolve, or define their own that implements the type
///                 definition `Scalar` and all the public member functions as in
///                 DenseSymShiftSolve.
///
/// Below is an example that illustrates the use of the shift-and-invert mode:
///
/// \code{.cpp}
/// #include <Eigen/Core>
/// #include <Spectra/SymEigsShiftSolver.h>
/// // <Spectra/MatOp/DenseSymShiftSolve.h> is implicitly included
/// #include <iostream>
///
/// using namespace Spectra;
///
/// int main()
/// {
///     // A size-10 diagonal matrix with elements 1, 2, ..., 10
///     Eigen::MatrixXd M = Eigen::MatrixXd::Zero(10, 10);
///     for (int i = 0; i < M.rows(); i++)
///         M(i, i) = i + 1;
///
///     // Construct matrix operation object using the wrapper class
///     DenseSymShiftSolve<double> op(M);
///
///     // Construct eigen solver object with shift 0
///     // This will find eigenvalues that are closest to 0
///     SymEigsShiftSolver<DenseSymShiftSolve<double>> eigs(op, 3, 6, 0.0);
///
///     eigs.init();
///     eigs.compute(SortRule::LargestMagn);
///     if (eigs.info() == CompInfo::Successful)
///     {
///         Eigen::VectorXd evalues = eigs.eigenvalues();
///         // Will get (3.0, 2.0, 1.0)
///         std::cout << "Eigenvalues found:\n" << evalues << std::endl;
///     }
///
///     return 0;
/// }
/// \endcode
///
/// Also an example for user-supplied matrix shift-solve operation class:
///
/// \code{.cpp}
/// #include <Eigen/Core>
/// #include <Spectra/SymEigsShiftSolver.h>
/// #include <iostream>
///
/// using namespace Spectra;
///
/// // M = diag(1, 2, ..., 10)
/// class MyDiagonalTenShiftSolve
/// {
/// private:
///     double sigma_;
/// public:
///     using Scalar = double;  // A typedef named "Scalar" is required
///     int rows() const { return 10; }
///     int cols() const { return 10; }
///     void set_shift(double sigma) { sigma_ = sigma; }
///     // y_out = inv(A - sigma * I) * x_in
///     // inv(A - sigma * I) = diag(1/(1-sigma), 1/(2-sigma), ...)
///     void perform_op(double *x_in, double *y_out) const
///     {
///         for (int i = 0; i < rows(); i++)
///         {
///             y_out[i] = x_in[i] / (i + 1 - sigma_);
///         }
///     }
/// };
///
/// int main()
/// {
///     MyDiagonalTenShiftSolve op;
///     // Find three eigenvalues that are closest to 3.14
///     SymEigsShiftSolver<MyDiagonalTenShiftSolve> eigs(op, 3, 6, 3.14);
///     eigs.init();
///     eigs.compute(SortRule::LargestMagn);
///     if (eigs.info() == CompInfo::Successful)
///     {
///         Eigen::VectorXd evalues = eigs.eigenvalues();
///         // Will get (4.0, 3.0, 2.0)
///         std::cout << "Eigenvalues found:\n" << evalues << std::endl;
///     }
///
///     return 0;
/// }
/// \endcode
///
template <typename OpType = DenseSymShiftSolve<double>>
class SymEigsShiftSolver : public SymEigsBase<OpType, IdentityBOp>
{
private:
    using Scalar = typename OpType::Scalar;
    using Index = Eigen::Index;
    using Array = Eigen::Array<Scalar, Eigen::Dynamic, 1>;

    using Base = SymEigsBase<OpType, IdentityBOp>;
    using Base::m_nev;
    using Base::m_ritz_val;

    const Scalar m_sigma;

    // First transform back the Ritz values, and then sort
    void sort_ritzpair(SortRule sort_rule) override
    {
        // The eigenvalues we get from the iteration is nu = 1 / (lambda - sigma)
        // So the eigenvalues of the original problem is lambda = 1 / nu + sigma
        m_ritz_val.head(m_nev).array() = Scalar(1) / m_ritz_val.head(m_nev).array() + m_sigma;
        Base::sort_ritzpair(sort_rule);
    }

public:
    ///
    /// Constructor to create a eigen solver object using the shift-and-invert mode.
    ///
    /// \param op     The matrix operation object that implements
    ///               the shift-solve operation of \f$A\f$: calculating
    ///               \f$(A-\sigma I)^{-1}v\f$ for any vector \f$v\f$. Users could either
    ///               create the object from the wrapper class such as DenseSymShiftSolve, or
    ///               define their own that implements all the public members
    ///               as in DenseSymShiftSolve.
    /// \param nev    Number of eigenvalues requested. This should satisfy \f$1\le nev \le n-1\f$,
    ///               where \f$n\f$ is the size of matrix.
    /// \param ncv    Parameter that controls the convergence speed of the algorithm.
    ///               Typically a larger `ncv_` means faster convergence, but it may
    ///               also result in greater memory use and more matrix operations
    ///               in each iteration. This parameter must satisfy \f$nev < ncv \le n\f$,
    ///               and is advised to take \f$ncv \ge 2\cdot nev\f$.
    /// \param sigma  The value of the shift.
    ///
    SymEigsShiftSolver(OpType& op, Index nev, Index ncv, const Scalar& sigma) :
        Base(op, IdentityBOp(), nev, ncv),
        m_sigma(sigma)
    {
        op.set_shift(m_sigma);
    }
};

}  // namespace Spectra

#endif  // SPECTRA_SYM_EIGS_SHIFT_SOLVER_H
