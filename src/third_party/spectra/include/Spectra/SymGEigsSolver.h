// Copyright (C) 2016-2022 Yixuan Qiu <yixuan.qiu@cos.name>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SPECTRA_SYM_GEIGS_SOLVER_H
#define SPECTRA_SYM_GEIGS_SOLVER_H

#include "SymEigsBase.h"
#include "Util/GEigsMode.h"
#include "MatOp/internal/SymGEigsCholeskyOp.h"
#include "MatOp/internal/SymGEigsRegInvOp.h"

namespace Spectra {

///
/// \defgroup GEigenSolver Generalized Eigen Solvers
///
/// Generalized eigen solvers for different types of problems.
///

///
/// \ingroup GEigenSolver
///
/// This class implements the generalized eigen solver for real symmetric
/// matrices, i.e., to solve \f$Ax=\lambda Bx\f$ where \f$A\f$ is symmetric and
/// \f$B\f$ is positive definite.
///
/// There are two modes of this solver, specified by the template parameter `Mode`.
/// See the pages for the specialized classes for details.
/// - The Cholesky mode assumes that \f$B\f$ can be factorized using Cholesky
///   decomposition, which is the preferred mode when the decomposition is
///   available. (This can be easily done in Eigen using the dense or sparse
///   Cholesky solver.)
///   See \ref SymGEigsSolver<OpType, BOpType, GEigsMode::Cholesky> "SymGEigsSolver (Cholesky mode)" for more details.
/// - The regular inverse mode requires the matrix-vector product \f$Bv\f$ and the
///   linear equation solving operation \f$B^{-1}v\f$. This mode should only be
///   used when the Cholesky decomposition of \f$B\f$ is hard to implement, or
///   when computing \f$B^{-1}v\f$ is much faster than the Cholesky decomposition.
///   See \ref SymGEigsSolver<OpType, BOpType, GEigsMode::RegularInverse> "SymGEigsSolver (Regular inverse mode)" for more details.

// Empty class template
template <typename OpType, typename BOpType, GEigsMode Mode>
class SymGEigsSolver
{};

///
/// \ingroup GEigenSolver
///
/// This class implements the generalized eigen solver for real symmetric
/// matrices using Cholesky decomposition, i.e., to solve \f$Ax=\lambda Bx\f$
/// where \f$A\f$ is symmetric and \f$B\f$ is positive definite with the Cholesky
/// decomposition \f$B=LL'\f$.
///
/// This solver requires two matrix operation objects: one for \f$A\f$ that implements
/// the matrix multiplication \f$Av\f$, and one for \f$B\f$ that implements the lower
/// and upper triangular solving \f$L^{-1}v\f$ and \f$(L')^{-1}v\f$.
///
/// If \f$A\f$ and \f$B\f$ are stored as Eigen matrices, then the first operation
/// can be created using the DenseSymMatProd or SparseSymMatProd classes, and
/// the second operation can be created using the DenseCholesky or SparseCholesky
/// classes. If the users need to define their own operation classes, then they
/// should implement all the public member functions as in those built-in classes.
///
/// \tparam OpType   The name of the matrix operation class for \f$A\f$. Users could either
///                  use the wrapper classes such as DenseSymMatProd and
///                  SparseSymMatProd, or define their own that implements the type
///                  definition `Scalar` and all the public member functions as in
///                  DenseSymMatProd.
/// \tparam BOpType  The name of the matrix operation class for \f$B\f$. Users could either
///                  use the wrapper classes such as DenseCholesky and
///                  SparseCholesky, or define their own that implements all the
///                  public member functions as in DenseCholesky.
/// \tparam Mode     Mode of the generalized eigen solver. In this solver
///                  it is Spectra::GEigsMode::Cholesky.
///
/// Below is an example that demonstrates the usage of this class.
///
/// \code{.cpp}
/// #include <Eigen/Core>
/// #include <Eigen/SparseCore>
/// #include <Eigen/Eigenvalues>
/// #include <Spectra/SymGEigsSolver.h>
/// #include <Spectra/MatOp/DenseSymMatProd.h>
/// #include <Spectra/MatOp/SparseCholesky.h>
/// #include <iostream>
///
/// using namespace Spectra;
///
/// int main()
/// {
///     // We are going to solve the generalized eigenvalue problem A * x = lambda * B * x
///     const int n = 100;
///
///     // Define the A matrix
///     Eigen::MatrixXd M = Eigen::MatrixXd::Random(n, n);
///     Eigen::MatrixXd A = M + M.transpose();
///
///     // Define the B matrix, a band matrix with 2 on the diagonal and 1 on the subdiagonals
///     Eigen::SparseMatrix<double> B(n, n);
///     B.reserve(Eigen::VectorXi::Constant(n, 3));
///     for (int i = 0; i < n; i++)
///     {
///         B.insert(i, i) = 2.0;
///         if (i > 0)
///             B.insert(i - 1, i) = 1.0;
///         if (i < n - 1)
///             B.insert(i + 1, i) = 1.0;
///     }
///
///     // Construct matrix operation objects using the wrapper classes
///     DenseSymMatProd<double> op(A);
///     SparseCholesky<double>  Bop(B);
///
///     // Construct generalized eigen solver object, requesting the largest three generalized eigenvalues
///     SymGEigsSolver<DenseSymMatProd<double>, SparseCholesky<double>, GEigsMode::Cholesky>
///         geigs(op, Bop, 3, 6);
///
///     // Initialize and compute
///     geigs.init();
///     int nconv = geigs.compute(SortRule::LargestAlge);
///
///     // Retrieve results
///     Eigen::VectorXd evalues;
///     Eigen::MatrixXd evecs;
///     if (geigs.info() == CompInfo::Successful)
///     {
///         evalues = geigs.eigenvalues();
///         evecs = geigs.eigenvectors();
///     }
///
///     std::cout << "Generalized eigenvalues found:\n" << evalues << std::endl;
///     std::cout << "Generalized eigenvectors found:\n" << evecs.topRows(10) << std::endl;
///
///     // Verify results using the generalized eigen solver in Eigen
///     Eigen::MatrixXd Bdense = B;
///     Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> es(A, Bdense);
///
///     std::cout << "Generalized eigenvalues:\n" << es.eigenvalues().tail(3) << std::endl;
///     std::cout << "Generalized eigenvectors:\n" << es.eigenvectors().rightCols(3).topRows(10) << std::endl;
///
///     return 0;
/// }
/// \endcode

// Partial specialization for mode = GEigsMode::Cholesky
template <typename OpType, typename BOpType>
class SymGEigsSolver<OpType, BOpType, GEigsMode::Cholesky> :
    public SymEigsBase<SymGEigsCholeskyOp<OpType, BOpType>, IdentityBOp>
{
private:
    using Scalar = typename OpType::Scalar;
    using Index = Eigen::Index;
    using Matrix = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
    using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

    using ModeMatOp = SymGEigsCholeskyOp<OpType, BOpType>;
    using Base = SymEigsBase<ModeMatOp, IdentityBOp>;

    const BOpType& m_Bop;

public:
    ///
    /// Constructor to create a solver object.
    ///
    /// \param op   The \f$A\f$ matrix operation object that implements the matrix-vector
    ///             multiplication operation of \f$A\f$:
    ///             calculating \f$Av\f$ for any vector \f$v\f$. Users could either
    ///             create the object from the wrapper classes such as DenseSymMatProd, or
    ///             define their own that implements all the public members
    ///             as in DenseSymMatProd.
    /// \param Bop  The \f$B\f$ matrix operation object that represents a Cholesky decomposition of \f$B\f$.
    ///             It should implement the lower and upper triangular solving operations:
    ///             calculating \f$L^{-1}v\f$ and \f$(L')^{-1}v\f$ for any vector
    ///             \f$v\f$, where \f$LL'=B\f$. Users could either
    ///             create the object from the wrapper classes such as DenseCholesky, or
    ///             define their own that implements all the public member functions
    ///             as in DenseCholesky. \f$B\f$ needs to be positive definite.
    /// \param nev  Number of eigenvalues requested. This should satisfy \f$1\le nev \le n-1\f$,
    ///             where \f$n\f$ is the size of matrix.
    /// \param ncv  Parameter that controls the convergence speed of the algorithm.
    ///             Typically a larger `ncv` means faster convergence, but it may
    ///             also result in greater memory use and more matrix operations
    ///             in each iteration. This parameter must satisfy \f$nev < ncv \le n\f$,
    ///             and is advised to take \f$ncv \ge 2\cdot nev\f$.
    ///
    SymGEigsSolver(OpType& op, BOpType& Bop, Index nev, Index ncv) :
        Base(ModeMatOp(op, Bop), IdentityBOp(), nev, ncv),
        m_Bop(Bop)
    {}

    /// \cond

    Matrix eigenvectors(Index nvec) const override
    {
        Matrix res = Base::eigenvectors(nvec);
        Vector tmp(res.rows());
        const Index nconv = res.cols();
        for (Index i = 0; i < nconv; i++)
        {
            m_Bop.upper_triangular_solve(&res(0, i), tmp.data());
            res.col(i).noalias() = tmp;
        }

        return res;
    }

    Matrix eigenvectors() const override
    {
        return SymGEigsSolver<OpType, BOpType, GEigsMode::Cholesky>::eigenvectors(this->m_nev);
    }

    /// \endcond
};

///
/// \ingroup GEigenSolver
///
/// This class implements the generalized eigen solver for real symmetric
/// matrices in the regular inverse mode, i.e., to solve \f$Ax=\lambda Bx\f$
/// where \f$A\f$ is symmetric, and \f$B\f$ is positive definite with the operations
/// defined below.
///
/// This solver requires two matrix operation objects: one for \f$A\f$ that implements
/// the matrix multiplication \f$Av\f$, and one for \f$B\f$ that implements the
/// matrix-vector product \f$Bv\f$ and the linear equation solving operation \f$B^{-1}v\f$.
///
/// If \f$A\f$ and \f$B\f$ are stored as Eigen matrices, then the first operation
/// can be created using the DenseSymMatProd or SparseSymMatProd classes, and
/// the second operation can be created using the SparseRegularInverse class. There is no
/// wrapper class for a dense \f$B\f$ matrix since in this case the Cholesky mode
/// is always preferred. If the users need to define their own operation classes, then they
/// should implement all the public member functions as in those built-in classes.
///
/// \tparam OpType   The name of the matrix operation class for \f$A\f$. Users could either
///                  use the wrapper classes such as DenseSymMatProd and
///                  SparseSymMatProd, or define their own that implements the type
///                  definition `Scalar` and all the public member functions as in
///                  DenseSymMatProd.
/// \tparam BOpType  The name of the matrix operation class for \f$B\f$. Users could either
///                  use the wrapper class SparseRegularInverse, or define their
///                  own that implements all the public member functions as in
///                  SparseRegularInverse.
/// \tparam Mode     Mode of the generalized eigen solver. In this solver
///                  it is Spectra::GEigsMode::RegularInverse.
///

// Partial specialization for mode = GEigsMode::RegularInverse
template <typename OpType, typename BOpType>
class SymGEigsSolver<OpType, BOpType, GEigsMode::RegularInverse> :
    public SymEigsBase<SymGEigsRegInvOp<OpType, BOpType>, BOpType>
{
private:
    using Index = Eigen::Index;

    using ModeMatOp = SymGEigsRegInvOp<OpType, BOpType>;
    using Base = SymEigsBase<ModeMatOp, BOpType>;

public:
    ///
    /// Constructor to create a solver object.
    ///
    /// \param op   The \f$A\f$ matrix operation object that implements the matrix-vector
    ///             multiplication operation of \f$A\f$:
    ///             calculating \f$Av\f$ for any vector \f$v\f$. Users could either
    ///             create the object from the wrapper classes such as DenseSymMatProd, or
    ///             define their own that implements all the public members
    ///             as in DenseSymMatProd.
    /// \param Bop  The \f$B\f$ matrix operation object that implements the multiplication operation
    ///             \f$Bv\f$ and the linear equation solving operation \f$B^{-1}v\f$ for any vector \f$v\f$.
    ///             Users could either create the object from the wrapper class SparseRegularInverse, or
    ///             define their own that implements all the public member functions
    ///             as in SparseRegularInverse. \f$B\f$ needs to be positive definite.
    /// \param nev  Number of eigenvalues requested. This should satisfy \f$1\le nev \le n-1\f$,
    ///             where \f$n\f$ is the size of matrix.
    /// \param ncv  Parameter that controls the convergence speed of the algorithm.
    ///             Typically a larger `ncv` means faster convergence, but it may
    ///             also result in greater memory use and more matrix operations
    ///             in each iteration. This parameter must satisfy \f$nev < ncv \le n\f$,
    ///             and is advised to take \f$ncv \ge 2\cdot nev\f$.
    ///
    SymGEigsSolver(OpType& op, BOpType& Bop, Index nev, Index ncv) :
        Base(ModeMatOp(op, Bop), Bop, nev, ncv)
    {}
};

}  // namespace Spectra

#endif  // SPECTRA_SYM_GEIGS_SOLVER_H
