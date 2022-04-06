// Copyright (C) 2016-2022 Yixuan Qiu <yixuan.qiu@cos.name>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SPECTRA_GEN_EIGS_SOLVER_H
#define SPECTRA_GEN_EIGS_SOLVER_H

#include <Eigen/Core>

#include "GenEigsBase.h"
#include "Util/SelectionRule.h"
#include "MatOp/DenseGenMatProd.h"

namespace Spectra {

///
/// \ingroup EigenSolver
///
/// This class implements the eigen solver for general real matrices, i.e.,
/// to solve \f$Ax=\lambda x\f$ for a possibly non-symmetric \f$A\f$ matrix.
///
/// Most of the background information documented in the SymEigsSolver class
/// also applies to the GenEigsSolver class here, except that the eigenvalues
/// and eigenvectors of a general matrix can now be complex-valued.
///
/// \tparam OpType  The name of the matrix operation class. Users could either
///                 use the wrapper classes such as DenseGenMatProd and
///                 SparseGenMatProd, or define their own that implements the type
///                 definition `Scalar` and all the public member functions as in
///                 DenseGenMatProd.
///
/// An example that illustrates the usage of GenEigsSolver is give below:
///
/// \code{.cpp}
/// #include <Eigen/Core>
/// #include <Spectra/GenEigsSolver.h>
/// // <Spectra/MatOp/DenseGenMatProd.h> is implicitly included
/// #include <iostream>
///
/// using namespace Spectra;
///
/// int main()
/// {
///     // We are going to calculate the eigenvalues of M
///     Eigen::MatrixXd M = Eigen::MatrixXd::Random(10, 10);
///
///     // Construct matrix operation object using the wrapper class
///     DenseGenMatProd<double> op(M);
///
///     // Construct eigen solver object, requesting the largest
///     // (in magnitude, or norm) three eigenvalues
///     GenEigsSolver<DenseGenMatProd<double>> eigs(op, 3, 6);
///
///     // Initialize and compute
///     eigs.init();
///     int nconv = eigs.compute(SortRule::LargestMagn);
///
///     // Retrieve results
///     Eigen::VectorXcd evalues;
///     if (eigs.info() == CompInfo::Successful)
///         evalues = eigs.eigenvalues();
///
///     std::cout << "Eigenvalues found:\n" << evalues << std::endl;
///
///     return 0;
/// }
/// \endcode
///
/// And also an example for sparse matrices:
///
/// \code{.cpp}
/// #include <Eigen/Core>
/// #include <Eigen/SparseCore>
/// #include <Spectra/GenEigsSolver.h>
/// #include <Spectra/MatOp/SparseGenMatProd.h>
/// #include <iostream>
///
/// using namespace Spectra;
///
/// int main()
/// {
///     // A band matrix with 1 on the main diagonal, 2 on the below-main subdiagonal,
///     // and 3 on the above-main subdiagonal
///     const int n = 10;
///     Eigen::SparseMatrix<double> M(n, n);
///     M.reserve(Eigen::VectorXi::Constant(n, 3));
///     for (int i = 0; i < n; i++)
///     {
///         M.insert(i, i) = 1.0;
///         if (i > 0)
///             M.insert(i - 1, i) = 3.0;
///         if (i < n - 1)
///             M.insert(i + 1, i) = 2.0;
///     }
///
///     // Construct matrix operation object using the wrapper class SparseGenMatProd
///     SparseGenMatProd<double> op(M);
///
///     // Construct eigen solver object, requesting the largest three eigenvalues
///     GenEigsSolver<SparseGenMatProd<double>> eigs(op, 3, 6);
///
///     // Initialize and compute
///     eigs.init();
///     int nconv = eigs.compute(SortRule::LargestMagn);
///
///     // Retrieve results
///     Eigen::VectorXcd evalues;
///     if (eigs.info() == CompInfo::Successful)
///         evalues = eigs.eigenvalues();
///
///     std::cout << "Eigenvalues found:\n" << evalues << std::endl;
///
///     return 0;
/// }
/// \endcode
template <typename OpType = DenseGenMatProd<double>>
class GenEigsSolver : public GenEigsBase<OpType, IdentityBOp>
{
private:
    using Index = Eigen::Index;

public:
    ///
    /// Constructor to create a solver object.
    ///
    /// \param op   The matrix operation object that implements
    ///             the matrix-vector multiplication operation of \f$A\f$:
    ///             calculating \f$Av\f$ for any vector \f$v\f$. Users could either
    ///             create the object from the wrapper class such as DenseGenMatProd, or
    ///             define their own that implements all the public members
    ///             as in DenseGenMatProd.
    /// \param nev  Number of eigenvalues requested. This should satisfy \f$1\le nev \le n-2\f$,
    ///             where \f$n\f$ is the size of matrix.
    /// \param ncv  Parameter that controls the convergence speed of the algorithm.
    ///             Typically a larger `ncv` means faster convergence, but it may
    ///             also result in greater memory use and more matrix operations
    ///             in each iteration. This parameter must satisfy \f$nev+2 \le ncv \le n\f$,
    ///             and is advised to take \f$ncv \ge 2\cdot nev + 1\f$.
    ///
    GenEigsSolver(OpType& op, Index nev, Index ncv) :
        GenEigsBase<OpType, IdentityBOp>(op, IdentityBOp(), nev, ncv)
    {}
};

}  // namespace Spectra

#endif  // SPECTRA_GEN_EIGS_SOLVER_H
