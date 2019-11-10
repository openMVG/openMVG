// Copyright (C) 2016-2019 Yixuan Qiu <yixuan.qiu@cos.name>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SYM_EIGS_SOLVER_H
#define SYM_EIGS_SOLVER_H

#include <Eigen/Core>

#include "SymEigsBase.h"
#include "Util/SelectionRule.h"
#include "MatOp/DenseSymMatProd.h"

namespace Spectra {


///
/// \ingroup EigenSolver
///
/// This class implements the eigen solver for real symmetric matrices, i.e.,
/// to solve \f$Ax=\lambda x\f$ where \f$A\f$ is symmetric.
///
/// **Spectra** is designed to calculate a specified number (\f$k\f$)
/// of eigenvalues of a large square matrix (\f$A\f$). Usually \f$k\f$ is much
/// less than the size of the matrix (\f$n\f$), so that only a few eigenvalues
/// and eigenvectors are computed.
///
/// Rather than providing the whole \f$A\f$ matrix, the algorithm only requires
/// the matrix-vector multiplication operation of \f$A\f$. Therefore, users of
/// this solver need to supply a class that computes the result of \f$Av\f$
/// for any given vector \f$v\f$. The name of this class should be given to
/// the template parameter `OpType`, and instance of this class passed to
/// the constructor of SymEigsSolver.
///
/// If the matrix \f$A\f$ is already stored as a matrix object in **Eigen**,
/// for example `Eigen::MatrixXd`, then there is an easy way to construct such
/// matrix operation class, by using the built-in wrapper class DenseSymMatProd
/// which wraps an existing matrix object in **Eigen**. This is also the
/// default template parameter for SymEigsSolver. For sparse matrices, the
/// wrapper class SparseSymMatProd can be used similarly.
///
/// If the users need to define their own matrix-vector multiplication operation
/// class, it should implement all the public member functions as in DenseSymMatProd.
///
/// \tparam Scalar        The element type of the matrix.
///                       Currently supported types are `float`, `double` and `long double`.
/// \tparam SelectionRule An enumeration value indicating the selection rule of
///                       the requested eigenvalues, for example `LARGEST_MAGN`
///                       to retrieve eigenvalues with the largest magnitude.
///                       The full list of enumeration values can be found in
///                       \ref Enumerations.
/// \tparam OpType        The name of the matrix operation class. Users could either
///                       use the wrapper classes such as DenseSymMatProd and
///                       SparseSymMatProd, or define their
///                       own that implements all the public member functions as in
///                       DenseSymMatProd.
///
/// Below is an example that demonstrates the usage of this class.
///
/// \code{.cpp}
/// #include <Eigen/Core>
/// #include <Spectra/SymEigsSolver.h>
/// // <Spectra/MatOp/DenseSymMatProd.h> is implicitly included
/// #include <iostream>
///
/// using namespace Spectra;
///
/// int main()
/// {
///     // We are going to calculate the eigenvalues of M
///     Eigen::MatrixXd A = Eigen::MatrixXd::Random(10, 10);
///     Eigen::MatrixXd M = A + A.transpose();
///
///     // Construct matrix operation object using the wrapper class DenseSymMatProd
///     DenseSymMatProd<double> op(M);
///
///     // Construct eigen solver object, requesting the largest three eigenvalues
///     SymEigsSolver< double, LARGEST_ALGE, DenseSymMatProd<double> > eigs(&op, 3, 6);
///
///     // Initialize and compute
///     eigs.init();
///     int nconv = eigs.compute();
///
///     // Retrieve results
///     Eigen::VectorXd evalues;
///     if(eigs.info() == SUCCESSFUL)
///         evalues = eigs.eigenvalues();
///
///     std::cout << "Eigenvalues found:\n" << evalues << std::endl;
///
///     return 0;
/// }
/// \endcode
///
/// And here is an example for user-supplied matrix operation class.
///
/// \code{.cpp}
/// #include <Eigen/Core>
/// #include <Spectra/SymEigsSolver.h>
/// #include <iostream>
///
/// using namespace Spectra;
///
/// // M = diag(1, 2, ..., 10)
/// class MyDiagonalTen
/// {
/// public:
///     int rows() { return 10; }
///     int cols() { return 10; }
///     // y_out = M * x_in
///     void perform_op(double *x_in, double *y_out)
///     {
///         for(int i = 0; i < rows(); i++)
///         {
///             y_out[i] = x_in[i] * (i + 1);
///         }
///     }
/// };
///
/// int main()
/// {
///     MyDiagonalTen op;
///     SymEigsSolver<double, LARGEST_ALGE, MyDiagonalTen> eigs(&op, 3, 6);
///     eigs.init();
///     eigs.compute();
///     if(eigs.info() == SUCCESSFUL)
///     {
///         Eigen::VectorXd evalues = eigs.eigenvalues();
///         // Will get (10, 9, 8)
///         std::cout << "Eigenvalues found:\n" << evalues << std::endl;
///     }
///
///     return 0;
/// }
/// \endcode
///
template < typename Scalar = double,
           int SelectionRule = LARGEST_MAGN,
           typename OpType = DenseSymMatProd<double> >
class SymEigsSolver: public SymEigsBase<Scalar, SelectionRule, OpType, IdentityBOp>
{
private:
    typedef Eigen::Index Index;

public:
    ///
    /// Constructor to create a solver object.
    ///
    /// \param op   Pointer to the matrix operation object, which should implement
    ///             the matrix-vector multiplication operation of \f$A\f$:
    ///             calculating \f$Av\f$ for any vector \f$v\f$. Users could either
    ///             create the object from the wrapper class such as DenseSymMatProd, or
    ///             define their own that implements all the public member functions
    ///             as in DenseSymMatProd.
    /// \param nev  Number of eigenvalues requested. This should satisfy \f$1\le nev \le n-1\f$,
    ///             where \f$n\f$ is the size of matrix.
    /// \param ncv  Parameter that controls the convergence speed of the algorithm.
    ///             Typically a larger `ncv` means faster convergence, but it may
    ///             also result in greater memory use and more matrix operations
    ///             in each iteration. This parameter must satisfy \f$nev < ncv \le n\f$,
    ///             and is advised to take \f$ncv \ge 2\cdot nev\f$.
    ///
    SymEigsSolver(OpType* op, Index nev, Index ncv) :
        SymEigsBase<Scalar, SelectionRule, OpType, IdentityBOp>(op, NULL, nev, ncv)
    {}

};


} // namespace Spectra

#endif // SYM_EIGS_SOLVER_H
