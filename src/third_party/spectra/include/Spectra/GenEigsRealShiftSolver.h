// Copyright (C) 2016-2019 Yixuan Qiu <yixuan.qiu@cos.name>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef GEN_EIGS_REAL_SHIFT_SOLVER_H
#define GEN_EIGS_REAL_SHIFT_SOLVER_H

#include <Eigen/Core>

#include "GenEigsBase.h"
#include "Util/SelectionRule.h"
#include "MatOp/DenseGenRealShiftSolve.h"

namespace Spectra {


///
/// \ingroup EigenSolver
///
/// This class implements the eigen solver for general real matrices with
/// a real shift value in the **shift-and-invert mode**. The background
/// knowledge of the shift-and-invert mode can be found in the documentation
/// of the SymEigsShiftSolver class.
///
/// \tparam Scalar        The element type of the matrix.
///                       Currently supported types are `float`, `double` and `long double`.
/// \tparam SelectionRule An enumeration value indicating the selection rule of
///                       the shifted-and-inverted eigenvalues.
///                       The full list of enumeration values can be found in
///                       \ref Enumerations.
/// \tparam OpType        The name of the matrix operation class. Users could either
///                       use the wrapper classes such as DenseGenRealShiftSolve and
///                       SparseGenRealShiftSolve, or define their
///                       own that implements all the public member functions as in
///                       DenseGenRealShiftSolve.
///
template <typename Scalar = double,
          int SelectionRule = LARGEST_MAGN,
          typename OpType = DenseGenRealShiftSolve<double> >
class GenEigsRealShiftSolver: public GenEigsBase<Scalar, SelectionRule, OpType, IdentityBOp>
{
private:
    typedef Eigen::Index Index;
    typedef std::complex<Scalar> Complex;
    typedef Eigen::Array<Complex, Eigen::Dynamic, 1> ComplexArray;

    const Scalar m_sigma;

    // First transform back the Ritz values, and then sort
    void sort_ritzpair(int sort_rule)
    {
        // The eigenvalues we get from the iteration is nu = 1 / (lambda - sigma)
        // So the eigenvalues of the original problem is lambda = 1 / nu + sigma
        ComplexArray ritz_val_org = Scalar(1.0) / this->m_ritz_val.head(this->m_nev).array() + m_sigma;
        this->m_ritz_val.head(this->m_nev) = ritz_val_org;
        GenEigsBase<Scalar, SelectionRule, OpType, IdentityBOp>::sort_ritzpair(sort_rule);
    }
public:
    ///
    /// Constructor to create a eigen solver object using the shift-and-invert mode.
    ///
    /// \param op     Pointer to the matrix operation object. This class should implement
    ///               the shift-solve operation of \f$A\f$: calculating
    ///               \f$(A-\sigma I)^{-1}v\f$ for any vector \f$v\f$. Users could either
    ///               create the object from the wrapper class such as DenseGenRealShiftSolve, or
    ///               define their own that implements all the public member functions
    ///               as in DenseGenRealShiftSolve.
    /// \param nev    Number of eigenvalues requested. This should satisfy \f$1\le nev \le n-2\f$,
    ///               where \f$n\f$ is the size of matrix.
    /// \param ncv    Parameter that controls the convergence speed of the algorithm.
    ///               Typically a larger `ncv` means faster convergence, but it may
    ///               also result in greater memory use and more matrix operations
    ///               in each iteration. This parameter must satisfy \f$nev+2 \le ncv \le n\f$,
    ///               and is advised to take \f$ncv \ge 2\cdot nev + 1\f$.
    /// \param sigma  The real-valued shift.
    ///
    GenEigsRealShiftSolver(OpType* op, Index nev, Index ncv, Scalar sigma) :
        GenEigsBase<Scalar, SelectionRule, OpType, IdentityBOp>(op, NULL, nev, ncv),
        m_sigma(sigma)
    {
        this->m_op->set_shift(m_sigma);
    }
};


} // namespace Spectra

#endif // GEN_EIGS_REAL_SHIFT_SOLVER_H
