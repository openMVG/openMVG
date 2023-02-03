// Copyright (C) 2016-2022 Yixuan Qiu <yixuan.qiu@cos.name>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SPECTRA_GEN_EIGS_REAL_SHIFT_SOLVER_H
#define SPECTRA_GEN_EIGS_REAL_SHIFT_SOLVER_H

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
/// \tparam OpType  The name of the matrix operation class. Users could either
///                 use the wrapper classes such as DenseGenRealShiftSolve and
///                 SparseGenRealShiftSolve, or define their own that implements the type
///                 definition `Scalar` and all the public member functions as in
///                 DenseGenRealShiftSolve.
///
template <typename OpType = DenseGenRealShiftSolve<double>>
class GenEigsRealShiftSolver : public GenEigsBase<OpType, IdentityBOp>
{
private:
    using Scalar = typename OpType::Scalar;
    using Index = Eigen::Index;
    using Complex = std::complex<Scalar>;
    using ComplexArray = Eigen::Array<Complex, Eigen::Dynamic, 1>;

    using Base = GenEigsBase<OpType, IdentityBOp>;
    using Base::m_nev;
    using Base::m_ritz_val;

    const Scalar m_sigma;

    // First transform back the Ritz values, and then sort
    void sort_ritzpair(SortRule sort_rule) override
    {
        // The eigenvalues we get from the iteration is nu = 1 / (lambda - sigma)
        // So the eigenvalues of the original problem is lambda = 1 / nu + sigma
        m_ritz_val.head(m_nev) = Scalar(1) / m_ritz_val.head(m_nev).array() + m_sigma;
        Base::sort_ritzpair(sort_rule);
    }

public:
    ///
    /// Constructor to create a eigen solver object using the shift-and-invert mode.
    ///
    /// \param op     The matrix operation object that implements
    ///               the shift-solve operation of \f$A\f$: calculating
    ///               \f$(A-\sigma I)^{-1}v\f$ for any vector \f$v\f$. Users could either
    ///               create the object from the wrapper class such as DenseGenRealShiftSolve, or
    ///               define their own that implements all the public members
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
    GenEigsRealShiftSolver(OpType& op, Index nev, Index ncv, const Scalar& sigma) :
        Base(op, IdentityBOp(), nev, ncv),
        m_sigma(sigma)
    {
        op.set_shift(m_sigma);
    }
};

}  // namespace Spectra

#endif  // SPECTRA_GEN_EIGS_REAL_SHIFT_SOLVER_H
