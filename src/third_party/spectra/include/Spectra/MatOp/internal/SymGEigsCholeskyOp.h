// Copyright (C) 2016-2019 Yixuan Qiu <yixuan.qiu@cos.name>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SYM_GEIGS_CHOLESKY_OP_H
#define SYM_GEIGS_CHOLESKY_OP_H

#include <Eigen/Core>
#include "../DenseSymMatProd.h"
#include "../DenseCholesky.h"

namespace Spectra {


///
/// \ingroup Operators
///
/// This class defines the matrix operation for generalized eigen solver in the
/// Cholesky decomposition mode. It calculates \f$y=L^{-1}A(L')^{-1}x\f$ for any
/// vector \f$x\f$, where \f$L\f$ is the Cholesky decomposition of \f$B\f$.
/// This class is intended for internal use.
///
template < typename Scalar = double,
           typename OpType = DenseSymMatProd<double>,
           typename BOpType = DenseCholesky<double> >
class SymGEigsCholeskyOp
{
private:
    typedef Eigen::Index Index;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Matrix;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;

    OpType&  m_op;
    BOpType& m_Bop;
    Vector   m_cache;  // temporary working space

public:
    ///
    /// Constructor to create the matrix operation object.
    ///
    /// \param op   Pointer to the \f$A\f$ matrix operation object.
    /// \param Bop  Pointer to the \f$B\f$ matrix operation object.
    ///
    SymGEigsCholeskyOp(OpType& op, BOpType& Bop) :
        m_op(op), m_Bop(Bop), m_cache(op.rows())
    {}

    ///
    /// Return the number of rows of the underlying matrix.
    ///
    Index rows() const { return m_Bop.rows(); }
    ///
    /// Return the number of columns of the underlying matrix.
    ///
    Index cols() const { return m_Bop.rows(); }

    ///
    /// Perform the matrix operation \f$y=L^{-1}A(L')^{-1}x\f$.
    ///
    /// \param x_in  Pointer to the \f$x\f$ vector.
    /// \param y_out Pointer to the \f$y\f$ vector.
    ///
    // y_out = inv(L) * A * inv(L') * x_in
    void perform_op(const Scalar* x_in, Scalar* y_out)
    {
        m_Bop.upper_triangular_solve(x_in, y_out);
        m_op.perform_op(y_out, m_cache.data());
        m_Bop.lower_triangular_solve(m_cache.data(), y_out);
    }
};


} // namespace Spectra

#endif // SYM_GEIGS_CHOLESKY_OP_H
