// Copyright (C) 2017-2022 Yixuan Qiu <yixuan.qiu@cos.name>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SPECTRA_SYM_GEIGS_REG_INV_OP_H
#define SPECTRA_SYM_GEIGS_REG_INV_OP_H

#include <Eigen/Core>

#include "../SparseSymMatProd.h"
#include "../SparseRegularInverse.h"

namespace Spectra {

///
/// \ingroup Operators
///
/// This class defines the matrix operation for generalized eigen solver in the
/// regular inverse mode. This class is intended for internal use.
///
template <typename OpType = SparseSymMatProd<double>,
          typename BOpType = SparseRegularInverse<double>>
class SymGEigsRegInvOp
{
public:
    using Scalar = typename OpType::Scalar;

private:
    using Index = Eigen::Index;
    using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

    const OpType& m_op;
    const BOpType& m_Bop;
    mutable Vector m_cache;  // temporary working space

public:
    ///
    /// Constructor to create the matrix operation object.
    ///
    /// \param op   The \f$A\f$ matrix operation object.
    /// \param Bop  The \f$B\f$ matrix operation object.
    ///
    SymGEigsRegInvOp(const OpType& op, const BOpType& Bop) :
        m_op(op), m_Bop(Bop), m_cache(op.rows())
    {}

    ///
    /// Move constructor.
    ///
    SymGEigsRegInvOp(SymGEigsRegInvOp&& other) :
        m_op(other.m_op), m_Bop(other.m_Bop)
    {
        // We emulate the move constructor for Vector using Vector::swap()
        m_cache.swap(other.m_cache);
    }

    ///
    /// Return the number of rows of the underlying matrix.
    ///
    Index rows() const { return m_Bop.rows(); }
    ///
    /// Return the number of columns of the underlying matrix.
    ///
    Index cols() const { return m_Bop.rows(); }

    ///
    /// Perform the matrix operation \f$y=B^{-1}Ax\f$.
    ///
    /// \param x_in  Pointer to the \f$x\f$ vector.
    /// \param y_out Pointer to the \f$y\f$ vector.
    ///
    // y_out = inv(B) * A * x_in
    void perform_op(const Scalar* x_in, Scalar* y_out) const
    {
        m_op.perform_op(x_in, m_cache.data());
        m_Bop.solve(m_cache.data(), y_out);
    }
};

}  // namespace Spectra

#endif  // SPECTRA_SYM_GEIGS_REG_INV_OP_H
