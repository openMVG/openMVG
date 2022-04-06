// Copyright (C) 2020-2022 Yixuan Qiu <yixuan.qiu@cos.name>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SPECTRA_SYM_GEIGS_BUCKLING_OP_H
#define SPECTRA_SYM_GEIGS_BUCKLING_OP_H

#include <Eigen/Core>

#include "../SymShiftInvert.h"
#include "../SparseSymMatProd.h"

namespace Spectra {

///
/// \ingroup Operators
///
/// This class defines the matrix operation for generalized eigen solver in the
/// buckling mode. It computes \f$y=(K-\sigma K_G)^{-1}Kx\f$ for any
/// vector \f$x\f$, where \f$K\f$ is positive definite, \f$K_G\f$ is symmetric,
/// and \f$\sigma\f$ is a real shift.
/// This class is intended for internal use.
///
template <typename OpType = SymShiftInvert<double>,
          typename BOpType = SparseSymMatProd<double>>
class SymGEigsBucklingOp
{
public:
    using Scalar = typename OpType::Scalar;

private:
    using Index = Eigen::Index;
    using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

    OpType& m_op;
    const BOpType& m_Bop;
    mutable Vector m_cache;  // temporary working space

public:
    ///
    /// Constructor to create the matrix operation object.
    ///
    /// \param op   The \f$(K-\sigma K_G)^{-1}\f$ matrix operation object.
    /// \param Bop  The \f$K\f$ matrix operation object.
    ///
    SymGEigsBucklingOp(OpType& op, const BOpType& Bop) :
        m_op(op), m_Bop(Bop), m_cache(op.rows())
    {}

    ///
    /// Move constructor.
    ///
    SymGEigsBucklingOp(SymGEigsBucklingOp&& other) :
        m_op(other.m_op), m_Bop(other.m_Bop)
    {
        // We emulate the move constructor for Vector using Vector::swap()
        m_cache.swap(other.m_cache);
    }

    ///
    /// Return the number of rows of the underlying matrix.
    ///
    Index rows() const { return m_op.rows(); }
    ///
    /// Return the number of columns of the underlying matrix.
    ///
    Index cols() const { return m_op.rows(); }

    ///
    /// Set the real shift \f$\sigma\f$.
    ///
    void set_shift(const Scalar& sigma)
    {
        m_op.set_shift(sigma);
    }

    ///
    /// Perform the matrix operation \f$y=(K-\sigma K_G)^{-1}Kx\f$.
    ///
    /// \param x_in  Pointer to the \f$x\f$ vector.
    /// \param y_out Pointer to the \f$y\f$ vector.
    ///
    // y_out = inv(K - sigma * K_G) * K * x_in
    void perform_op(const Scalar* x_in, Scalar* y_out) const
    {
        m_Bop.perform_op(x_in, m_cache.data());
        m_op.perform_op(m_cache.data(), y_out);
    }
};

}  // namespace Spectra

#endif  // SPECTRA_SYM_GEIGS_BUCKLING_OP_H
