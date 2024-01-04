// Copyright (C) 2020-2022 Yixuan Qiu <yixuan.qiu@cos.name>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SPECTRA_SYM_GEIGS_CAYLEY_OP_H
#define SPECTRA_SYM_GEIGS_CAYLEY_OP_H

#include <Eigen/Core>

#include "../SymShiftInvert.h"
#include "../SparseSymMatProd.h"

namespace Spectra {

///
/// \ingroup Operators
///
/// This class defines the matrix operation for generalized eigen solver in the
/// Cayley mode. It computes \f$y=(A-\sigma B)^{-1}(A+\sigma B)x\f$ for any
/// vector \f$x\f$, where \f$A\f$ is a symmetric matrix, \f$B\f$ is positive definite,
/// and \f$\sigma\f$ is a real shift.
/// This class is intended for internal use.
///
template <typename OpType = SymShiftInvert<double>,
          typename BOpType = SparseSymMatProd<double>>
class SymGEigsCayleyOp
{
public:
    using Scalar = typename OpType::Scalar;

private:
    using Index = Eigen::Index;
    using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
    using MapConstVec = Eigen::Map<const Vector>;
    using MapVec = Eigen::Map<Vector>;

    OpType& m_op;
    const BOpType& m_Bop;
    mutable Vector m_cache;  // temporary working space
    Scalar m_sigma;

public:
    ///
    /// Constructor to create the matrix operation object.
    ///
    /// \param op   The \f$(A-\sigma B)^{-1}\f$ matrix operation object.
    /// \param Bop  The \f$B\f$ matrix operation object.
    ///
    SymGEigsCayleyOp(OpType& op, const BOpType& Bop) :
        m_op(op), m_Bop(Bop), m_cache(op.rows())
    {}

    ///
    /// Move constructor.
    ///
    SymGEigsCayleyOp(SymGEigsCayleyOp&& other) :
        m_op(other.m_op), m_Bop(other.m_Bop), m_sigma(other.m_sigma)
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
        m_sigma = sigma;
    }

    ///
    /// Perform the matrix operation \f$y=(A-\sigma B)^{-1}(A+\sigma B)x\f$.
    ///
    /// \param x_in  Pointer to the \f$x\f$ vector.
    /// \param y_out Pointer to the \f$y\f$ vector.
    ///
    // y_out = inv(A - sigma * B) * (A + sigma * B) * x_in
    void perform_op(const Scalar* x_in, Scalar* y_out) const
    {
        //   inv(A - sigma * B) * (A + sigma * B) * x
        // = inv(A - sigma * B) * (A - sigma * B + 2 * sigma * B) * x
        // = x + 2 * sigma * inv(A - sigma * B) * B * x
        m_Bop.perform_op(x_in, m_cache.data());
        m_op.perform_op(m_cache.data(), y_out);
        MapConstVec x(x_in, this->rows());
        MapVec y(y_out, this->rows());
        y.noalias() = x + (Scalar(2) * m_sigma) * y;
    }
};

}  // namespace Spectra

#endif  // SPECTRA_SYM_GEIGS_CAYLEY_OP_H
