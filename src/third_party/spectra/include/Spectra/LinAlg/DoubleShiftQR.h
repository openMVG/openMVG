// Copyright (C) 2016-2022 Yixuan Qiu <yixuan.qiu@cos.name>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SPECTRA_DOUBLE_SHIFT_QR_H
#define SPECTRA_DOUBLE_SHIFT_QR_H

#include <Eigen/Core>
#include <vector>     // std::vector
#include <algorithm>  // std::min, std::fill, std::copy
#include <utility>    // std::swap
#include <cmath>      // std::abs, std::sqrt, std::pow
#include <stdexcept>  // std::invalid_argument, std::logic_error

#include "../Util/TypeTraits.h"

namespace Spectra {

template <typename Scalar = double>
class DoubleShiftQR
{
private:
    using Index = Eigen::Index;
    using Matrix = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
    using Matrix3X = Eigen::Matrix<Scalar, 3, Eigen::Dynamic>;
    using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
    using IntArray = Eigen::Array<unsigned char, Eigen::Dynamic, 1>;

    using GenericMatrix = Eigen::Ref<Matrix>;
    using ConstGenericMatrix = const Eigen::Ref<const Matrix>;

    // A very small value, but 1.0 / m_near_0 does not overflow
    // ~= 1e-307 for the "double" type
    static constexpr Scalar m_near_0 = TypeTraits<Scalar>::min() * Scalar(10);
    // The machine precision, ~= 1e-16 for the "double" type
    static constexpr Scalar m_eps = TypeTraits<Scalar>::epsilon();

    Index m_n;          // Dimension of the matrix
    Matrix m_mat_H;     // A copy of the matrix to be factorized
    Scalar m_shift_s;   // Shift constant
    Scalar m_shift_t;   // Shift constant
    Matrix3X m_ref_u;   // Householder reflectors
    IntArray m_ref_nr;  // How many rows does each reflector affects
                        // 3 - A general reflector
                        // 2 - A Givens rotation
                        // 1 - An identity transformation
    bool m_computed;    // Whether matrix has been factorized

    // Compute sqrt(x1^2 + x2^2 + x3^2) wit high precision
    static Scalar stable_norm3(Scalar x1, Scalar x2, Scalar x3)
    {
        using std::abs;
        using std::sqrt;

        x1 = abs(x1);
        x2 = abs(x2);
        x3 = abs(x3);
        // Make x1 >= {x2, x3}
        if (x1 < x2)
            std::swap(x1, x2);
        if (x1 < x3)
            std::swap(x1, x3);
        // If x1 is too small, return 0
        if (x1 < m_near_0)
            return Scalar(0);

        const Scalar r2 = x2 / x1, r3 = x3 / x1;
        // We choose a cutoff such that cutoff^4 < eps
        // If max(r2, r3) > cutoff, use the standard way; otherwise use Taylor series expansion
        // to avoid an explicit sqrt() call that may lose precision
        const Scalar cutoff = Scalar(0.1) * pow(m_eps, Scalar(0.25));
        Scalar r = r2 * r2 + r3 * r3;
        r = (r2 >= cutoff || r3 >= cutoff) ?
            sqrt(Scalar(1) + r) :
            (Scalar(1) + r * (Scalar(0.5) - Scalar(0.125) * r));  // sqrt(1 + t) ~= 1 + t/2 - t^2/8
        return x1 * r;
    }

    // x[i] <- x[i] / r, r = sqrt(x1^2 + x2^2 + x3^2)
    // Assume |x1| >= {|x2|, |x3|}, x1 != 0
    static void stable_scaling(Scalar& x1, Scalar& x2, Scalar& x3)
    {
        using std::abs;
        using std::pow;
        using std::sqrt;

        const Scalar x1sign = (x1 > Scalar(0)) ? Scalar(1) : Scalar(-1);
        x1 = abs(x1);
        // Use the same method as in stable_norm3()
        const Scalar r2 = x2 / x1, r3 = x3 / x1;
        const Scalar cutoff = Scalar(0.1) * pow(m_eps, Scalar(0.25));
        Scalar r = r2 * r2 + r3 * r3;
        // r = 1/sqrt(1 + r2^2 + r3^2)
        r = (abs(r2) >= cutoff || abs(r3) >= cutoff) ?
            Scalar(1) / sqrt(Scalar(1) + r) :
            (Scalar(1) - r * (Scalar(0.5) - Scalar(0.375) * r));  // 1/sqrt(1 + t) ~= 1 - t * (1/2 - (3/8) * t)
        x1 = x1sign * r;
        x2 = r2 * r;
        x3 = r3 * r;
    }

    void compute_reflector(const Scalar& x1, const Scalar& x2, const Scalar& x3, Index ind)
    {
        using std::abs;

        Scalar* u = &m_ref_u.coeffRef(0, ind);
        unsigned char* nr = m_ref_nr.data();
        const Scalar x2m = abs(x2), x3m = abs(x3);
        // If both x2 and x3 are zero, nr is 1, and we early exit
        if (x2m < m_near_0 && x3m < m_near_0)
        {
            nr[ind] = 1;
            return;
        }

        // In general case the reflector affects 3 rows
        // If x3 is zero, decrease nr by 1
        nr[ind] = (x3m < m_near_0) ? 2 : 3;
        const Scalar x_norm = (x3m < m_near_0) ? Eigen::numext::hypot(x1, x2) : stable_norm3(x1, x2, x3);

        // x1' = x1 - rho * ||x||
        // rho = -sign(x1), if x1 == 0, we choose rho = 1
        const Scalar rho = (x1 <= Scalar(0)) - (x1 > Scalar(0));
        const Scalar x1_new = x1 - rho * x_norm, x1m = abs(x1_new);
        // Copy x to u
        u[0] = x1_new;
        u[1] = x2;
        u[2] = x3;
        if (x1m >= x2m && x1m >= x3m)
        {
            stable_scaling(u[0], u[1], u[2]);
        }
        else if (x2m >= x1m && x2m >= x3m)
        {
            stable_scaling(u[1], u[0], u[2]);
        }
        else
        {
            stable_scaling(u[2], u[0], u[1]);
        }
    }

    void compute_reflector(const Scalar* x, Index ind)
    {
        compute_reflector(x[0], x[1], x[2], ind);
    }

    // Update the block X = H(il:iu, il:iu)
    void update_block(Index il, Index iu)
    {
        // Block size
        const Index bsize = iu - il + 1;

        // If block size == 1, there is no need to apply reflectors
        if (bsize == 1)
        {
            m_ref_nr.coeffRef(il) = 1;
            return;
        }

        const Scalar x00 = m_mat_H.coeff(il, il),
                     x01 = m_mat_H.coeff(il, il + 1),
                     x10 = m_mat_H.coeff(il + 1, il),
                     x11 = m_mat_H.coeff(il + 1, il + 1);
        // m00 = x00 * (x00 - s) + x01 * x10 + t
        const Scalar m00 = x00 * (x00 - m_shift_s) + x01 * x10 + m_shift_t;
        // m10 = x10 * (x00 + x11 - s)
        const Scalar m10 = x10 * (x00 + x11 - m_shift_s);

        // For block size == 2, do a Givens rotation on M = X * X - s * X + t * I
        if (bsize == 2)
        {
            // This causes nr=2
            compute_reflector(m00, m10, 0, il);
            // Apply the reflector to X
            apply_PX(m_mat_H.block(il, il, 2, m_n - il), m_n, il);
            apply_XP(m_mat_H.block(0, il, il + 2, 2), m_n, il);

            m_ref_nr.coeffRef(il + 1) = 1;
            return;
        }

        // For block size >=3, use the regular strategy
        // m20 = x21 * x10
        const Scalar m20 = m_mat_H.coeff(il + 2, il + 1) * m_mat_H.coeff(il + 1, il);
        compute_reflector(m00, m10, m20, il);

        // Apply the first reflector
        apply_PX(m_mat_H.block(il, il, 3, m_n - il), m_n, il);
        apply_XP(m_mat_H.block(0, il, il + (std::min)(bsize, Index(4)), 3), m_n, il);

        // Calculate the following reflectors
        // If entering this loop, block size is at least 4.
        for (Index i = 1; i < bsize - 2; i++)
        {
            compute_reflector(&m_mat_H.coeffRef(il + i, il + i - 1), il + i);
            // Apply the reflector to X
            apply_PX(m_mat_H.block(il + i, il + i - 1, 3, m_n - il - i + 1), m_n, il + i);
            apply_XP(m_mat_H.block(0, il + i, il + (std::min)(bsize, Index(i + 4)), 3), m_n, il + i);
        }

        // The last reflector
        // This causes nr=2
        compute_reflector(m_mat_H.coeff(iu - 1, iu - 2), m_mat_H.coeff(iu, iu - 2), 0, iu - 1);
        // Apply the reflector to X
        apply_PX(m_mat_H.block(iu - 1, iu - 2, 2, m_n - iu + 2), m_n, iu - 1);
        apply_XP(m_mat_H.block(0, iu - 1, il + bsize, 2), m_n, iu - 1);

        m_ref_nr.coeffRef(iu) = 1;
    }

    // P = I - 2 * u * u' = P'
    // PX = X - 2 * u * (u'X)
    void apply_PX(GenericMatrix X, Index stride, Index u_ind) const
    {
        const Index nr = m_ref_nr.coeff(u_ind);
        if (nr == 1)
            return;

        const Scalar u0 = m_ref_u.coeff(0, u_ind), u1 = m_ref_u.coeff(1, u_ind);
        const Scalar u0_2 = Scalar(2) * u0, u1_2 = Scalar(2) * u1;

        const Index nrow = X.rows();
        const Index ncol = X.cols();

        Scalar* xptr = X.data();
        if (nr == 2 || nrow == 2)
        {
            for (Index i = 0; i < ncol; i++, xptr += stride)
            {
                const Scalar tmp = u0_2 * xptr[0] + u1_2 * xptr[1];
                xptr[0] -= tmp * u0;
                xptr[1] -= tmp * u1;
            }
        }
        else
        {
            const Scalar u2 = m_ref_u.coeff(2, u_ind);
            const Scalar u2_2 = Scalar(2) * u2;
            for (Index i = 0; i < ncol; i++, xptr += stride)
            {
                const Scalar tmp = u0_2 * xptr[0] + u1_2 * xptr[1] + u2_2 * xptr[2];
                xptr[0] -= tmp * u0;
                xptr[1] -= tmp * u1;
                xptr[2] -= tmp * u2;
            }
        }
    }

    // x is a pointer to a vector
    // Px = x - 2 * dot(x, u) * u
    void apply_PX(Scalar* x, Index u_ind) const
    {
        const Index nr = m_ref_nr.coeff(u_ind);
        if (nr == 1)
            return;

        const Scalar u0 = m_ref_u.coeff(0, u_ind),
                     u1 = m_ref_u.coeff(1, u_ind),
                     u2 = m_ref_u.coeff(2, u_ind);

        // When the reflector only contains two elements, u2 has been set to zero
        const bool nr_is_2 = (nr == 2);
        const Scalar dot2 = Scalar(2) * (x[0] * u0 + x[1] * u1 + (nr_is_2 ? 0 : (x[2] * u2)));
        x[0] -= dot2 * u0;
        x[1] -= dot2 * u1;
        if (!nr_is_2)
            x[2] -= dot2 * u2;
    }

    // XP = X - 2 * (X * u) * u'
    void apply_XP(GenericMatrix X, Index stride, Index u_ind) const
    {
        const Index nr = m_ref_nr.coeff(u_ind);
        if (nr == 1)
            return;

        const Scalar u0 = m_ref_u.coeff(0, u_ind), u1 = m_ref_u.coeff(1, u_ind);
        const Scalar u0_2 = Scalar(2) * u0, u1_2 = Scalar(2) * u1;

        const int nrow = X.rows();
        const int ncol = X.cols();
        Scalar *X0 = X.data(), *X1 = X0 + stride;  // X0 => X.col(0), X1 => X.col(1)

        if (nr == 2 || ncol == 2)
        {
            // tmp = 2 * u0 * X0 + 2 * u1 * X1
            // X0 => X0 - u0 * tmp
            // X1 => X1 - u1 * tmp
            for (Index i = 0; i < nrow; i++)
            {
                const Scalar tmp = u0_2 * X0[i] + u1_2 * X1[i];
                X0[i] -= tmp * u0;
                X1[i] -= tmp * u1;
            }
        }
        else
        {
            Scalar* X2 = X1 + stride;  // X2 => X.col(2)
            const Scalar u2 = m_ref_u.coeff(2, u_ind);
            const Scalar u2_2 = Scalar(2) * u2;
            for (Index i = 0; i < nrow; i++)
            {
                const Scalar tmp = u0_2 * X0[i] + u1_2 * X1[i] + u2_2 * X2[i];
                X0[i] -= tmp * u0;
                X1[i] -= tmp * u1;
                X2[i] -= tmp * u2;
            }
        }
    }

public:
    DoubleShiftQR(Index size) :
        m_n(size),
        m_computed(false)
    {}

    DoubleShiftQR(ConstGenericMatrix& mat, const Scalar& s, const Scalar& t) :
        m_n(mat.rows()),
        m_mat_H(m_n, m_n),
        m_shift_s(s),
        m_shift_t(t),
        m_ref_u(3, m_n),
        m_ref_nr(m_n),
        m_computed(false)
    {
        compute(mat, s, t);
    }

    void compute(ConstGenericMatrix& mat, const Scalar& s, const Scalar& t)
    {
        using std::abs;

        m_n = mat.rows();
        if (m_n != mat.cols())
            throw std::invalid_argument("DoubleShiftQR: matrix must be square");

        m_mat_H.resize(m_n, m_n);
        m_shift_s = s;
        m_shift_t = t;
        m_ref_u.resize(3, m_n);
        m_ref_nr.resize(m_n);

        // Make a copy of mat
        m_mat_H.noalias() = mat;

        // Obtain the indices of zero elements in the subdiagonal,
        // so that H can be divided into several blocks
        const Scalar eps_abs = m_near_0 * (m_n / m_eps);
        constexpr Scalar eps_rel = m_eps;
        std::vector<int> zero_ind;
        zero_ind.reserve(m_n - 1);
        zero_ind.push_back(0);
        Scalar* Hii = m_mat_H.data();
        for (Index i = 0; i < m_n - 1; i++, Hii += (m_n + 1))
        {
            // Hii[0] => m_mat_H(i, i)
            // Hii[1] => m_mat_H(i + 1, i)
            // Hii[m_n + 1] => m_mat_H(i + 1, i + 1)
            const Scalar h = abs(Hii[1]);
            // Deflate small sub-diagonal elements
            const Scalar diag = abs(Hii[0]) + abs(Hii[m_n + 1]);
            if (h <= eps_abs || h <= eps_rel * diag)
            {
                Hii[1] = 0;
                zero_ind.push_back(i + 1);
            }
            // Make sure m_mat_H is upper Hessenberg
            // Zero the elements below m_mat_H(i + 1, i)
            std::fill(Hii + 2, Hii + m_n - i, Scalar(0));
        }
        zero_ind.push_back(m_n);

        const Index len = zero_ind.size() - 1;
        for (Index i = 0; i < len; i++)
        {
            const Index start = zero_ind[i];
            const Index end = zero_ind[i + 1] - 1;
            // Compute refelctors and update each block
            update_block(start, end);
        }

        // Deflation on the computed result
        Hii = m_mat_H.data();
        for (Index i = 0; i < m_n - 1; i++, Hii += (m_n + 1))
        {
            const Scalar h = abs(Hii[1]);
            const Scalar diag = abs(Hii[0]) + abs(Hii[m_n + 1]);
            if (h <= eps_abs || h <= eps_rel * diag)
                Hii[1] = 0;
        }

        m_computed = true;
    }

    void matrix_QtHQ(Matrix& dest) const
    {
        if (!m_computed)
            throw std::logic_error("DoubleShiftQR: need to call compute() first");

        dest.noalias() = m_mat_H;
    }

    // Q = P0 * P1 * ...
    // Q'y = P_{n-2} * ... * P1 * P0 * y
    void apply_QtY(Vector& y) const
    {
        if (!m_computed)
            throw std::logic_error("DoubleShiftQR: need to call compute() first");

        Scalar* y_ptr = y.data();
        const Index n1 = m_n - 1;
        for (Index i = 0; i < n1; i++, y_ptr++)
        {
            apply_PX(y_ptr, i);
        }
    }

    // Q = P0 * P1 * ...
    // YQ = Y * P0 * P1 * ...
    void apply_YQ(GenericMatrix Y) const
    {
        if (!m_computed)
            throw std::logic_error("DoubleShiftQR: need to call compute() first");

        const Index nrow = Y.rows();
        const Index n2 = m_n - 2;
        for (Index i = 0; i < n2; i++)
        {
            apply_XP(Y.block(0, i, nrow, 3), nrow, i);
        }
        apply_XP(Y.block(0, n2, nrow, 2), nrow, n2);
    }
};

}  // namespace Spectra

#endif  // SPECTRA_DOUBLE_SHIFT_QR_H
