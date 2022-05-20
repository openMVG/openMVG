// The code was adapted from Eigen/src/Eigenvaleus/RealSchur.h
//
// Copyright (C) 2008 Gael Guennebaud <gael.guennebaud@inria.fr>
// Copyright (C) 2010,2012 Jitse Niesen <jitse@maths.leeds.ac.uk>
// Copyright (C) 2021-2022 Yixuan Qiu <yixuan.qiu@cos.name>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SPECTRA_UPPER_HESSENBERG_SCHUR_H
#define SPECTRA_UPPER_HESSENBERG_SCHUR_H

#include <Eigen/Core>
#include <Eigen/Jacobi>
#include <Eigen/Householder>
#include <stdexcept>

#include "../Util/TypeTraits.h"

namespace Spectra {

template <typename Scalar = double>
class UpperHessenbergSchur
{
private:
    using Index = Eigen::Index;
    using Matrix = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
    using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
    using Vector2s = Eigen::Matrix<Scalar, 2, 1>;
    using Vector3s = Eigen::Matrix<Scalar, 3, 1>;
    using GenericMatrix = Eigen::Ref<Matrix>;
    using ConstGenericMatrix = const Eigen::Ref<const Matrix>;

    Index m_n;   // Size of the matrix
    Matrix m_T;  // T matrix, A = UTU'
    Matrix m_U;  // U matrix, A = UTU'
    bool m_computed;

    // L1 norm of an upper Hessenberg matrix
    static Scalar upper_hessenberg_l1_norm(ConstGenericMatrix& x)
    {
        const Index n = x.cols();
        Scalar norm(0);
        for (Index j = 0; j < n; j++)
            norm += x.col(j).segment(0, (std::min)(n, j + 2)).cwiseAbs().sum();
        return norm;
    }

    // Look for single small sub-diagonal element and returns its index
    Index find_small_subdiag(Index iu, const Scalar& near_0) const
    {
        using std::abs;

        const Scalar eps = Eigen::NumTraits<Scalar>::epsilon();
        Index res = iu;
        while (res > 0)
        {
            Scalar s = abs(m_T.coeff(res - 1, res - 1)) + abs(m_T.coeff(res, res));
            s = Eigen::numext::maxi<Scalar>(s * eps, near_0);
            if (abs(m_T.coeff(res, res - 1)) <= s)
                break;
            res--;
        }

        return res;
    }

    // Update T given that rows iu-1 and iu decouple from the rest
    void split_off_two_rows(Index iu, const Scalar& ex_shift)
    {
        using std::sqrt;
        using std::abs;

        // The eigenvalues of the 2x2 matrix [a b; c d] are
        // trace +/- sqrt(discr/4) where discr = tr^2 - 4*det, tr = a + d, det = ad - bc
        Scalar p = Scalar(0.5) * (m_T.coeff(iu - 1, iu - 1) - m_T.coeff(iu, iu));
        Scalar q = p * p + m_T.coeff(iu, iu - 1) * m_T.coeff(iu - 1, iu);  // q = tr^2 / 4 - det = discr/4
        m_T.coeffRef(iu, iu) += ex_shift;
        m_T.coeffRef(iu - 1, iu - 1) += ex_shift;

        if (q >= Scalar(0))  // Two real eigenvalues
        {
            Scalar z = sqrt(abs(q));
            Eigen::JacobiRotation<Scalar> rot;
            rot.makeGivens((p >= Scalar(0)) ? (p + z) : (p - z), m_T.coeff(iu, iu - 1));
            m_T.rightCols(m_n - iu + 1).applyOnTheLeft(iu - 1, iu, rot.adjoint());
            m_T.topRows(iu + 1).applyOnTheRight(iu - 1, iu, rot);
            m_T.coeffRef(iu, iu - 1) = Scalar(0);
            m_U.applyOnTheRight(iu - 1, iu, rot);
        }
        if (iu > 1)
            m_T.coeffRef(iu - 1, iu - 2) = Scalar(0);
    }

    // Form shift in shift_info, and update ex_shift if an exceptional shift is performed
    void compute_shift(Index iu, Index iter, Scalar& ex_shift, Vector3s& shift_info)
    {
        using std::sqrt;
        using std::abs;

        shift_info.coeffRef(0) = m_T.coeff(iu, iu);
        shift_info.coeffRef(1) = m_T.coeff(iu - 1, iu - 1);
        shift_info.coeffRef(2) = m_T.coeff(iu, iu - 1) * m_T.coeff(iu - 1, iu);

        // Wilkinson's original ad hoc shift
        if (iter == 10)
        {
            ex_shift += shift_info.coeff(0);
            for (Index i = 0; i <= iu; ++i)
                m_T.coeffRef(i, i) -= shift_info.coeff(0);
            Scalar s = abs(m_T.coeff(iu, iu - 1)) + abs(m_T.coeff(iu - 1, iu - 2));
            shift_info.coeffRef(0) = Scalar(0.75) * s;
            shift_info.coeffRef(1) = Scalar(0.75) * s;
            shift_info.coeffRef(2) = Scalar(-0.4375) * s * s;
        }

        // MATLAB's new ad hoc shift
        if (iter == 30)
        {
            Scalar s = (shift_info.coeff(1) - shift_info.coeff(0)) / Scalar(2);
            s = s * s + shift_info.coeff(2);
            if (s > Scalar(0))
            {
                s = sqrt(s);
                if (shift_info.coeff(1) < shift_info.coeff(0))
                    s = -s;
                s = s + (shift_info.coeff(1) - shift_info.coeff(0)) / Scalar(2);
                s = shift_info.coeff(0) - shift_info.coeff(2) / s;
                ex_shift += s;
                for (Index i = 0; i <= iu; ++i)
                    m_T.coeffRef(i, i) -= s;
                shift_info.setConstant(Scalar(0.964));
            }
        }
    }

    // Compute index im at which Francis QR step starts and the first Householder vector
    void init_francis_qr_step(Index il, Index iu, const Vector3s& shift_info, Index& im, Vector3s& first_householder_vec) const
    {
        using std::abs;

        const Scalar eps = Eigen::NumTraits<Scalar>::epsilon();
        Vector3s& v = first_householder_vec;  // alias to save typing
        for (im = iu - 2; im >= il; --im)
        {
            const Scalar Tmm = m_T.coeff(im, im);
            const Scalar r = shift_info.coeff(0) - Tmm;
            const Scalar s = shift_info.coeff(1) - Tmm;
            v.coeffRef(0) = (r * s - shift_info.coeff(2)) / m_T.coeff(im + 1, im) + m_T.coeff(im, im + 1);
            v.coeffRef(1) = m_T.coeff(im + 1, im + 1) - Tmm - r - s;
            v.coeffRef(2) = m_T.coeff(im + 2, im + 1);
            if (im == il)
                break;
            const Scalar lhs = m_T.coeff(im, im - 1) * (abs(v.coeff(1)) + abs(v.coeff(2)));
            const Scalar rhs = v.coeff(0) * (abs(m_T.coeff(im - 1, im - 1)) + abs(Tmm) + abs(m_T.coeff(im + 1, im + 1)));
            if (abs(lhs) < eps * rhs)
                break;
        }
    }

    // P = I - tau * v * v' = P'
    // PX = X - tau * v * (v'X), X [3 x c]
    static void apply_householder_left(const Vector2s& ess, const Scalar& tau, Scalar* x, Index ncol, Index stride)
    {
        const Scalar v1 = ess.coeff(0), v2 = ess.coeff(1);
        const Scalar* const x_end = x + ncol * stride;
        for (; x < x_end; x += stride)
        {
            const Scalar tvx = tau * (x[0] + v1 * x[1] + v2 * x[2]);
            x[0] -= tvx;
            x[1] -= tvx * v1;
            x[2] -= tvx * v2;
        }
    }

    // P = I - tau * v * v' = P'
    // XP = X - tau * (X * v) * v', X [r x 3]
    static void apply_householder_right(const Vector2s& ess, const Scalar& tau, Scalar* x, Index nrow, Index stride)
    {
        const Scalar v1 = ess.coeff(0), v2 = ess.coeff(1);
        Scalar* x0 = x;
        Scalar* x1 = x + stride;
        Scalar* x2 = x1 + stride;
        for (Index i = 0; i < nrow; i++)
        {
            const Scalar txv = tau * (x0[i] + v1 * x1[i] + v2 * x2[i]);
            x0[i] -= txv;
            x1[i] -= txv * v1;
            x2[i] -= txv * v2;
        }
    }

    // SIMD version of apply_householder_right()
    // Inspired by apply_rotation_in_the_plane_selector() in Eigen/src/Jacobi/Jacobi.h
    static void apply_householder_right_simd(const Vector2s& ess, const Scalar& tau, Scalar* x, Index nrow, Index stride)
    {
        // Packet type
        using Eigen::internal::ploadu;
        using Eigen::internal::pstoreu;
        using Eigen::internal::pset1;
        using Eigen::internal::padd;
        using Eigen::internal::psub;
        using Eigen::internal::pmul;
        using Packet = typename Eigen::internal::packet_traits<Scalar>::type;
        constexpr unsigned char PacketSize = Eigen::internal::packet_traits<Scalar>::size;
        constexpr unsigned char Peeling = 2;
        constexpr unsigned char Increment = Peeling * PacketSize;

        // Column heads
        Scalar* x0 = x;
        Scalar* x1 = x + stride;
        Scalar* x2 = x1 + stride;
        // Pointers for the current row
        Scalar* px0 = x0;
        Scalar* px1 = x1;
        Scalar* px2 = x2;

        // Householder reflectors
        const Scalar v1 = ess.coeff(0), v2 = ess.coeff(1);
        // Vectorized versions
        const Packet vtau = pset1<Packet>(tau);
        const Packet vv1 = pset1<Packet>(v1);
        const Packet vv2 = pset1<Packet>(v2);

        // n % (2^k) == n & (2^k-1), see https://stackoverflow.com/q/3072665
        // const Index peeling_end = nrow - nrow % Increment;
        const Index aligned_end = nrow - (nrow & (PacketSize - 1));
        const Index peeling_end = nrow - (nrow & (Increment - 1));
        for (Index i = 0; i < peeling_end; i += Increment)
        {
            Packet vx01 = ploadu<Packet>(px0);
            Packet vx02 = ploadu<Packet>(px0 + PacketSize);
            Packet vx11 = ploadu<Packet>(px1);
            Packet vx12 = ploadu<Packet>(px1 + PacketSize);
            Packet vx21 = ploadu<Packet>(px2);
            Packet vx22 = ploadu<Packet>(px2 + PacketSize);

            // Packet txv1 = vtau * (vx01 + vv1 * vx11 + vv2 * vx21);
            Packet txv1 = pmul(vtau, padd(padd(vx01, pmul(vv1, vx11)), pmul(vv2, vx21)));
            Packet txv2 = pmul(vtau, padd(padd(vx02, pmul(vv1, vx12)), pmul(vv2, vx22)));

            pstoreu(px0, psub(vx01, txv1));
            pstoreu(px0 + PacketSize, psub(vx02, txv2));
            pstoreu(px1, psub(vx11, pmul(txv1, vv1)));
            pstoreu(px1 + PacketSize, psub(vx12, pmul(txv2, vv1)));
            pstoreu(px2, psub(vx21, pmul(txv1, vv2)));
            pstoreu(px2 + PacketSize, psub(vx22, pmul(txv2, vv2)));

            px0 += Increment;
            px1 += Increment;
            px2 += Increment;
        }
        if (aligned_end != peeling_end)
        {
            px0 = x0 + peeling_end;
            px1 = x1 + peeling_end;
            px2 = x2 + peeling_end;

            Packet x0_p = ploadu<Packet>(px0);
            Packet x1_p = ploadu<Packet>(px1);
            Packet x2_p = ploadu<Packet>(px2);
            Packet txv = pmul(vtau, padd(padd(x0_p, pmul(vv1, x1_p)), pmul(vv2, x2_p)));

            pstoreu(px0, psub(x0_p, txv));
            pstoreu(px1, psub(x1_p, pmul(txv, vv1)));
            pstoreu(px2, psub(x2_p, pmul(txv, vv2)));
        }

        // Remaining rows
        for (Index i = aligned_end; i < nrow; i++)
        {
            const Scalar txv = tau * (x0[i] + v1 * x1[i] + v2 * x2[i]);
            x0[i] -= txv;
            x1[i] -= txv * v1;
            x2[i] -= txv * v2;
        }
    }

    // Perform a Francis QR step involving rows il:iu and columns im:iu
    void perform_francis_qr_step(Index il, Index im, Index iu, const Vector3s& first_householder_vec, const Scalar& near_0)
    {
        using std::abs;

        for (Index k = im; k <= iu - 2; ++k)
        {
            const bool first_iter = (k == im);
            Vector3s v;
            if (first_iter)
                v = first_householder_vec;
            else
                v = m_T.template block<3, 1>(k, k - 1);

            Scalar tau, beta;
            Vector2s ess;
            v.makeHouseholder(ess, tau, beta);

            if (abs(beta) > near_0)  // if v is not zero
            {
                if (first_iter && k > il)
                    m_T.coeffRef(k, k - 1) = -m_T.coeff(k, k - 1);
                else if (!first_iter)
                    m_T.coeffRef(k, k - 1) = beta;

                // These Householder transformations form the O(n^3) part of the algorithm
                // m_T.block(k, k, 3, m_n - k).applyHouseholderOnTheLeft(ess, tau, workspace);
                // m_T.block(0, k, (std::min)(iu, k + 3) + 1, 3).applyHouseholderOnTheRight(ess, tau, workspace);
                // m_U.block(0, k, m_n, 3).applyHouseholderOnTheRight(ess, tau, workspace);
                apply_householder_left(ess, tau, &m_T.coeffRef(k, k), m_n - k, m_n);
                apply_householder_right_simd(ess, tau, &m_T.coeffRef(0, k), (std::min)(iu, k + 3) + 1, m_n);
                apply_householder_right_simd(ess, tau, &m_U.coeffRef(0, k), m_n, m_n);
            }
        }

        // The last 2-row block
        Eigen::JacobiRotation<Scalar> rot;
        Scalar beta;
        rot.makeGivens(m_T.coeff(iu - 1, iu - 2), m_T.coeff(iu, iu - 2), &beta);

        if (abs(beta) > near_0)  // if v is not zero
        {
            m_T.coeffRef(iu - 1, iu - 2) = beta;
            m_T.rightCols(m_n - iu + 1).applyOnTheLeft(iu - 1, iu, rot.adjoint());
            m_T.topRows(iu + 1).applyOnTheRight(iu - 1, iu, rot);
            m_U.applyOnTheRight(iu - 1, iu, rot);
        }

        // clean up pollution due to round-off errors
        for (Index i = im + 2; i <= iu; ++i)
        {
            m_T.coeffRef(i, i - 2) = Scalar(0);
            if (i > im + 2)
                m_T.coeffRef(i, i - 3) = Scalar(0);
        }
    }

public:
    UpperHessenbergSchur() :
        m_n(0), m_computed(false)
    {}

    UpperHessenbergSchur(ConstGenericMatrix& mat) :
        m_n(mat.rows()), m_computed(false)
    {
        compute(mat);
    }

    void compute(ConstGenericMatrix& mat)
    {
        using std::abs;
        using std::sqrt;

        if (mat.rows() != mat.cols())
            throw std::invalid_argument("UpperHessenbergSchur: matrix must be square");

        m_n = mat.rows();
        m_T.resize(m_n, m_n);
        m_U.resize(m_n, m_n);
        constexpr Index max_iter_per_row = 40;
        const Index max_iter = m_n * max_iter_per_row;

        m_T.noalias() = mat;
        m_U.setIdentity();

        // The matrix m_T is divided in three parts.
        // Rows 0,...,il-1 are decoupled from the rest because m_T(il,il-1) is zero.
        // Rows il,...,iu is the part we are working on (the active window).
        // Rows iu+1,...,end are already brought in triangular form.
        Index iu = m_n - 1;
        Index iter = 0;        // iteration count for current eigenvalue
        Index total_iter = 0;  // iteration count for whole matrix
        Scalar ex_shift(0);    // sum of exceptional shifts
        const Scalar norm = upper_hessenberg_l1_norm(m_T);
        // sub-diagonal entries smaller than near_0 will be treated as zero.
        // We use eps^2 to enable more precision in small eigenvalues.
        const Scalar eps = Eigen::NumTraits<Scalar>::epsilon();
        const Scalar near_0 = Eigen::numext::maxi<Scalar>(norm * eps * eps, TypeTraits<Scalar>::min());

        if (norm != Scalar(0))
        {
            while (iu >= 0)
            {
                Index il = find_small_subdiag(iu, near_0);

                // Check for convergence
                if (il == iu)  // One root found
                {
                    m_T.coeffRef(iu, iu) += ex_shift;
                    if (iu > 0)
                        m_T.coeffRef(iu, iu - 1) = Scalar(0);
                    iu--;
                    iter = 0;
                }
                else if (il == iu - 1)  // Two roots found
                {
                    split_off_two_rows(iu, ex_shift);
                    iu -= 2;
                    iter = 0;
                }
                else  // No convergence yet
                {
                    Vector3s first_householder_vec = Vector3s::Zero(), shift_info;
                    compute_shift(iu, iter, ex_shift, shift_info);
                    iter++;
                    total_iter++;
                    if (total_iter > max_iter)
                        break;
                    Index im;
                    init_francis_qr_step(il, iu, shift_info, im, first_householder_vec);
                    perform_francis_qr_step(il, im, iu, first_householder_vec, near_0);
                }
            }
        }

        if (total_iter > max_iter)
            throw std::runtime_error("UpperHessenbergSchur: Schur decomposition failed");

        m_computed = true;
    }

    const Matrix& matrix_T() const
    {
        if (!m_computed)
            throw std::logic_error("UpperHessenbergSchur: need to call compute() first");

        return m_T;
    }

    const Matrix& matrix_U() const
    {
        if (!m_computed)
            throw std::logic_error("UpperHessenbergSchur: need to call compute() first");

        return m_U;
    }

    void swap_T(Matrix& other)
    {
        m_T.swap(other);
    }

    void swap_U(Matrix& other)
    {
        m_U.swap(other);
    }
};

}  // namespace Spectra

#endif  // SPECTRA_UPPER_HESSENBERG_SCHUR_H
