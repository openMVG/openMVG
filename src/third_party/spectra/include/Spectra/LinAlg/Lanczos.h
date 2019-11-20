// Copyright (C) 2018-2019 Yixuan Qiu <yixuan.qiu@cos.name>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef LANCZOS_H
#define LANCZOS_H

#include <Eigen/Core>
#include <cmath>      // std::sqrt
#include <stdexcept>  // std::invalid_argument
#include <sstream>    // std::stringstream

#include "Arnoldi.h"

namespace Spectra {


// Lanczos factorization A * V = V * H + f * e'
// A: n x n
// V: n x k
// H: k x k
// f: n x 1
// e: [0, ..., 0, 1]
// V and H are allocated of dimension m, so the maximum value of k is m
template <typename Scalar, typename ArnoldiOpType>
class Lanczos: public Arnoldi<Scalar, ArnoldiOpType>
{
private:
    typedef Eigen::Index Index;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Matrix;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;
    typedef Eigen::Map<Matrix> MapMat;
    typedef Eigen::Map<Vector> MapVec;
    typedef Eigen::Map<const Matrix> MapConstMat;
    typedef Eigen::Map<const Vector> MapConstVec;

    using Arnoldi<Scalar, ArnoldiOpType>::m_op;
    using Arnoldi<Scalar, ArnoldiOpType>::m_n;
    using Arnoldi<Scalar, ArnoldiOpType>::m_m;
    using Arnoldi<Scalar, ArnoldiOpType>::m_k;
    using Arnoldi<Scalar, ArnoldiOpType>::m_fac_V;
    using Arnoldi<Scalar, ArnoldiOpType>::m_fac_H;
    using Arnoldi<Scalar, ArnoldiOpType>::m_fac_f;
    using Arnoldi<Scalar, ArnoldiOpType>::m_beta;
    using Arnoldi<Scalar, ArnoldiOpType>::m_near_0;
    using Arnoldi<Scalar, ArnoldiOpType>::m_eps;

public:
    Lanczos(const ArnoldiOpType& op, Index m) :
        Arnoldi<Scalar, ArnoldiOpType>(op, m)
    {}

    // Lanczos factorization starting from step-k
    void factorize_from(Index from_k, Index to_m, Index& op_counter)
    {
        using std::sqrt;

        if(to_m <= from_k) return;

        if(from_k > m_k)
        {
            std::stringstream msg;
            msg << "Lanczos: from_k (= " << from_k <<
                   ") is larger than the current subspace dimension (= " <<
                   m_k << ")";
            throw std::invalid_argument(msg.str());
        }

        const Scalar beta_thresh = m_eps * sqrt(Scalar(m_n));

        // Pre-allocate vectors
        Vector Vf(to_m);
        Vector w(m_n);

        // Keep the upperleft k x k submatrix of H and set other elements to 0
        m_fac_H.rightCols(m_m - from_k).setZero();
        m_fac_H.block(from_k, 0, m_m - from_k, from_k).setZero();

        for(Index i = from_k; i <= to_m - 1; i++)
        {
            bool restart = false;
            // If beta = 0, then the next V is not full rank
            // We need to generate a new residual vector that is orthogonal
            // to the current V, which we call a restart
            if(m_beta < m_near_0)
            {
                MapConstMat V(m_fac_V.data(), m_n, i); // The first i columns
                this->expand_basis(V, 2 * i, m_fac_f, m_beta);
                restart = true;
            }

            // v <- f / ||f||
            MapVec v(&m_fac_V(0, i), m_n); // The (i+1)-th column
            v.noalias() = m_fac_f / m_beta;

            // Note that H[i+1, i] equals to the unrestarted beta
            m_fac_H(i, i - 1) = restart ? Scalar(0) : m_beta;

            // w <- A * v
            m_op.perform_op(v.data(), w.data());
            op_counter++;

            // H[i+1, i+1] = <v, w> = v'Bw
            m_fac_H(i - 1, i) = m_fac_H(i, i - 1); // Due to symmetry
            m_fac_H(i, i) = m_op.inner_product(v, w);

            // f <- w - V * V'Bw = w - H[i+1, i] * V{i} - H[i+1, i+1] * V{i+1}
            // If restarting, we know that H[i+1, i] = 0
            if(restart)
                m_fac_f.noalias() = w - m_fac_H(i, i) * v;
            else
                m_fac_f.noalias() = w - m_fac_H(i, i - 1) * m_fac_V.col(i - 1) - m_fac_H(i, i) * v;

            m_beta = m_op.norm(m_fac_f);

            // f/||f|| is going to be the next column of V, so we need to test
            // whether V'B(f/||f||) ~= 0
            const Index i1 = i + 1;
            MapMat Vs(m_fac_V.data(), m_n, i1); // The first (i+1) columns
            m_op.trans_product(Vs, m_fac_f, Vf.head(i1));
            Scalar ortho_err = Vf.head(i1).cwiseAbs().maxCoeff();
            // If not, iteratively correct the residual
            int count = 0;
            while(count < 5 && ortho_err > m_eps * m_beta)
            {
                // There is an edge case: when beta=||f|| is close to zero, f mostly consists
                // of noises of rounding errors, so the test [ortho_err < eps * beta] is very
                // likely to fail. In particular, if beta=0, then the test is ensured to fail.
                // Hence when this happens, we force f to be zero, and then restart in the
                // next iteration.
                if(m_beta < beta_thresh)
                {
                    m_fac_f.setZero();
                    m_beta = Scalar(0);
                    break;
                }

                // f <- f - V * Vf
                m_fac_f.noalias() -= Vs * Vf.head(i1);
                // h <- h + Vf
                m_fac_H(i - 1, i) += Vf[i - 1];
                m_fac_H(i, i - 1) = m_fac_H(i - 1, i);
                m_fac_H(i, i) += Vf[i];
                // beta <- ||f||
                m_beta = m_op.norm(m_fac_f);

                m_op.trans_product(Vs, m_fac_f, Vf.head(i1));
                ortho_err = Vf.head(i1).cwiseAbs().maxCoeff();
                count++;
            }
        }

        // Indicate that this is a step-m factorization
        m_k = to_m;
    }

    // Apply H -> Q'HQ, where Q is from a tridiagonal QR decomposition
    void compress_H(const TridiagQR<Scalar>& decomp)
    {
        decomp.matrix_QtHQ(m_fac_H);
        m_k--;
    }
};


} // namespace Spectra

#endif // LANCZOS_H
