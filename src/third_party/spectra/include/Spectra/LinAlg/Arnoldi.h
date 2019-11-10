// Copyright (C) 2018-2019 Yixuan Qiu <yixuan.qiu@cos.name>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef ARNOLDI_H
#define ARNOLDI_H

#include <Eigen/Core>
#include <cmath>      // std::sqrt
#include <stdexcept>  // std::invalid_argument
#include <sstream>    // std::stringstream

#include "../MatOp/internal/ArnoldiOp.h"
#include "../Util/TypeTraits.h"
#include "../Util/SimpleRandom.h"
#include "UpperHessenbergQR.h"
#include "DoubleShiftQR.h"

namespace Spectra {


// Arnoldi factorization A * V = V * H + f * e'
// A: n x n
// V: n x k
// H: k x k
// f: n x 1
// e: [0, ..., 0, 1]
// V and H are allocated of dimension m, so the maximum value of k is m
template <typename Scalar, typename ArnoldiOpType>
class Arnoldi
{
private:
    typedef Eigen::Index Index;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Matrix;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;
    typedef Eigen::Map<Matrix> MapMat;
    typedef Eigen::Map<Vector> MapVec;
    typedef Eigen::Map<const Matrix> MapConstMat;
    typedef Eigen::Map<const Vector> MapConstVec;

protected:
    ArnoldiOpType m_op;       // Operators for the Arnoldi factorization

    const Index m_n;          // dimension of A
    const Index m_m;          // maximum dimension of subspace V
    Index       m_k;          // current dimension of subspace V

    Matrix m_fac_V;           // V matrix in the Arnoldi factorization
    Matrix m_fac_H;           // H matrix in the Arnoldi factorization
    Vector m_fac_f;           // residual in the Arnoldi factorization
    Scalar m_beta;            // ||f||, B-norm of f

    const Scalar m_near_0;    // a very small value, but 1.0 / m_near_0 does not overflow
                              // ~= 1e-307 for the "double" type
    const Scalar m_eps;       // the machine precision, ~= 1e-16 for the "double" type

    // Given orthonormal basis functions V, find a nonzero vector f such that V'Bf = 0
    // Assume that f has been properly allocated
    void expand_basis(MapConstMat& V, const Index seed, Vector& f, Scalar& fnorm)
    {
        using std::sqrt;

        const Scalar thresh = m_eps * sqrt(Scalar(m_n));
        Vector Vf(V.cols());
        for(Index iter = 0; iter < 5; iter++)
        {
            // Randomly generate a new vector and orthogonalize it against V
            SimpleRandom<Scalar> rng(seed + 123 * iter);
            f.noalias() = rng.random_vec(m_n);
            // f <- f - V * V'Bf, so that f is orthogonal to V in B-norm
            m_op.trans_product(V, f, Vf);
            f.noalias() -= V * Vf;
            // fnorm <- ||f||
            fnorm = m_op.norm(f);

            // If fnorm is too close to zero, we try a new random vector,
            // otherwise return the result
            if(fnorm >= thresh)
                return;
        }
    }

public:
    Arnoldi(const ArnoldiOpType& op, Index m) :
        m_op(op), m_n(op.rows()), m_m(m), m_k(0),
        m_near_0(TypeTraits<Scalar>::min() * Scalar(10)),
        m_eps(Eigen::NumTraits<Scalar>::epsilon())
    {}

    virtual ~Arnoldi() {}

    // Const-reference to internal structures
    const Matrix& matrix_V() const { return m_fac_V; }
    const Matrix& matrix_H() const { return m_fac_H; }
    const Vector& vector_f() const { return m_fac_f; }
    Scalar f_norm() const { return m_beta; }
    Index subspace_dim() const { return m_k; }

    // Initialize with an operator and an initial vector
    void init(MapConstVec& v0, Index& op_counter)
    {
        m_fac_V.resize(m_n, m_m);
        m_fac_H.resize(m_m, m_m);
        m_fac_f.resize(m_n);
        m_fac_H.setZero();

        // Verify the initial vector
        const Scalar v0norm = m_op.norm(v0);
        if(v0norm < m_near_0)
            throw std::invalid_argument("initial residual vector cannot be zero");

        // Points to the first column of V
        MapVec v(m_fac_V.data(), m_n);

        // Normalize
        v.noalias() = v0 / v0norm;

        // Compute H and f
        Vector w(m_n);
        m_op.perform_op(v.data(), w.data());
        op_counter++;

        m_fac_H(0, 0) = m_op.inner_product(v, w);
        m_fac_f.noalias() = w - v * m_fac_H(0, 0);

        // In some cases f is zero in exact arithmetics, but due to rounding errors
        // it may contain tiny fluctuations. When this happens, we force f to be zero
        if(m_fac_f.cwiseAbs().maxCoeff() < m_eps)
        {
            m_fac_f.setZero();
            m_beta = Scalar(0);
        } else {
            m_beta = m_op.norm(m_fac_f);
        }

        // Indicate that this is a step-1 factorization
        m_k = 1;
    }

    // Arnoldi factorization starting from step-k
    virtual void factorize_from(Index from_k, Index to_m, Index& op_counter)
    {
        using std::sqrt;

        if(to_m <= from_k) return;

        if(from_k > m_k)
        {
            std::stringstream msg;
            msg << "Arnoldi: from_k (= " << from_k <<
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
                expand_basis(V, 2 * i, m_fac_f, m_beta);
                restart = true;
            }

            // v <- f / ||f||
            m_fac_V.col(i).noalias() = m_fac_f / m_beta; // The (i+1)-th column

            // Note that H[i+1, i] equals to the unrestarted beta
            m_fac_H(i, i - 1) = restart ? Scalar(0) : m_beta;

            // w <- A * v, v = m_fac_V.col(i)
            m_op.perform_op(&m_fac_V(0, i), w.data());
            op_counter++;

            const Index i1 = i + 1;
            // First i+1 columns of V
            MapConstMat Vs(m_fac_V.data(), m_n, i1);
            // h = m_fac_H(0:i, i)
            MapVec h(&m_fac_H(0, i), i1);
            // h <- V'Bw
            m_op.trans_product(Vs, w, h);

            // f <- w - V * h
            m_fac_f.noalias() = w - Vs * h;
            m_beta = m_op.norm(m_fac_f);

            if(m_beta > Scalar(0.717) * m_op.norm(h))
                continue;

            // f/||f|| is going to be the next column of V, so we need to test
            // whether V'B(f/||f||) ~= 0
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
                h.noalias() += Vf.head(i1);
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

    // Apply H -> Q'HQ, where Q is from a double shift QR decomposition
    void compress_H(const DoubleShiftQR<Scalar>& decomp)
    {
        decomp.matrix_QtHQ(m_fac_H);
        m_k -= 2;
    }

    // Apply H -> Q'HQ, where Q is from an upper Hessenberg QR decomposition
    void compress_H(const UpperHessenbergQR<Scalar>& decomp)
    {
        decomp.matrix_QtHQ(m_fac_H);
        m_k--;
    }

    // Apply V -> VQ and compute the new f.
    // Should be called after compress_H(), since m_k is updated there.
    // Only need to update the first k+1 columns of V
    // The first (m - k + i) elements of the i-th column of Q are non-zero,
    // and the rest are zero
    void compress_V(const Matrix& Q)
    {
        Matrix Vs(m_n, m_k + 1);
        for(Index i = 0; i < m_k; i++)
        {
            const Index nnz = m_m - m_k + i + 1;
            MapConstVec q(&Q(0, i), nnz);
            Vs.col(i).noalias() = m_fac_V.leftCols(nnz) * q;
        }
        Vs.col(m_k).noalias() = m_fac_V * Q.col(m_k);
        m_fac_V.leftCols(m_k + 1).noalias() = Vs;

        Vector fk = m_fac_f * Q(m_m - 1, m_k - 1) + m_fac_V.col(m_k) * m_fac_H(m_k, m_k - 1);
        m_fac_f.swap(fk);
        m_beta = m_op.norm(m_fac_f);
    }
};


} // namespace Spectra

#endif // ARNOLDI_H
