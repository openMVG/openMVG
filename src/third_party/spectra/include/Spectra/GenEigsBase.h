// Copyright (C) 2018-2019 Yixuan Qiu <yixuan.qiu@cos.name>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef GEN_EIGS_BASE_H
#define GEN_EIGS_BASE_H

#include <Eigen/Core>
#include <vector>     // std::vector
#include <cmath>      // std::abs, std::pow, std::sqrt
#include <algorithm>  // std::min, std::copy
#include <complex>    // std::complex, std::conj, std::norm, std::abs
#include <stdexcept>  // std::invalid_argument

#include "Util/TypeTraits.h"
#include "Util/SelectionRule.h"
#include "Util/CompInfo.h"
#include "Util/SimpleRandom.h"
#include "MatOp/internal/ArnoldiOp.h"
#include "LinAlg/UpperHessenbergQR.h"
#include "LinAlg/DoubleShiftQR.h"
#include "LinAlg/UpperHessenbergEigen.h"
#include "LinAlg/Arnoldi.h"

namespace Spectra {


///
/// \ingroup EigenSolver
///
/// This is the base class for general eigen solvers, mainly for internal use.
/// It is kept here to provide the documentation for member functions of concrete eigen solvers
/// such as GenEigsSolver and GenEigsRealShiftSolver.
///
template < typename Scalar,
           int      SelectionRule,
           typename OpType,
           typename BOpType >
class GenEigsBase
{
private:
    typedef Eigen::Index Index;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Matrix;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;
    typedef Eigen::Array<Scalar, Eigen::Dynamic, 1> Array;
    typedef Eigen::Array<bool, Eigen::Dynamic, 1> BoolArray;
    typedef Eigen::Map<Matrix> MapMat;
    typedef Eigen::Map<Vector> MapVec;
    typedef Eigen::Map<const Vector> MapConstVec;

    typedef std::complex<Scalar> Complex;
    typedef Eigen::Matrix<Complex, Eigen::Dynamic, Eigen::Dynamic> ComplexMatrix;
    typedef Eigen::Matrix<Complex, Eigen::Dynamic, 1> ComplexVector;

    typedef ArnoldiOp<Scalar, OpType, BOpType> ArnoldiOpType;
    typedef Arnoldi<Scalar, ArnoldiOpType> ArnoldiFac;

protected:
    OpType*       m_op;        // object to conduct matrix operation,
                               // e.g. matrix-vector product
    const Index   m_n;         // dimension of matrix A
    const Index   m_nev;       // number of eigenvalues requested
    const Index   m_ncv;       // dimension of Krylov subspace in the Arnoldi method
    Index         m_nmatop;    // number of matrix operations called
    Index         m_niter;     // number of restarting iterations

    ArnoldiFac    m_fac;       // Arnoldi factorization

    ComplexVector m_ritz_val;  // Ritz values
    ComplexMatrix m_ritz_vec;  // Ritz vectors
    ComplexVector m_ritz_est;  // last row of m_ritz_vec

private:
    BoolArray     m_ritz_conv; // indicator of the convergence of Ritz values
    int           m_info;      // status of the computation

    const Scalar  m_near_0;    // a very small value, but 1.0 / m_near_0 does not overflow
                               // ~= 1e-307 for the "double" type
    const Scalar  m_eps;       // the machine precision, ~= 1e-16 for the "double" type
    const Scalar  m_eps23;     // m_eps^(2/3), used to test the convergence

    // Real Ritz values calculated from UpperHessenbergEigen have exact zero imaginary part
    // Complex Ritz values have exact conjugate pairs
    // So we use exact tests here
    static bool is_complex(const Complex& v) { return v.imag() != Scalar(0); }
    static bool is_conj(const Complex& v1, const Complex& v2) { return v1 == Eigen::numext::conj(v2); }

    // Implicitly restarted Arnoldi factorization
    void restart(Index k)
    {
        using std::norm;

        if(k >= m_ncv)
            return;

        DoubleShiftQR<Scalar> decomp_ds(m_ncv);
        UpperHessenbergQR<Scalar> decomp_hb(m_ncv);
        Matrix Q = Matrix::Identity(m_ncv, m_ncv);

        for(Index i = k; i < m_ncv; i++)
        {
            if(is_complex(m_ritz_val[i]) && is_conj(m_ritz_val[i], m_ritz_val[i + 1]))
            {
                // H - mu * I = Q1 * R1
                // H <- R1 * Q1 + mu * I = Q1' * H * Q1
                // H - conj(mu) * I = Q2 * R2
                // H <- R2 * Q2 + conj(mu) * I = Q2' * H * Q2
                //
                // (H - mu * I) * (H - conj(mu) * I) = Q1 * Q2 * R2 * R1 = Q * R
                const Scalar s = Scalar(2) * m_ritz_val[i].real();
                const Scalar t = norm(m_ritz_val[i]);

                decomp_ds.compute(m_fac.matrix_H(), s, t);

                // Q -> Q * Qi
                decomp_ds.apply_YQ(Q);
                // H -> Q'HQ
                // Matrix Q = Matrix::Identity(m_ncv, m_ncv);
                // decomp_ds.apply_YQ(Q);
                // m_fac_H = Q.transpose() * m_fac_H * Q;
                m_fac.compress_H(decomp_ds);

                i++;
            } else {
                // QR decomposition of H - mu * I, mu is real
                decomp_hb.compute(m_fac.matrix_H(), m_ritz_val[i].real());

                // Q -> Q * Qi
                decomp_hb.apply_YQ(Q);
                // H -> Q'HQ = RQ + mu * I
                m_fac.compress_H(decomp_hb);
            }
        }

        m_fac.compress_V(Q);
        m_fac.factorize_from(k, m_ncv, m_nmatop);

        retrieve_ritzpair();
    }

    // Calculates the number of converged Ritz values
    Index num_converged(Scalar tol)
    {
        // thresh = tol * max(m_eps23, abs(theta)), theta for Ritz value
        Array thresh = tol * m_ritz_val.head(m_nev).array().abs().max(m_eps23);
        Array resid = m_ritz_est.head(m_nev).array().abs() * m_fac.f_norm();
        // Converged "wanted" Ritz values
        m_ritz_conv = (resid < thresh);

        return m_ritz_conv.cast<Index>().sum();
    }

    // Returns the adjusted nev for restarting
    Index nev_adjusted(Index nconv)
    {
        using std::abs;

        Index nev_new = m_nev;
        for(Index i = m_nev; i < m_ncv; i++)
            if(abs(m_ritz_est[i]) < m_near_0)  nev_new++;

        // Adjust nev_new, according to dnaup2.f line 660~674 in ARPACK
        nev_new += std::min(nconv, (m_ncv - nev_new) / 2);
        if(nev_new == 1 && m_ncv >= 6)
            nev_new = m_ncv / 2;
        else if(nev_new == 1 && m_ncv > 3)
            nev_new = 2;

        if(nev_new > m_ncv - 2)
            nev_new = m_ncv - 2;

        // Increase nev by one if ritz_val[nev - 1] and
        // ritz_val[nev] are conjugate pairs
        if(is_complex(m_ritz_val[nev_new - 1]) &&
           is_conj(m_ritz_val[nev_new - 1], m_ritz_val[nev_new]))
        {
            nev_new++;
        }

        return nev_new;
    }

    // Retrieves and sorts Ritz values and Ritz vectors
    void retrieve_ritzpair()
    {
        UpperHessenbergEigen<Scalar> decomp(m_fac.matrix_H());
        const ComplexVector& evals = decomp.eigenvalues();
        ComplexMatrix evecs = decomp.eigenvectors();

        SortEigenvalue<Complex, SelectionRule> sorting(evals.data(), evals.size());
        std::vector<int> ind = sorting.index();

        // Copy the Ritz values and vectors to m_ritz_val and m_ritz_vec, respectively
        for(Index i = 0; i < m_ncv; i++)
        {
            m_ritz_val[i] = evals[ind[i]];
            m_ritz_est[i] = evecs(m_ncv - 1, ind[i]);
        }
        for(Index i = 0; i < m_nev; i++)
        {
            m_ritz_vec.col(i).noalias() = evecs.col(ind[i]);
        }
    }

protected:
    // Sorts the first nev Ritz pairs in the specified order
    // This is used to return the final results
    virtual void sort_ritzpair(int sort_rule)
    {
        // First make sure that we have a valid index vector
        SortEigenvalue<Complex, LARGEST_MAGN> sorting(m_ritz_val.data(), m_nev);
        std::vector<int> ind = sorting.index();

        switch(sort_rule)
        {
            case LARGEST_MAGN:
                break;
            case LARGEST_REAL:
            {
                SortEigenvalue<Complex, LARGEST_REAL> sorting(m_ritz_val.data(), m_nev);
                ind = sorting.index();
            }
                break;
            case LARGEST_IMAG:
            {
                SortEigenvalue<Complex, LARGEST_IMAG> sorting(m_ritz_val.data(), m_nev);
                ind = sorting.index();
            }
                break;
            case SMALLEST_MAGN:
            {
                SortEigenvalue<Complex, SMALLEST_MAGN> sorting(m_ritz_val.data(), m_nev);
                ind = sorting.index();
            }
                break;
            case SMALLEST_REAL:
            {
                SortEigenvalue<Complex, SMALLEST_REAL> sorting(m_ritz_val.data(), m_nev);
                ind = sorting.index();
            }
                break;
            case SMALLEST_IMAG:
            {
                SortEigenvalue<Complex, SMALLEST_IMAG> sorting(m_ritz_val.data(), m_nev);
                ind = sorting.index();
            }
                break;
            default:
                throw std::invalid_argument("unsupported sorting rule");
        }

        ComplexVector new_ritz_val(m_ncv);
        ComplexMatrix new_ritz_vec(m_ncv, m_nev);
        BoolArray new_ritz_conv(m_nev);

        for(Index i = 0; i < m_nev; i++)
        {
            new_ritz_val[i] = m_ritz_val[ind[i]];
            new_ritz_vec.col(i).noalias() = m_ritz_vec.col(ind[i]);
            new_ritz_conv[i] = m_ritz_conv[ind[i]];
        }

        m_ritz_val.swap(new_ritz_val);
        m_ritz_vec.swap(new_ritz_vec);
        m_ritz_conv.swap(new_ritz_conv);
    }

public:
    /// \cond

    GenEigsBase(OpType* op, BOpType* Bop, Index nev, Index ncv) :
        m_op(op),
        m_n(m_op->rows()),
        m_nev(nev),
        m_ncv(ncv > m_n ? m_n : ncv),
        m_nmatop(0),
        m_niter(0),
        m_fac(ArnoldiOpType(op, Bop), m_ncv),
        m_info(NOT_COMPUTED),
        m_near_0(TypeTraits<Scalar>::min() * Scalar(10)),
        m_eps(Eigen::NumTraits<Scalar>::epsilon()),
        m_eps23(Eigen::numext::pow(m_eps, Scalar(2.0) / 3))
    {
        if(nev < 1 || nev > m_n - 2)
            throw std::invalid_argument("nev must satisfy 1 <= nev <= n - 2, n is the size of matrix");

        if(ncv < nev + 2 || ncv > m_n)
            throw std::invalid_argument("ncv must satisfy nev + 2 <= ncv <= n, n is the size of matrix");
    }

    ///
    /// Virtual destructor
    ///
    virtual ~GenEigsBase() {}

    /// \endcond

    ///
    /// Initializes the solver by providing an initial residual vector.
    ///
    /// \param init_resid Pointer to the initial residual vector.
    ///
    /// **Spectra** (and also **ARPACK**) uses an iterative algorithm
    /// to find eigenvalues. This function allows the user to provide the initial
    /// residual vector.
    ///
    void init(const Scalar* init_resid)
    {
        // Reset all matrices/vectors to zero
        m_ritz_val.resize(m_ncv);
        m_ritz_vec.resize(m_ncv, m_nev);
        m_ritz_est.resize(m_ncv);
        m_ritz_conv.resize(m_nev);

        m_ritz_val.setZero();
        m_ritz_vec.setZero();
        m_ritz_est.setZero();
        m_ritz_conv.setZero();

        m_nmatop = 0;
        m_niter = 0;

        // Initialize the Arnoldi factorization
        MapConstVec v0(init_resid, m_n);
        m_fac.init(v0, m_nmatop);
    }

    ///
    /// Initializes the solver by providing a random initial residual vector.
    ///
    /// This overloaded function generates a random initial residual vector
    /// (with a fixed random seed) for the algorithm. Elements in the vector
    /// follow independent Uniform(-0.5, 0.5) distribution.
    ///
    void init()
    {
        SimpleRandom<Scalar> rng(0);
        Vector init_resid = rng.random_vec(m_n);
        init(init_resid.data());
    }

    ///
    /// Conducts the major computation procedure.
    ///
    /// \param maxit      Maximum number of iterations allowed in the algorithm.
    /// \param tol        Precision parameter for the calculated eigenvalues.
    /// \param sort_rule  Rule to sort the eigenvalues and eigenvectors.
    ///                   Supported values are
    ///                   `Spectra::LARGEST_MAGN`, `Spectra::LARGEST_REAL`,
    ///                   `Spectra::LARGEST_IMAG`, `Spectra::SMALLEST_MAGN`,
    ///                   `Spectra::SMALLEST_REAL` and `Spectra::SMALLEST_IMAG`,
    ///                   for example `LARGEST_MAGN` indicates that eigenvalues
    ///                   with largest magnitude come first.
    ///                   Note that this argument is only used to
    ///                   **sort** the final result, and the **selection** rule
    ///                   (e.g. selecting the largest or smallest eigenvalues in the
    ///                   full spectrum) is specified by the template parameter
    ///                   `SelectionRule` of GenEigsSolver.
    ///
    /// \return Number of converged eigenvalues.
    ///
    Index compute(Index maxit = 1000, Scalar tol = 1e-10, int sort_rule = LARGEST_MAGN)
    {
        // The m-step Arnoldi factorization
        m_fac.factorize_from(1, m_ncv, m_nmatop);
        retrieve_ritzpair();
        // Restarting
        Index i, nconv = 0, nev_adj;
        for(i = 0; i < maxit; i++)
        {
            nconv = num_converged(tol);
            if(nconv >= m_nev)
                break;

            nev_adj = nev_adjusted(nconv);
            restart(nev_adj);
        }
        // Sorting results
        sort_ritzpair(sort_rule);

        m_niter += i + 1;
        m_info = (nconv >= m_nev) ? SUCCESSFUL : NOT_CONVERGING;

        return std::min(m_nev, nconv);
    }

    ///
    /// Returns the status of the computation.
    /// The full list of enumeration values can be found in \ref Enumerations.
    ///
    int info() const { return m_info; }

    ///
    /// Returns the number of iterations used in the computation.
    ///
    Index num_iterations() const { return m_niter; }

    ///
    /// Returns the number of matrix operations used in the computation.
    ///
    Index num_operations() const { return m_nmatop; }

    ///
    /// Returns the converged eigenvalues.
    ///
    /// \return A complex-valued vector containing the eigenvalues.
    /// Returned vector type will be `Eigen::Vector<std::complex<Scalar>, ...>`, depending on
    /// the template parameter `Scalar` defined.
    ///
    ComplexVector eigenvalues() const
    {
        const Index nconv = m_ritz_conv.cast<Index>().sum();
        ComplexVector res(nconv);

        if(!nconv)
            return res;

        Index j = 0;
        for(Index i = 0; i < m_nev; i++)
        {
            if(m_ritz_conv[i])
            {
                res[j] = m_ritz_val[i];
                j++;
            }
        }

        return res;
    }

    ///
    /// Returns the eigenvectors associated with the converged eigenvalues.
    ///
    /// \param nvec The number of eigenvectors to return.
    ///
    /// \return A complex-valued matrix containing the eigenvectors.
    /// Returned matrix type will be `Eigen::Matrix<std::complex<Scalar>, ...>`,
    /// depending on the template parameter `Scalar` defined.
    ///
    ComplexMatrix eigenvectors(Index nvec) const
    {
        const Index nconv = m_ritz_conv.cast<Index>().sum();
        nvec = std::min(nvec, nconv);
        ComplexMatrix res(m_n, nvec);

        if(!nvec)
            return res;

        ComplexMatrix ritz_vec_conv(m_ncv, nvec);
        Index j = 0;
        for(Index i = 0; i < m_nev && j < nvec; i++)
        {
            if(m_ritz_conv[i])
            {
                ritz_vec_conv.col(j).noalias() = m_ritz_vec.col(i);
                j++;
            }
        }

        res.noalias() = m_fac.matrix_V() * ritz_vec_conv;

        return res;
    }

    ///
    /// Returns all converged eigenvectors.
    ///
    ComplexMatrix eigenvectors() const
    {
        return eigenvectors(m_nev);
    }
};


} // namespace Spectra

#endif // GEN_EIGS_BASE_H
