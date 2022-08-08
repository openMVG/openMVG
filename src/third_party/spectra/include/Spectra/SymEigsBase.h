// Copyright (C) 2018-2022 Yixuan Qiu <yixuan.qiu@cos.name>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SPECTRA_SYM_EIGS_BASE_H
#define SPECTRA_SYM_EIGS_BASE_H

#include <Eigen/Core>
#include <vector>     // std::vector
#include <cmath>      // std::abs, std::pow
#include <algorithm>  // std::min
#include <stdexcept>  // std::invalid_argument
#include <utility>    // std::move

#include "Util/Version.h"
#include "Util/TypeTraits.h"
#include "Util/SelectionRule.h"
#include "Util/CompInfo.h"
#include "Util/SimpleRandom.h"
#include "MatOp/internal/ArnoldiOp.h"
#include "LinAlg/UpperHessenbergQR.h"
#include "LinAlg/TridiagEigen.h"
#include "LinAlg/Lanczos.h"

namespace Spectra {

///
/// \defgroup EigenSolver Eigen Solvers
///
/// Eigen solvers for different types of problems.
///

///
/// \ingroup EigenSolver
///
/// This is the base class for symmetric eigen solvers, mainly for internal use.
/// It is kept here to provide the documentation for member functions of concrete eigen solvers
/// such as SymEigsSolver and SymEigsShiftSolver.
///
template <typename OpType, typename BOpType>
class SymEigsBase
{
private:
    using Scalar = typename OpType::Scalar;
    using Index = Eigen::Index;
    using Matrix = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
    using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
    using Array = Eigen::Array<Scalar, Eigen::Dynamic, 1>;
    using BoolArray = Eigen::Array<bool, Eigen::Dynamic, 1>;
    using MapMat = Eigen::Map<Matrix>;
    using MapVec = Eigen::Map<Vector>;
    using MapConstVec = Eigen::Map<const Vector>;

    using ArnoldiOpType = ArnoldiOp<Scalar, OpType, BOpType>;
    using LanczosFac = Lanczos<Scalar, ArnoldiOpType>;

protected:
    // clang-format off

    // In SymEigsSolver and SymEigsShiftSolver, the A operator is an lvalue provided by
    // the user. In SymGEigsSolver, the A operator is an rvalue. To avoid copying objects,
    // we use the following scheme:
    // 1. If the op parameter in the constructor is an lvalue, make m_op a const reference to op
    // 2. If op is an rvalue, move op to m_op_container, and then make m_op a const
    //    reference to m_op_container[0]
    std::vector<OpType> m_op_container;
    const OpType& m_op;         // matrix operator for A
    const Index   m_n;          // dimension of matrix A
    const Index   m_nev;        // number of eigenvalues requested
    const Index   m_ncv;        // dimension of Krylov subspace in the Lanczos method
    Index         m_nmatop;     // number of matrix operations called
    Index         m_niter;      // number of restarting iterations

    LanczosFac    m_fac;        // Lanczos factorization
    Vector        m_ritz_val;   // Ritz values

private:
    Matrix        m_ritz_vec;   // Ritz vectors
    Vector        m_ritz_est;   // last row of m_ritz_vec, also called the Ritz estimates
    BoolArray     m_ritz_conv;  // indicator of the convergence of Ritz values
    CompInfo      m_info;       // status of the computation
    // clang-format on

    // Move rvalue object to the container
    static std::vector<OpType> create_op_container(OpType&& rval)
    {
        std::vector<OpType> container;
        container.emplace_back(std::move(rval));
        return container;
    }

    // Implicitly restarted Lanczos factorization
    void restart(Index k, SortRule selection)
    {
        using std::abs;

        if (k >= m_ncv)
            return;

        TridiagQR<Scalar> decomp(m_ncv);
        Matrix Q = Matrix::Identity(m_ncv, m_ncv);

        // Apply large shifts first
        const int nshift = m_ncv - k;
        Vector shifts = m_ritz_val.tail(nshift);
        std::sort(shifts.data(), shifts.data() + nshift, [](const Scalar& v1, const Scalar& v2) { return abs(v1) > abs(v2); });

        for (Index i = 0; i < nshift; i++)
        {
            // QR decomposition of H-mu*I, mu is the shift
            decomp.compute(m_fac.matrix_H(), shifts[i]);

            // Q -> Q * Qi
            decomp.apply_YQ(Q);
            // H -> Q'HQ
            // Since QR = H - mu * I, we have H = QR + mu * I
            // and therefore Q'HQ = RQ + mu * I
            m_fac.compress_H(decomp);
        }

        m_fac.compress_V(Q);
        m_fac.factorize_from(k, m_ncv, m_nmatop);

        retrieve_ritzpair(selection);
    }

    // Calculates the number of converged Ritz values
    Index num_converged(const Scalar& tol)
    {
        using std::pow;

        // The machine precision, ~= 1e-16 for the "double" type
        constexpr Scalar eps = TypeTraits<Scalar>::epsilon();
        // std::pow() is not constexpr, so we do not declare eps23 to be constexpr
        // But most compilers should be able to compute eps23 at compile time
        const Scalar eps23 = pow(eps, Scalar(2) / 3);

        // thresh = tol * max(eps23, abs(theta)), theta for Ritz value
        Array thresh = tol * m_ritz_val.head(m_nev).array().abs().max(eps23);
        Array resid = m_ritz_est.head(m_nev).array().abs() * m_fac.f_norm();
        // Converged "wanted" Ritz values
        m_ritz_conv = (resid < thresh);

        return m_ritz_conv.count();
    }

    // Returns the adjusted nev for restarting
    Index nev_adjusted(Index nconv)
    {
        using std::abs;

        // A very small value, but 1.0 / near_0 does not overflow
        // ~= 1e-307 for the "double" type
        constexpr Scalar near_0 = TypeTraits<Scalar>::min() * Scalar(10);

        Index nev_new = m_nev;
        for (Index i = m_nev; i < m_ncv; i++)
            if (abs(m_ritz_est[i]) < near_0)
                nev_new++;

        // Adjust nev_new, according to dsaup2.f line 677~684 in ARPACK
        nev_new += (std::min)(nconv, (m_ncv - nev_new) / 2);
        if (nev_new == 1 && m_ncv >= 6)
            nev_new = m_ncv / 2;
        else if (nev_new == 1 && m_ncv > 2)
            nev_new = 2;

        if (nev_new > m_ncv - 1)
            nev_new = m_ncv - 1;

        return nev_new;
    }

    // Retrieves and sorts Ritz values and Ritz vectors
    void retrieve_ritzpair(SortRule selection)
    {
        TridiagEigen<Scalar> decomp(m_fac.matrix_H());
        const Vector& evals = decomp.eigenvalues();
        const Matrix& evecs = decomp.eigenvectors();

        // Sort Ritz values and put the wanted ones at the beginning
        std::vector<Index> ind = argsort(selection, evals, m_ncv);

        // Copy the Ritz values and vectors to m_ritz_val and m_ritz_vec, respectively
        for (Index i = 0; i < m_ncv; i++)
        {
            m_ritz_val[i] = evals[ind[i]];
            m_ritz_est[i] = evecs(m_ncv - 1, ind[i]);
        }
        for (Index i = 0; i < m_nev; i++)
        {
            m_ritz_vec.col(i).noalias() = evecs.col(ind[i]);
        }
    }

protected:
    // Sorts the first nev Ritz pairs in the specified order
    // This is used to return the final results
    virtual void sort_ritzpair(SortRule sort_rule)
    {
        if ((sort_rule != SortRule::LargestAlge) && (sort_rule != SortRule::LargestMagn) &&
            (sort_rule != SortRule::SmallestAlge) && (sort_rule != SortRule::SmallestMagn))
            throw std::invalid_argument("unsupported sorting rule");

        std::vector<Index> ind = argsort(sort_rule, m_ritz_val, m_nev);

        Vector new_ritz_val(m_ncv);
        Matrix new_ritz_vec(m_ncv, m_nev);
        BoolArray new_ritz_conv(m_nev);

        for (Index i = 0; i < m_nev; i++)
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

    // If op is an lvalue
    SymEigsBase(OpType& op, const BOpType& Bop, Index nev, Index ncv) :
        m_op(op),
        m_n(op.rows()),
        m_nev(nev),
        m_ncv(ncv > m_n ? m_n : ncv),
        m_nmatop(0),
        m_niter(0),
        m_fac(ArnoldiOpType(op, Bop), m_ncv),
        m_info(CompInfo::NotComputed)
    {
        if (nev < 1 || nev > m_n - 1)
            throw std::invalid_argument("nev must satisfy 1 <= nev <= n - 1, n is the size of matrix");

        if (ncv <= nev || ncv > m_n)
            throw std::invalid_argument("ncv must satisfy nev < ncv <= n, n is the size of matrix");
    }

    // If op is an rvalue
    SymEigsBase(OpType&& op, const BOpType& Bop, Index nev, Index ncv) :
        m_op_container(create_op_container(std::move(op))),
        m_op(m_op_container.front()),
        m_n(m_op.rows()),
        m_nev(nev),
        m_ncv(ncv > m_n ? m_n : ncv),
        m_nmatop(0),
        m_niter(0),
        m_fac(ArnoldiOpType(m_op, Bop), m_ncv),
        m_info(CompInfo::NotComputed)
    {
        if (nev < 1 || nev > m_n - 1)
            throw std::invalid_argument("nev must satisfy 1 <= nev <= n - 1, n is the size of matrix");

        if (ncv <= nev || ncv > m_n)
            throw std::invalid_argument("ncv must satisfy nev < ncv <= n, n is the size of matrix");
    }

    ///
    /// Virtual destructor
    ///
    virtual ~SymEigsBase() {}

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

        // Initialize the Lanczos factorization
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
    /// \param selection  An enumeration value indicating the selection rule of
    ///                   the requested eigenvalues, for example `SortRule::LargestMagn`
    ///                   to retrieve eigenvalues with the largest magnitude.
    ///                   The full list of enumeration values can be found in
    ///                   \ref Enumerations.
    /// \param maxit      Maximum number of iterations allowed in the algorithm.
    /// \param tol        Precision parameter for the calculated eigenvalues.
    /// \param sorting    Rule to sort the eigenvalues and eigenvectors.
    ///                   Supported values are
    ///                   `SortRule::LargestAlge`, `SortRule::LargestMagn`,
    ///                   `SortRule::SmallestAlge`, and `SortRule::SmallestMagn`.
    ///                   For example, `SortRule::LargestAlge` indicates that largest eigenvalues
    ///                   come first. Note that this argument is only used to
    ///                   **sort** the final result, and the **selection** rule
    ///                   (e.g. selecting the largest or smallest eigenvalues in the
    ///                   full spectrum) is specified by the parameter `selection`.
    ///
    /// \return Number of converged eigenvalues.
    ///
    Index compute(SortRule selection = SortRule::LargestMagn, Index maxit = 1000,
                  Scalar tol = 1e-10, SortRule sorting = SortRule::LargestAlge)
    {
        // The m-step Lanczos factorization
        m_fac.factorize_from(1, m_ncv, m_nmatop);
        retrieve_ritzpair(selection);
        // Restarting
        Index i, nconv = 0, nev_adj;
        for (i = 0; i < maxit; i++)
        {
            nconv = num_converged(tol);
            if (nconv >= m_nev)
                break;

            nev_adj = nev_adjusted(nconv);
            restart(nev_adj, selection);
        }
        // Sorting results
        sort_ritzpair(sorting);

        m_niter += i + 1;
        m_info = (nconv >= m_nev) ? CompInfo::Successful : CompInfo::NotConverging;

        return (std::min)(m_nev, nconv);
    }

    ///
    /// Returns the status of the computation.
    /// The full list of enumeration values can be found in \ref Enumerations.
    ///
    CompInfo info() const { return m_info; }

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
    /// \return A vector containing the eigenvalues.
    /// Returned vector type will be `Eigen::Vector<Scalar, ...>`, depending on
    /// the template parameter `Scalar` defined.
    ///
    Vector eigenvalues() const
    {
        const Index nconv = m_ritz_conv.count();
        Vector res(nconv);

        if (!nconv)
            return res;

        Index j = 0;
        for (Index i = 0; i < m_nev; i++)
        {
            if (m_ritz_conv[i])
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
    /// \return A matrix containing the eigenvectors.
    /// Returned matrix type will be `Eigen::Matrix<Scalar, ...>`,
    /// depending on the template parameter `Scalar` defined.
    ///
    virtual Matrix eigenvectors(Index nvec) const
    {
        const Index nconv = m_ritz_conv.count();
        nvec = (std::min)(nvec, nconv);
        Matrix res(m_n, nvec);

        if (!nvec)
            return res;

        Matrix ritz_vec_conv(m_ncv, nvec);
        Index j = 0;
        for (Index i = 0; i < m_nev && j < nvec; i++)
        {
            if (m_ritz_conv[i])
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
    virtual Matrix eigenvectors() const
    {
        return eigenvectors(m_nev);
    }
};

}  // namespace Spectra

#endif  // SPECTRA_SYM_EIGS_BASE_H
