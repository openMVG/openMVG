// Copyright (C) 2020 Netherlands eScience Center <f.zapata@esciencecenter.nl>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SPECTRA_DAVIDSON_SYM_EIGS_SOLVER_H
#define SPECTRA_DAVIDSON_SYM_EIGS_SOLVER_H

#include <Eigen/Core>

#include "JDSymEigsBase.h"
#include "Util/SelectionRule.h"

namespace Spectra {

///
/// \ingroup EigenSolver
///
/// This class implement the DPR correction for the Davidson algorithms.
/// The algorithms in the Davidson family only differ in how the correction
/// vectors are computed and optionally in the initial orthogonal basis set.
///
/// the DPR correction compute the new correction vector using the following expression:
/// \f[ correction = -(\boldsymbol{D} - \rho \boldsymbol{I})^{-1} \boldsymbol{r} \f]
/// where
/// \f$D\f$ is the diagonal of the target matrix, \f$\rho\f$ the Ritz eigenvalue,
/// \f$I\f$ the identity matrix and \f$r\f$ the residue vector.
///
template <typename OpType>
class DavidsonSymEigsSolver : public JDSymEigsBase<DavidsonSymEigsSolver<OpType>, OpType>
{
private:
    using Index = Eigen::Index;
    using Scalar = typename OpType::Scalar;
    using Matrix = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
    using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

    Vector m_diagonal;

public:
    DavidsonSymEigsSolver(OpType& op, Index nev, Index nvec_init, Index nvec_max) :
        JDSymEigsBase<DavidsonSymEigsSolver<OpType>, OpType>(op, nev, nvec_init, nvec_max)
    {
        m_diagonal.resize(this->m_matrix_operator.rows());
        for (Index i = 0; i < op.rows(); i++)
        {
            m_diagonal(i) = op(i, i);
        }
    }

    DavidsonSymEigsSolver(OpType& op, Index nev) :
        DavidsonSymEigsSolver(op, nev, 2 * nev, 10 * nev) {}

    /// Create initial search space based on the diagonal
    /// and the spectrum'target (highest or lowest)
    ///
    /// \param selection Spectrum section to target (e.g. lowest, etc.)
    /// \return Matrix with the initial orthonormal basis
    Matrix setup_initial_search_space(SortRule selection) const
    {
        std::vector<Eigen::Index> indices_sorted = argsort(selection, m_diagonal);

        Matrix initial_basis = Matrix::Zero(this->m_matrix_operator.rows(), this->m_initial_search_space_size);

        for (Index k = 0; k < this->m_initial_search_space_size; k++)
        {
            Index row = indices_sorted[k];
            initial_basis(row, k) = 1.0;
        }
        return initial_basis;
    }

    /// Compute the corrections using the DPR method.
    ///
    /// \return New correction vectors.
    Matrix calculate_correction_vector() const
    {
        const Matrix& residues = this->m_ritz_pairs.residues();
        const Vector& eigvals = this->m_ritz_pairs.ritz_values();
        Matrix correction = Matrix::Zero(this->m_matrix_operator.rows(), this->m_correction_size);
        for (Index k = 0; k < this->m_correction_size; k++)
        {
            Vector tmp = eigvals(k) - m_diagonal.array();
            correction.col(k) = residues.col(k).array() / tmp.array();
        }
        return correction;
    }
};

}  // namespace Spectra

#endif  // SPECTRA_DAVIDSON_SYM_EIGS_SOLVER_H
