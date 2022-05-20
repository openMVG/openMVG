// Copyright (C) 2020 Netherlands eScience Center <n.renauld@esciencecenter.nl>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SPECTRA_RITZ_PAIRS_H
#define SPECTRA_RITZ_PAIRS_H

#include <Eigen/Core>
#include <Eigen/Eigenvalues>

#include "../Util/SelectionRule.h"

namespace Spectra {

template <typename Scalar>
class SearchSpace;

/// This class handles the creation and manipulation of Ritz eigen pairs
/// for iterative eigensolvers such as Davidson, Jacobi-Davidson, etc.
template <typename Scalar>
class RitzPairs
{
private:
    using Index = Eigen::Index;
    using Matrix = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
    using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
    using Array = Eigen::Array<Scalar, Eigen::Dynamic, 1>;
    using BoolArray = Eigen::Array<bool, Eigen::Dynamic, 1>;

    Vector m_values;         // eigenvalues
    Matrix m_small_vectors;  // eigenvectors of the small problem, makes restart cheaper.
    Matrix m_vectors;        // Ritz (or harmonic Ritz) eigenvectors
    Matrix m_residues;       // residues of the pairs
    BoolArray m_root_converged;

public:
    RitzPairs() = default;

    /// Compute the eigen values/vectors
    ///
    /// \param search_space Instance of the class handling the search space
    /// \return Eigen::ComputationalInfo Whether small eigenvalue problem worked
    Eigen::ComputationInfo compute_eigen_pairs(const SearchSpace<Scalar>& search_space);

    /// Returns the size of the ritz eigen pairs
    ///
    /// \return Eigen::Index Number of pairs
    Index size() const { return m_values.size(); }

    /// Sort the eigen pairs according to the selection rule
    ///
    /// \param selection Sorting rule
    void sort(SortRule selection)
    {
        std::vector<Index> ind = argsort(selection, m_values);
        RitzPairs<Scalar> temp = *this;
        for (Index i = 0; i < size(); i++)
        {
            m_values[i] = temp.m_values[ind[i]];
            m_vectors.col(i) = temp.m_vectors.col(ind[i]);
            m_residues.col(i) = temp.m_residues.col(ind[i]);
            m_small_vectors.col(i) = temp.m_small_vectors.col(ind[i]);
        }
    }

    /// Checks if the algorithm has converged and updates root_converged
    ///
    /// \param tol Tolerance for convergence
    /// \param number_eigenvalue Number of request eigenvalues
    /// \return bool true if all eigenvalues are converged
    bool check_convergence(Scalar tol, Index number_eigenvalues)
    {
        const Array norms = m_residues.colwise().norm();
        bool converged = true;
        m_root_converged = BoolArray::Zero(norms.size());
        for (Index j = 0; j < norms.size(); j++)
        {
            m_root_converged[j] = (norms[j] < tol);
            if (j < number_eigenvalues)
            {
                converged &= (norms[j] < tol);
            }
        }
        return converged;
    }

    const Matrix& ritz_vectors() const { return m_vectors; }
    const Vector& ritz_values() const { return m_values; }
    const Matrix& small_ritz_vectors() const { return m_small_vectors; }
    const Matrix& residues() const { return m_residues; }
    const BoolArray& converged_eigenvalues() const { return m_root_converged; }
};

}  // namespace Spectra

#include "SearchSpace.h"

namespace Spectra {

/// Creates the small space matrix and computes its eigen pairs
/// Also computes the ritz vectors and residues
///
/// \param search_space Instance of the SearchSpace class
template <typename Scalar>
Eigen::ComputationInfo RitzPairs<Scalar>::compute_eigen_pairs(const SearchSpace<Scalar>& search_space)
{
    const Matrix& basis_vectors = search_space.basis_vectors();
    const Matrix& op_basis_prod = search_space.operator_basis_product();

    // Form the small eigenvalue
    Matrix small_matrix = basis_vectors.transpose() * op_basis_prod;

    // Small eigenvalue problem
    Eigen::SelfAdjointEigenSolver<Matrix> eigen_solver(small_matrix);
    m_values = eigen_solver.eigenvalues();
    m_small_vectors = eigen_solver.eigenvectors();

    // Ritz vectors
    m_vectors = basis_vectors * m_small_vectors;

    // Residues
    m_residues = op_basis_prod * m_small_vectors - m_vectors * m_values.asDiagonal();
    return eigen_solver.info();
}

}  // namespace Spectra

#endif  // SPECTRA_RITZ_PAIRS_H
