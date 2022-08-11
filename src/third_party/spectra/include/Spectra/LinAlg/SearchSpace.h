// Copyright (C) 2020 Netherlands eScience Center <n.renauld@esciencecenter.nl>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SPECTRA_SEARCH_SPACE_H
#define SPECTRA_SEARCH_SPACE_H

#include <Eigen/Core>

#include "RitzPairs.h"
#include "Orthogonalization.h"

namespace Spectra {

/// This class handles the creation and manipulation of the search space
/// for iterative eigensolvers such as Davidson, Jacobi-Davidson, etc.
template <typename Scalar>
class SearchSpace
{
private:
    using Index = Eigen::Index;
    using Matrix = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;

    Matrix m_basis_vectors;
    Matrix m_op_basis_product;

    /// Append new vector to the basis
    ///
    /// \param new_vect Matrix of new correction vectors
    void append_new_vectors_to_basis(const Matrix& new_vect)
    {
        Index num_update = new_vect.cols();
        m_basis_vectors.conservativeResize(Eigen::NoChange, m_basis_vectors.cols() + num_update);
        m_basis_vectors.rightCols(num_update).noalias() = new_vect;
    }

public:
    SearchSpace() = default;

    /// Returns the current size of the search space
    Index size() const { return m_basis_vectors.cols(); }

    void initialize_search_space(const Eigen::Ref<const Matrix>& initial_vectors)
    {
        m_basis_vectors = initial_vectors;
        m_op_basis_product = Matrix(initial_vectors.rows(), 0);
    }

    /// Updates the matrix formed by the operator applied to the search space
    /// after the addition of new vectors in the search space. Only the product
    /// of the operator with the new vectors is computed and the result is appended
    /// to the op_basis_product member variable
    ///
    /// \param OpType Operator representing the matrix
    template <typename OpType>
    void update_operator_basis_product(OpType& op)
    {
        Index nvec = m_basis_vectors.cols() - m_op_basis_product.cols();
        m_op_basis_product.conservativeResize(Eigen::NoChange, m_basis_vectors.cols());
        m_op_basis_product.rightCols(nvec).noalias() = op * m_basis_vectors.rightCols(nvec);
    }

    /// Restart the search space by reducing the basis vector to the last
    /// Ritz eigenvector
    ///
    /// \param ritz_pair Instance of a RitzPair class
    /// \param size Size of the restart
    void restart(const RitzPairs<Scalar>& ritz_pairs, Index size)
    {
        m_basis_vectors = ritz_pairs.ritz_vectors().leftCols(size);
        m_op_basis_product = m_op_basis_product * ritz_pairs.small_ritz_vectors().leftCols(size);
    }

    /// Append new vectors to the search space and
    /// orthogonalize the resulting matrix
    ///
    /// \param new_vect Matrix of new correction vectors
    void extend_basis(const Matrix& new_vect)
    {
        Index left_cols_to_skip = size();
        append_new_vectors_to_basis(new_vect);
        twice_is_enough_orthogonalisation(m_basis_vectors, left_cols_to_skip);
    }

    /// Returns the basis vectors
    const Matrix& basis_vectors() const { return m_basis_vectors; }

    /// Returns the operator applied to basis vector
    const Matrix& operator_basis_product() const { return m_op_basis_product; }
};

}  // namespace Spectra

#endif  // SPECTRA_SEARCH_SPACE_H
