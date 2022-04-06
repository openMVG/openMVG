// Copyright (C) 2020 Netherlands eScience Center <f.zapata@esciencecenter.nl>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SPECTRA_ORTHOGONALIZATION_H
#define SPECTRA_ORTHOGONALIZATION_H

#include <Eigen/Core>
#include <Eigen/QR>

namespace Spectra {

/// Check if the number of columns to skip is
/// larger than 0 but smaller than the total number
/// of columns of the matrix
/// \param in_output Matrix to be orthogonalized
/// \param left_cols_to_skip Number of left columns to be left untouched
template <typename Matrix>
void assert_left_cols_to_skip(Matrix& in_output, Eigen::Index left_cols_to_skip)
{
    assert(in_output.cols() > left_cols_to_skip && "left_cols_to_skip is larger than columns of matrix");
    assert(left_cols_to_skip >= 0 && "left_cols_to_skip is negative");
}

/// If the the number of columns to skip is null,
/// normalize the first column and set left_cols_to_skip=1
/// \param in_output Matrix to be orthogonalized
/// \param left_cols_to_skip Number of left columns to be left untouched
/// \return Actual number of left columns to skip
template <typename Matrix>
Eigen::Index treat_first_col(Matrix& in_output, Eigen::Index left_cols_to_skip)
{
    if (left_cols_to_skip == 0)
    {
        in_output.col(0).normalize();
        left_cols_to_skip = 1;
    }
    return left_cols_to_skip;
}

/// Orthogonalize the in_output matrix using a QR decomposition
/// \param in_output Matrix to be orthogonalized
template <typename Matrix>
void QR_orthogonalisation(Matrix& in_output)
{
    using InternalMatrix = Eigen::Matrix<typename Matrix::Scalar, Eigen::Dynamic, Eigen::Dynamic>;
    Eigen::Index nrows = in_output.rows();
    Eigen::Index ncols = in_output.cols();
    ncols = (std::min)(nrows, ncols);
    InternalMatrix I = InternalMatrix::Identity(nrows, ncols);
    Eigen::HouseholderQR<Matrix> qr(in_output);
    in_output.leftCols(ncols).noalias() = qr.householderQ() * I;
}

/// Orthogonalize the in_output matrix using a modified Gram Schmidt process
/// \param in_output matrix to be orthogonalized
/// \param left_cols_to_skip Number of left columns to be left untouched
template <typename Matrix>
void MGS_orthogonalisation(Matrix& in_output, Eigen::Index left_cols_to_skip = 0)
{
    assert_left_cols_to_skip(in_output, left_cols_to_skip);
    left_cols_to_skip = treat_first_col(in_output, left_cols_to_skip);

    for (Eigen::Index k = left_cols_to_skip; k < in_output.cols(); ++k)
    {
        for (Eigen::Index j = 0; j < k; j++)
        {
            in_output.col(k) -= in_output.col(j).dot(in_output.col(k)) * in_output.col(j);
        }
        in_output.col(k).normalize();
    }
}

/// Orthogonalize the in_output matrix using a Gram Schmidt process
/// \param in_output matrix to be orthogonalized
/// \param left_cols_to_skip Number of left columns to be left untouched
template <typename Matrix>
void GS_orthogonalisation(Matrix& in_output, Eigen::Index left_cols_to_skip = 0)
{
    assert_left_cols_to_skip(in_output, left_cols_to_skip);
    left_cols_to_skip = treat_first_col(in_output, left_cols_to_skip);

    for (Eigen::Index j = left_cols_to_skip; j < in_output.cols(); ++j)
    {
        in_output.col(j) -= in_output.leftCols(j) * (in_output.leftCols(j).transpose() * in_output.col(j));
        in_output.col(j).normalize();
    }
}

/// Orthogonalize the subspace spanned by right columns of in_output
/// against the subspace spanned by left columns
/// It assumes that the left columns are already orthogonal and normalized,
/// and it does not orthogonalize the left columns against each other
/// \param in_output Matrix to be orthogonalized
/// \param left_cols_to_skip Number of left columns to be left untouched
template <typename Matrix>
void subspace_orthogonalisation(Matrix& in_output, Eigen::Index left_cols_to_skip)
{
    assert_left_cols_to_skip(in_output, left_cols_to_skip);
    if (left_cols_to_skip == 0)
    {
        return;
    }

    Eigen::Index right_cols_to_ortho = in_output.cols() - left_cols_to_skip;
    in_output.rightCols(right_cols_to_ortho) -= in_output.leftCols(left_cols_to_skip) *
        (in_output.leftCols(left_cols_to_skip).transpose() * in_output.rightCols(right_cols_to_ortho));
}

/// Orthogonalize the in_output matrix using a Jens process
/// The subspace spanned by right columns are first orthogonalized
/// agains the left columns, and then a QR decomposition is applied on the right columns
/// to make them orthogonalized agains each other
/// \param in_output Matrix to be orthogonalized
/// \param left_cols_to_skip Number of left columns to be left untouched
template <typename Matrix>
void JensWehner_orthogonalisation(Matrix& in_output, Eigen::Index left_cols_to_skip = 0)
{
    assert_left_cols_to_skip(in_output, left_cols_to_skip);

    Eigen::Index right_cols_to_ortho = in_output.cols() - left_cols_to_skip;
    subspace_orthogonalisation(in_output, left_cols_to_skip);
    Eigen::Ref<Matrix> right_cols = in_output.rightCols(right_cols_to_ortho);
    QR_orthogonalisation(right_cols);
}

/// Orthogonalize the in_output matrix using a twice-is-enough Jens process
/// \param in_output Matrix to be orthogonalized
/// \param left_cols_to_skip Number of left columns to be left untouched
template <typename Matrix>
void twice_is_enough_orthogonalisation(Matrix& in_output, Eigen::Index left_cols_to_skip = 0)
{
    JensWehner_orthogonalisation(in_output, left_cols_to_skip);
    JensWehner_orthogonalisation(in_output, left_cols_to_skip);
}

}  // namespace Spectra

#endif  // SPECTRA_ORTHOGONALIZATION_H
