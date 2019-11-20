// Copyright (C) 2016-2019 Yixuan Qiu <yixuan.qiu@cos.name>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef DENSE_CHOLESKY_H
#define DENSE_CHOLESKY_H

#include <Eigen/Core>
#include <Eigen/Cholesky>
#include <stdexcept>
#include "../Util/CompInfo.h"

namespace Spectra {


///
/// \ingroup MatOp
///
/// This class defines the operations related to Cholesky decomposition on a
/// positive definite matrix, \f$B=LL'\f$, where \f$L\f$ is a lower triangular
/// matrix. It is mainly used in the SymGEigsSolver generalized eigen solver
/// in the Cholesky decomposition mode.
///
template <typename Scalar, int Uplo = Eigen::Lower>
class DenseCholesky
{
private:
    typedef Eigen::Index Index;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Matrix;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;
    typedef Eigen::Map<const Matrix> MapConstMat;
    typedef Eigen::Map<const Vector> MapConstVec;
    typedef Eigen::Map<Vector> MapVec;
    typedef const Eigen::Ref<const Matrix> ConstGenericMatrix;

    const Index m_n;
    Eigen::LLT<Matrix, Uplo> m_decomp;
    int m_info;  // status of the decomposition

public:
    ///
    /// Constructor to create the matrix operation object.
    ///
    /// \param mat An **Eigen** matrix object, whose type can be
    /// `Eigen::Matrix<Scalar, ...>` (e.g. `Eigen::MatrixXd` and
    /// `Eigen::MatrixXf`), or its mapped version
    /// (e.g. `Eigen::Map<Eigen::MatrixXd>`).
    ///
    DenseCholesky(ConstGenericMatrix& mat) :
        m_n(mat.rows()), m_info(NOT_COMPUTED)
    {
        if(mat.rows() != mat.cols())
            throw std::invalid_argument("DenseCholesky: matrix must be square");

        m_decomp.compute(mat);
        m_info = (m_decomp.info() == Eigen::Success) ?
                 SUCCESSFUL :
                 NUMERICAL_ISSUE;
    }

    ///
    /// Returns the number of rows of the underlying matrix.
    ///
    Index rows() const { return m_n; }
    ///
    /// Returns the number of columns of the underlying matrix.
    ///
    Index cols() const { return m_n; }

    ///
    /// Returns the status of the computation.
    /// The full list of enumeration values can be found in \ref Enumerations.
    ///
    int info() const { return m_info; }

    ///
    /// Performs the lower triangular solving operation \f$y=L^{-1}x\f$.
    ///
    /// \param x_in  Pointer to the \f$x\f$ vector.
    /// \param y_out Pointer to the \f$y\f$ vector.
    ///
    // y_out = inv(L) * x_in
    void lower_triangular_solve(const Scalar* x_in, Scalar* y_out) const
    {
        MapConstVec x(x_in,  m_n);
        MapVec      y(y_out, m_n);
        y.noalias() = m_decomp.matrixL().solve(x);
    }

    ///
    /// Performs the upper triangular solving operation \f$y=(L')^{-1}x\f$.
    ///
    /// \param x_in  Pointer to the \f$x\f$ vector.
    /// \param y_out Pointer to the \f$y\f$ vector.
    ///
    // y_out = inv(L') * x_in
    void upper_triangular_solve(const Scalar* x_in, Scalar* y_out) const
    {
        MapConstVec x(x_in,  m_n);
        MapVec      y(y_out, m_n);
        y.noalias() = m_decomp.matrixU().solve(x);
    }
};


} // namespace Spectra

#endif // DENSE_CHOLESKY_H
