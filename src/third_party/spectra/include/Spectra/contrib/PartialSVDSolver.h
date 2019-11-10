// Copyright (C) 2018 Yixuan Qiu <yixuan.qiu@cos.name>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef PARTIAL_SVD_SOLVER_H
#define PARTIAL_SVD_SOLVER_H

#include <Eigen/Core>
#include "../SymEigsSolver.h"


namespace Spectra {


// Abstract class for matrix operation
template <typename Scalar>
class SVDMatOp
{
public:
    virtual int rows() const = 0;
    virtual int cols() const = 0;

    // y_out = A' * A * x_in or y_out = A * A' * x_in
    virtual void perform_op(const Scalar* x_in, Scalar* y_out) = 0;

    virtual ~SVDMatOp() {}
};

// Operation of a tall matrix in SVD
// We compute the eigenvalues of A' * A
// MatrixType is either Eigen::Matrix<Scalar, ...> or Eigen::SparseMatrix<Scalar, ...>
template <typename Scalar, typename MatrixType>
class SVDTallMatOp: public SVDMatOp<Scalar>
{
private:
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;
    typedef Eigen::Map<const Vector> MapConstVec;
    typedef Eigen::Map<Vector> MapVec;
    typedef const Eigen::Ref<const MatrixType> ConstGenericMatrix;

    ConstGenericMatrix m_mat;
    const int m_dim;
    Vector m_cache;

public:
    // Constructor
    SVDTallMatOp(ConstGenericMatrix& mat) :
        m_mat(mat),
        m_dim(std::min(mat.rows(), mat.cols())),
        m_cache(mat.rows())
    {}

    // These are the rows and columns of A' * A
    int rows() const { return m_dim; }
    int cols() const { return m_dim; }

    // y_out = A' * A * x_in
    void perform_op(const Scalar* x_in, Scalar* y_out)
    {
        MapConstVec x(x_in,  m_mat.cols());
        MapVec      y(y_out, m_mat.cols());
        m_cache.noalias() = m_mat * x;
        y.noalias() = m_mat.transpose() * m_cache;
    }
};

// Operation of a wide matrix in SVD
// We compute the eigenvalues of A * A'
// MatrixType is either Eigen::Matrix<Scalar, ...> or Eigen::SparseMatrix<Scalar, ...>
template <typename Scalar, typename MatrixType>
class SVDWideMatOp: public SVDMatOp<Scalar>
{
private:
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;
    typedef Eigen::Map<const Vector> MapConstVec;
    typedef Eigen::Map<Vector> MapVec;
    typedef const Eigen::Ref<const MatrixType> ConstGenericMatrix;

    ConstGenericMatrix m_mat;
    const int m_dim;
    Vector m_cache;

public:
    // Constructor
    SVDWideMatOp(ConstGenericMatrix& mat) :
        m_mat(mat),
        m_dim(std::min(mat.rows(), mat.cols())),
        m_cache(mat.cols())
    {}

    // These are the rows and columns of A * A'
    int rows() const { return m_dim; }
    int cols() const { return m_dim; }

    // y_out = A * A' * x_in
    void perform_op(const Scalar* x_in, Scalar* y_out)
    {
        MapConstVec x(x_in,  m_mat.rows());
        MapVec      y(y_out, m_mat.rows());
        m_cache.noalias() = m_mat.transpose() * x;
        y.noalias() = m_mat * m_cache;
    }
};

// Partial SVD solver
// MatrixType is either Eigen::Matrix<Scalar, ...> or Eigen::SparseMatrix<Scalar, ...>
template < typename Scalar = double,
           typename MatrixType = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> >
class PartialSVDSolver
{
private:
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Matrix;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;
    typedef const Eigen::Ref<const MatrixType> ConstGenericMatrix;

    ConstGenericMatrix m_mat;
    const int m_m;
    const int m_n;
    SVDMatOp<Scalar>* m_op;
    SymEigsSolver< Scalar, LARGEST_ALGE, SVDMatOp<Scalar> >* m_eigs;
    int m_nconv;
    Matrix m_evecs;

public:
    // Constructor
    PartialSVDSolver(ConstGenericMatrix& mat, int ncomp, int ncv) :
        m_mat(mat), m_m(mat.rows()), m_n(mat.cols()), m_evecs(0, 0)
    {
        // Determine the matrix type, tall or wide
        if(m_m > m_n)
        {
            m_op = new SVDTallMatOp<Scalar, MatrixType>(mat);
        } else {
            m_op = new SVDWideMatOp<Scalar, MatrixType>(mat);
        }

        // Solver object
        m_eigs = new SymEigsSolver< Scalar, LARGEST_ALGE, SVDMatOp<Scalar> >(m_op, ncomp, ncv);
    }

    // Destructor
    virtual ~PartialSVDSolver()
    {
        delete m_eigs;
        delete m_op;
    }

    // Computation
    int compute(int maxit = 1000, Scalar tol = 1e-10)
    {
        m_eigs->init();
        m_nconv = m_eigs->compute(maxit, tol);

        return m_nconv;
    }

    // The converged singular values
    Vector singular_values() const
    {
        Vector svals = m_eigs->eigenvalues().cwiseSqrt();

        return svals;
    }

    // The converged left singular vectors
    Matrix matrix_U(int nu)
    {
        if(m_evecs.cols() < 1)
        {
            m_evecs = m_eigs->eigenvectors();
        }
        nu = std::min(nu, m_nconv);
        if(m_m <= m_n)
        {
            return m_evecs.leftCols(nu);
        }

        return m_mat * (m_evecs.leftCols(nu).array().rowwise() / m_eigs->eigenvalues().head(nu).transpose().array().sqrt()).matrix();
    }

    // The converged right singular vectors
    Matrix matrix_V(int nv)
    {
        if(m_evecs.cols() < 1)
        {
            m_evecs = m_eigs->eigenvectors();
        }
        nv = std::min(nv, m_nconv);
        if(m_m > m_n)
        {
            return m_evecs.leftCols(nv);
        }

        return m_mat.transpose() * (m_evecs.leftCols(nv).array().rowwise() / m_eigs->eigenvalues().head(nv).transpose().array().sqrt()).matrix();
    }
};


} // namespace Spectra

#endif // PARTIAL_SVD_SOLVER_H
