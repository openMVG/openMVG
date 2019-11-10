// Copyright (C) 2016-2019 Yixuan Qiu <yixuan.qiu@cos.name>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SPARSE_GEN_REAL_SHIFT_SOLVE_H
#define SPARSE_GEN_REAL_SHIFT_SOLVE_H

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <stdexcept>

namespace Spectra {


///
/// \ingroup MatOp
///
/// This class defines the shift-solve operation on a sparse real matrix \f$A\f$,
/// i.e., calculating \f$y=(A-\sigma I)^{-1}x\f$ for any real \f$\sigma\f$ and
/// vector \f$x\f$. It is mainly used in the GenEigsRealShiftSolver eigen solver.
///
template <typename Scalar,int Flags = 0, typename StorageIndex = int>
class SparseGenRealShiftSolve
{
private:
    typedef Eigen::Index Index;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;
    typedef Eigen::Map<const Vector> MapConstVec;
    typedef Eigen::Map<Vector> MapVec;
    typedef Eigen::SparseMatrix<Scalar, Flags, StorageIndex> SparseMatrix;
    typedef const Eigen::Ref<const SparseMatrix> ConstGenericSparseMatrix;

    ConstGenericSparseMatrix m_mat;
    const int m_n;
    Eigen::SparseLU<SparseMatrix> m_solver;

public:
    ///
    /// Constructor to create the matrix operation object.
    ///
    /// \param mat An **Eigen** sparse matrix object, whose type can be
    /// `Eigen::SparseMatrix<Scalar, ...>` or its mapped version
    /// `Eigen::Map<Eigen::SparseMatrix<Scalar, ...> >`.
    ///
    SparseGenRealShiftSolve(ConstGenericSparseMatrix& mat) :
        m_mat(mat), m_n(mat.rows())
    {
        if(mat.rows() != mat.cols())
            throw std::invalid_argument("SparseGenRealShiftSolve: matrix must be square");
    }

    ///
    /// Return the number of rows of the underlying matrix.
    ///
    Index rows() const { return m_n; }
    ///
    /// Return the number of columns of the underlying matrix.
    ///
    Index cols() const { return m_n; }

    ///
    /// Set the real shift \f$\sigma\f$.
    ///
    void set_shift(Scalar sigma)
    {
        SparseMatrix I(m_n, m_n);
        I.setIdentity();

        m_solver.compute(m_mat - sigma * I);
        if(m_solver.info() != Eigen::Success)
            throw std::invalid_argument("SparseGenRealShiftSolve: factorization failed with the given shift");
    }

    ///
    /// Perform the shift-solve operation \f$y=(A-\sigma I)^{-1}x\f$.
    ///
    /// \param x_in  Pointer to the \f$x\f$ vector.
    /// \param y_out Pointer to the \f$y\f$ vector.
    ///
    // y_out = inv(A - sigma * I) * x_in
    void perform_op(const Scalar* x_in, Scalar* y_out) const
    {
        MapConstVec x(x_in,  m_n);
        MapVec      y(y_out, m_n);
        y.noalias() = m_solver.solve(x);
    }
};


} // namespace Spectra

#endif // SPARSE_GEN_REAL_SHIFT_SOLVE_H
