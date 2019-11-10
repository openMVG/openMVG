// Copyright (C) 2016-2019 Yixuan Qiu <yixuan.qiu@cos.name>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef UPPER_HESSENBERG_QR_H
#define UPPER_HESSENBERG_QR_H

#include <Eigen/Core>
#include <cmath>      // std::sqrt
#include <algorithm>  // std::fill, std::copy
#include <stdexcept>  // std::logic_error

namespace Spectra {


///
/// \defgroup Internals Internal Classes
///
/// Classes for internal use. May be useful to developers.
///

///
/// \ingroup Internals
/// @{
///

///
/// \defgroup LinearAlgebra Linear Algebra
///
/// A number of classes for linear algebra operations.
///

///
/// \ingroup LinearAlgebra
///
/// Perform the QR decomposition of an upper Hessenberg matrix.
///
/// \tparam Scalar The element type of the matrix.
/// Currently supported types are `float`, `double` and `long double`.
///
template <typename Scalar = double>
class UpperHessenbergQR
{
private:
    typedef Eigen::Index Index;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Matrix;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;
    typedef Eigen::Matrix<Scalar, 1, Eigen::Dynamic> RowVector;
    typedef Eigen::Array<Scalar, Eigen::Dynamic, 1> Array;

    typedef Eigen::Ref<Matrix> GenericMatrix;
    typedef const Eigen::Ref<const Matrix> ConstGenericMatrix;

    Matrix m_mat_T;

protected:
    Index m_n;
    // Gi = [ cos[i]  sin[i]]
    //      [-sin[i]  cos[i]]
    // Q = G1 * G2 * ... * G_{n-1}
    Scalar m_shift;
    Array m_rot_cos;
    Array m_rot_sin;
    bool m_computed;

    // Given x and y, compute 1) r = sqrt(x^2 + y^2), 2) c = x / r, 3) s = -y / r
    // If both x and y are zero, set c = 1 and s = 0
    // We must implement it in a numerically stable way
    static void compute_rotation(const Scalar& x, const Scalar& y, Scalar& r, Scalar& c, Scalar& s)
    {
        using std::sqrt;

        const Scalar xsign = (x > Scalar(0)) - (x < Scalar(0));
        const Scalar ysign = (y > Scalar(0)) - (y < Scalar(0));
        const Scalar xabs = x * xsign;
        const Scalar yabs = y * ysign;
        if(xabs > yabs)
        {
            // In this case xabs != 0
            const Scalar ratio = yabs / xabs;  // so that 0 <= ratio < 1
            const Scalar common = sqrt(Scalar(1) + ratio * ratio);
            c = xsign / common;
            r = xabs * common;
            s = -y / r;
        } else {
            if(yabs == Scalar(0))
            {
                r = Scalar(0); c = Scalar(1); s = Scalar(0);
                return;
            }
            const Scalar ratio = xabs / yabs;  // so that 0 <= ratio <= 1
            const Scalar common = sqrt(Scalar(1) + ratio * ratio);
            s = -ysign / common;
            r = yabs * common;
            c = x / r;
        }
    }

public:
    ///
    /// Constructor to preallocate memory. Computation can
    /// be performed later by calling the compute() method.
    ///
    UpperHessenbergQR(Index size) :
        m_n(size),
        m_rot_cos(m_n - 1),
        m_rot_sin(m_n - 1),
        m_computed(false)
    {}

    ///
    /// Constructor to create an object that performs and stores the
    /// QR decomposition of an upper Hessenberg matrix `mat`, with an
    /// optional shift: \f$H-sI=QR\f$. Here \f$H\f$ stands for the matrix
    /// `mat`, and \f$s\f$ is the shift.
    ///
    /// \param mat Matrix type can be `Eigen::Matrix<Scalar, ...>` (e.g.
    /// `Eigen::MatrixXd` and `Eigen::MatrixXf`), or its mapped version
    /// (e.g. `Eigen::Map<Eigen::MatrixXd>`).
    /// Only the upper triangular and the lower subdiagonal parts of
    /// the matrix are used.
    ///
    UpperHessenbergQR(ConstGenericMatrix& mat, const Scalar& shift = Scalar(0)) :
        m_n(mat.rows()),
        m_shift(shift),
        m_rot_cos(m_n - 1),
        m_rot_sin(m_n - 1),
        m_computed(false)
    {
        compute(mat, shift);
    }

    ///
    /// Virtual destructor.
    ///
    virtual ~UpperHessenbergQR() {};

    ///
    /// Conduct the QR factorization of an upper Hessenberg matrix with
    /// an optional shift.
    ///
    /// \param mat Matrix type can be `Eigen::Matrix<Scalar, ...>` (e.g.
    /// `Eigen::MatrixXd` and `Eigen::MatrixXf`), or its mapped version
    /// (e.g. `Eigen::Map<Eigen::MatrixXd>`).
    /// Only the upper triangular and the lower subdiagonal parts of
    /// the matrix are used.
    ///
    virtual void compute(ConstGenericMatrix& mat, const Scalar& shift = Scalar(0))
    {
        m_n = mat.rows();
        if(m_n != mat.cols())
            throw std::invalid_argument("UpperHessenbergQR: matrix must be square");

        m_shift = shift;
        m_mat_T.resize(m_n, m_n);
        m_rot_cos.resize(m_n - 1);
        m_rot_sin.resize(m_n - 1);

        // Make a copy of mat - s * I
        std::copy(mat.data(), mat.data() + mat.size(), m_mat_T.data());
        m_mat_T.diagonal().array() -= m_shift;

        Scalar xi, xj, r, c, s;
        Scalar *Tii, *ptr;
        const Index n1 = m_n - 1;
        for(Index i = 0; i < n1; i++)
        {
            Tii = &m_mat_T.coeffRef(i, i);

            // Make sure mat_T is upper Hessenberg
            // Zero the elements below mat_T(i + 1, i)
            std::fill(Tii + 2, Tii + m_n - i, Scalar(0));

            xi = Tii[0];  // mat_T(i, i)
            xj = Tii[1];  // mat_T(i + 1, i)
            compute_rotation(xi, xj, r, c, s);
            m_rot_cos[i] = c;
            m_rot_sin[i] = s;

            // For a complete QR decomposition,
            // we first obtain the rotation matrix
            // G = [ cos  sin]
            //     [-sin  cos]
            // and then do T[i:(i + 1), i:(n - 1)] = G' * T[i:(i + 1), i:(n - 1)]

            // Gt << c, -s, s, c;
            // m_mat_T.block(i, i, 2, m_n - i) = Gt * m_mat_T.block(i, i, 2, m_n - i);
            Tii[0] = r;  // m_mat_T(i, i)     => r
            Tii[1] = 0;  // m_mat_T(i + 1, i) => 0
            ptr = Tii + m_n; // m_mat_T(i, k), k = i+1, i+2, ..., n-1
            for(Index j = i + 1; j < m_n; j++, ptr += m_n)
            {
                Scalar tmp = ptr[0];
                ptr[0] = c * tmp - s * ptr[1];
                ptr[1] = s * tmp + c * ptr[1];
            }

            // If we do not need to calculate the R matrix, then
            // only the cos and sin sequences are required.
            // In such case we only update T[i + 1, (i + 1):(n - 1)]
            // m_mat_T.block(i + 1, i + 1, 1, m_n - i - 1) *= c;
            // m_mat_T.block(i + 1, i + 1, 1, m_n - i - 1) += s * mat_T.block(i, i + 1, 1, m_n - i - 1);
        }

        m_computed = true;
    }

    ///
    /// Return the \f$R\f$ matrix in the QR decomposition, which is an
    /// upper triangular matrix.
    ///
    /// \return Returned matrix type will be `Eigen::Matrix<Scalar, ...>`, depending on
    /// the template parameter `Scalar` defined.
    ///
    virtual Matrix matrix_R() const
    {
        if(!m_computed)
            throw std::logic_error("UpperHessenbergQR: need to call compute() first");

        return m_mat_T;
    }

    ///
    /// Overwrite `dest` with \f$Q'HQ = RQ + sI\f$, where \f$H\f$ is the input matrix `mat`,
    /// and \f$s\f$ is the shift. The result is an upper Hessenberg matrix.
    ///
    /// \param mat The matrix to be overwritten, whose type should be `Eigen::Matrix<Scalar, ...>`,
    /// depending on the template parameter `Scalar` defined.
    ///
    virtual void matrix_QtHQ(Matrix& dest) const
    {
        if(!m_computed)
            throw std::logic_error("UpperHessenbergQR: need to call compute() first");

        // Make a copy of the R matrix
        dest.resize(m_n, m_n);
        std::copy(m_mat_T.data(), m_mat_T.data() + m_mat_T.size(), dest.data());

        // Compute the RQ matrix
        const Index n1 = m_n - 1;
        for(Index i = 0; i < n1; i++)
        {
            const Scalar c = m_rot_cos.coeff(i);
            const Scalar s = m_rot_sin.coeff(i);
            // RQ[, i:(i + 1)] = RQ[, i:(i + 1)] * Gi
            // Gi = [ cos[i]  sin[i]]
            //      [-sin[i]  cos[i]]
            Scalar *Yi, *Yi1;
            Yi = &dest.coeffRef(0, i);
            Yi1 = Yi + m_n;  // RQ(0, i + 1)
            const Index i2 = i + 2;
            for(Index j = 0; j < i2; j++)
            {
                const Scalar tmp = Yi[j];
                Yi[j]  = c * tmp - s * Yi1[j];
                Yi1[j] = s * tmp + c * Yi1[j];
            }

            /* Vector dest = RQ.block(0, i, i + 2, 1);
            dest.block(0, i, i + 2, 1)     = c * Yi - s * dest.block(0, i + 1, i + 2, 1);
            dest.block(0, i + 1, i + 2, 1) = s * Yi + c * dest.block(0, i + 1, i + 2, 1); */
        }

        // Add the shift to the diagonal
        dest.diagonal().array() += m_shift;
    }

    ///
    /// Apply the \f$Q\f$ matrix to a vector \f$y\f$.
    ///
    /// \param Y A vector that will be overwritten by the matrix product \f$Qy\f$.
    ///
    /// Vector type can be `Eigen::Vector<Scalar, ...>`, depending on
    /// the template parameter `Scalar` defined.
    ///
    // Y -> QY = G1 * G2 * ... * Y
    void apply_QY(Vector& Y) const
    {
        if(!m_computed)
            throw std::logic_error("UpperHessenbergQR: need to call compute() first");

        for(Index i = m_n - 2; i >= 0; i--)
        {
            const Scalar c = m_rot_cos.coeff(i);
            const Scalar s = m_rot_sin.coeff(i);
            // Y[i:(i + 1)] = Gi * Y[i:(i + 1)]
            // Gi = [ cos[i]  sin[i]]
            //      [-sin[i]  cos[i]]
            const Scalar tmp = Y[i];
            Y[i]     =  c * tmp + s * Y[i + 1];
            Y[i + 1] = -s * tmp + c * Y[i + 1];
        }
    }

    ///
    /// Apply the \f$Q\f$ matrix to a vector \f$y\f$.
    ///
    /// \param Y A vector that will be overwritten by the matrix product \f$Q'y\f$.
    ///
    /// Vector type can be `Eigen::Vector<Scalar, ...>`, depending on
    /// the template parameter `Scalar` defined.
    ///
    // Y -> Q'Y = G_{n-1}' * ... * G2' * G1' * Y
    void apply_QtY(Vector& Y) const
    {
        if(!m_computed)
            throw std::logic_error("UpperHessenbergQR: need to call compute() first");

        const Index n1 = m_n - 1;
        for(Index i = 0; i < n1; i++)
        {
            const Scalar c = m_rot_cos.coeff(i);
            const Scalar s = m_rot_sin.coeff(i);
            // Y[i:(i + 1)] = Gi' * Y[i:(i + 1)]
            // Gi = [ cos[i]  sin[i]]
            //      [-sin[i]  cos[i]]
            const Scalar tmp = Y[i];
            Y[i]     = c * tmp - s * Y[i + 1];
            Y[i + 1] = s * tmp + c * Y[i + 1];
        }
    }

    ///
    /// Apply the \f$Q\f$ matrix to another matrix \f$Y\f$.
    ///
    /// \param Y A matrix that will be overwritten by the matrix product \f$QY\f$.
    ///
    /// Matrix type can be `Eigen::Matrix<Scalar, ...>` (e.g.
    /// `Eigen::MatrixXd` and `Eigen::MatrixXf`), or its mapped version
    /// (e.g. `Eigen::Map<Eigen::MatrixXd>`).
    ///
    // Y -> QY = G1 * G2 * ... * Y
    void apply_QY(GenericMatrix Y) const
    {
        if(!m_computed)
            throw std::logic_error("UpperHessenbergQR: need to call compute() first");

        RowVector Yi(Y.cols()), Yi1(Y.cols());
        for(Index i = m_n - 2; i >= 0; i--)
        {
            const Scalar c = m_rot_cos.coeff(i);
            const Scalar s = m_rot_sin.coeff(i);
            // Y[i:(i + 1), ] = Gi * Y[i:(i + 1), ]
            // Gi = [ cos[i]  sin[i]]
            //      [-sin[i]  cos[i]]
            Yi.noalias()  = Y.row(i);
            Yi1.noalias() = Y.row(i + 1);
            Y.row(i)      =  c * Yi + s * Yi1;
            Y.row(i + 1)  = -s * Yi + c * Yi1;
        }
    }

    ///
    /// Apply the \f$Q\f$ matrix to another matrix \f$Y\f$.
    ///
    /// \param Y A matrix that will be overwritten by the matrix product \f$Q'Y\f$.
    ///
    /// Matrix type can be `Eigen::Matrix<Scalar, ...>` (e.g.
    /// `Eigen::MatrixXd` and `Eigen::MatrixXf`), or its mapped version
    /// (e.g. `Eigen::Map<Eigen::MatrixXd>`).
    ///
    // Y -> Q'Y = G_{n-1}' * ... * G2' * G1' * Y
    void apply_QtY(GenericMatrix Y) const
    {
        if(!m_computed)
            throw std::logic_error("UpperHessenbergQR: need to call compute() first");

        RowVector Yi(Y.cols()), Yi1(Y.cols());
        const Index n1 = m_n - 1;
        for(Index i = 0; i < n1; i++)
        {
            const Scalar c = m_rot_cos.coeff(i);
            const Scalar s = m_rot_sin.coeff(i);
            // Y[i:(i + 1), ] = Gi' * Y[i:(i + 1), ]
            // Gi = [ cos[i]  sin[i]]
            //      [-sin[i]  cos[i]]
            Yi.noalias()  = Y.row(i);
            Yi1.noalias() = Y.row(i + 1);
            Y.row(i)      = c * Yi - s * Yi1;
            Y.row(i + 1)  = s * Yi + c * Yi1;
        }
    }

    ///
    /// Apply the \f$Q\f$ matrix to another matrix \f$Y\f$.
    ///
    /// \param Y A matrix that will be overwritten by the matrix product \f$YQ\f$.
    ///
    /// Matrix type can be `Eigen::Matrix<Scalar, ...>` (e.g.
    /// `Eigen::MatrixXd` and `Eigen::MatrixXf`), or its mapped version
    /// (e.g. `Eigen::Map<Eigen::MatrixXd>`).
    ///
    // Y -> YQ = Y * G1 * G2 * ...
    void apply_YQ(GenericMatrix Y) const
    {
        if(!m_computed)
            throw std::logic_error("UpperHessenbergQR: need to call compute() first");

        /*Vector Yi(Y.rows());
        for(Index i = 0; i < m_n - 1; i++)
        {
            const Scalar c = m_rot_cos.coeff(i);
            const Scalar s = m_rot_sin.coeff(i);
            // Y[, i:(i + 1)] = Y[, i:(i + 1)] * Gi
            // Gi = [ cos[i]  sin[i]]
            //      [-sin[i]  cos[i]]
            Yi.noalias() = Y.col(i);
            Y.col(i)     = c * Yi - s * Y.col(i + 1);
            Y.col(i + 1) = s * Yi + c * Y.col(i + 1);
        }*/
        Scalar *Y_col_i, *Y_col_i1;
        const Index n1 = m_n - 1;
        const Index nrow = Y.rows();
        for(Index i = 0; i < n1; i++)
        {
            const Scalar c = m_rot_cos.coeff(i);
            const Scalar s = m_rot_sin.coeff(i);

            Y_col_i  = &Y.coeffRef(0, i);
            Y_col_i1 = &Y.coeffRef(0, i + 1);
            for(Index j = 0; j < nrow; j++)
            {
                Scalar tmp = Y_col_i[j];
                Y_col_i[j]  = c * tmp - s * Y_col_i1[j];
                Y_col_i1[j] = s * tmp + c * Y_col_i1[j];
            }
        }
    }

    ///
    /// Apply the \f$Q\f$ matrix to another matrix \f$Y\f$.
    ///
    /// \param Y A matrix that will be overwritten by the matrix product \f$YQ'\f$.
    ///
    /// Matrix type can be `Eigen::Matrix<Scalar, ...>` (e.g.
    /// `Eigen::MatrixXd` and `Eigen::MatrixXf`), or its mapped version
    /// (e.g. `Eigen::Map<Eigen::MatrixXd>`).
    ///
    // Y -> YQ' = Y * G_{n-1}' * ... * G2' * G1'
    void apply_YQt(GenericMatrix Y) const
    {
        if(!m_computed)
            throw std::logic_error("UpperHessenbergQR: need to call compute() first");

        Vector Yi(Y.rows());
        for(Index i = m_n - 2; i >= 0; i--)
        {
            const Scalar c = m_rot_cos.coeff(i);
            const Scalar s = m_rot_sin.coeff(i);
            // Y[, i:(i + 1)] = Y[, i:(i + 1)] * Gi'
            // Gi = [ cos[i]  sin[i]]
            //      [-sin[i]  cos[i]]
            Yi.noalias() = Y.col(i);
            Y.col(i)     =  c * Yi + s * Y.col(i + 1);
            Y.col(i + 1) = -s * Yi + c * Y.col(i + 1);
        }
    }
};



///
/// \ingroup LinearAlgebra
///
/// Perform the QR decomposition of a tridiagonal matrix, a special
/// case of upper Hessenberg matrices.
///
/// \tparam Scalar The element type of the matrix.
/// Currently supported types are `float`, `double` and `long double`.
///
template <typename Scalar = double>
class TridiagQR: public UpperHessenbergQR<Scalar>
{
private:
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Matrix;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;
    typedef const Eigen::Ref<const Matrix> ConstGenericMatrix;

    typedef typename Matrix::Index Index;

    Vector m_T_diag;   // diagonal elements of T
    Vector m_T_lsub;   // lower subdiagonal of T
    Vector m_T_usub;   // upper subdiagonal of T
    Vector m_T_usub2;  // 2nd upper subdiagonal of T

public:
    ///
    /// Constructor to preallocate memory. Computation can
    /// be performed later by calling the compute() method.
    ///
    TridiagQR(Index size) :
        UpperHessenbergQR<Scalar>(size)
    {}

    ///
    /// Constructor to create an object that performs and stores the
    /// QR decomposition of an upper Hessenberg matrix `mat`, with an
    /// optional shift: \f$H-sI=QR\f$. Here \f$H\f$ stands for the matrix
    /// `mat`, and \f$s\f$ is the shift.
    ///
    /// \param mat Matrix type can be `Eigen::Matrix<Scalar, ...>` (e.g.
    /// `Eigen::MatrixXd` and `Eigen::MatrixXf`), or its mapped version
    /// (e.g. `Eigen::Map<Eigen::MatrixXd>`).
    /// Only the major- and sub- diagonal parts of
    /// the matrix are used.
    ///
    TridiagQR(ConstGenericMatrix& mat, const Scalar& shift = Scalar(0)) :
        UpperHessenbergQR<Scalar>(mat.rows())
    {
        this->compute(mat, shift);
    }

    ///
    /// Conduct the QR factorization of a tridiagonal matrix with an
    /// optional shift.
    ///
    /// \param mat Matrix type can be `Eigen::Matrix<Scalar, ...>` (e.g.
    /// `Eigen::MatrixXd` and `Eigen::MatrixXf`), or its mapped version
    /// (e.g. `Eigen::Map<Eigen::MatrixXd>`).
    /// Only the major- and sub- diagonal parts of
    /// the matrix are used.
    ///
    void compute(ConstGenericMatrix& mat, const Scalar& shift = Scalar(0))
    {
        this->m_n = mat.rows();
        if(this->m_n != mat.cols())
            throw std::invalid_argument("TridiagQR: matrix must be square");

        this->m_shift = shift;
        m_T_diag.resize(this->m_n);
        m_T_lsub.resize(this->m_n - 1);
        m_T_usub.resize(this->m_n - 1);
        m_T_usub2.resize(this->m_n - 2);
        this->m_rot_cos.resize(this->m_n - 1);
        this->m_rot_sin.resize(this->m_n - 1);

        m_T_diag.array() = mat.diagonal().array() - this->m_shift;
        m_T_lsub.noalias() = mat.diagonal(-1);
        m_T_usub.noalias() = m_T_lsub;

        // A number of pointers to avoid repeated address calculation
        Scalar *c = this->m_rot_cos.data(),  // pointer to the cosine vector
               *s = this->m_rot_sin.data(),  // pointer to the sine vector
               r;
        const Index n1 = this->m_n - 1;
        for(Index i = 0; i < n1; i++)
        {
            // diag[i] == T[i, i]
            // lsub[i] == T[i + 1, i]
            // r = sqrt(T[i, i]^2 + T[i + 1, i]^2)
            // c = T[i, i] / r, s = -T[i + 1, i] / r
            this->compute_rotation(m_T_diag.coeff(i), m_T_lsub.coeff(i), r, *c, *s);

            // For a complete QR decomposition,
            // we first obtain the rotation matrix
            // G = [ cos  sin]
            //     [-sin  cos]
            // and then do T[i:(i + 1), i:(i + 2)] = G' * T[i:(i + 1), i:(i + 2)]

            // Update T[i, i] and T[i + 1, i]
            // The updated value of T[i, i] is known to be r
            // The updated value of T[i + 1, i] is known to be 0
            m_T_diag.coeffRef(i) = r;
            m_T_lsub.coeffRef(i) = Scalar(0);
            // Update T[i, i + 1] and T[i + 1, i + 1]
            // usub[i] == T[i, i + 1]
            // diag[i + 1] == T[i + 1, i + 1]
            const Scalar tmp = m_T_usub.coeff(i);
            m_T_usub.coeffRef(i)     = (*c) * tmp - (*s) * m_T_diag.coeff(i + 1);
            m_T_diag.coeffRef(i + 1) = (*s) * tmp + (*c) * m_T_diag.coeff(i + 1);
            // Update T[i, i + 2] and T[i + 1, i + 2]
            // usub2[i] == T[i, i + 2]
            // usub[i + 1] == T[i + 1, i + 2]
            if(i < n1 - 1)
            {
                m_T_usub2.coeffRef(i) = -(*s) * m_T_usub.coeff(i + 1);
                m_T_usub.coeffRef(i + 1) *= (*c);
            }

            c++;
            s++;

            // If we do not need to calculate the R matrix, then
            // only the cos and sin sequences are required.
            // In such case we only update T[i + 1, (i + 1):(i + 2)]
            // T[i + 1, i + 1] = c * T[i + 1, i + 1] + s * T[i, i + 1];
            // T[i + 1, i + 2] *= c;
        }

        this->m_computed = true;
    }

    ///
    /// Return the \f$R\f$ matrix in the QR decomposition, which is an
    /// upper triangular matrix.
    ///
    /// \return Returned matrix type will be `Eigen::Matrix<Scalar, ...>`, depending on
    /// the template parameter `Scalar` defined.
    ///
    Matrix matrix_R() const
    {
        if(!this->m_computed)
            throw std::logic_error("TridiagQR: need to call compute() first");

        Matrix R = Matrix::Zero(this->m_n, this->m_n);
        R.diagonal().noalias() = m_T_diag;
        R.diagonal(1).noalias() = m_T_usub;
        R.diagonal(2).noalias() = m_T_usub2;

        return R;
    }

    ///
    /// Overwrite `dest` with \f$Q'HQ = RQ + sI\f$, where \f$H\f$ is the input matrix `mat`,
    /// and \f$s\f$ is the shift. The result is a tridiagonal matrix.
    ///
    /// \param mat The matrix to be overwritten, whose type should be `Eigen::Matrix<Scalar, ...>`,
    /// depending on the template parameter `Scalar` defined.
    ///
    void matrix_QtHQ(Matrix& dest) const
    {
        if(!this->m_computed)
            throw std::logic_error("TridiagQR: need to call compute() first");

        // Make a copy of the R matrix
        dest.resize(this->m_n, this->m_n);
        dest.setZero();
        dest.diagonal().noalias() = m_T_diag;
        // The upper diagonal refers to m_T_usub
        // The 2nd upper subdiagonal will be zero in RQ

        // Compute the RQ matrix
        // [m11  m12] points to RQ[i:(i+1), i:(i+1)]
        // [0    m22]
        //
        // Gi = [ cos[i]  sin[i]]
        //      [-sin[i]  cos[i]]
        const Index n1 = this->m_n - 1;
        for(Index i = 0; i < n1; i++)
        {
            const Scalar c = this->m_rot_cos.coeff(i);
            const Scalar s = this->m_rot_sin.coeff(i);
            const Scalar m11 = dest.coeff(i, i),
                         m12 = m_T_usub.coeff(i),
                         m22 = m_T_diag.coeff(i + 1);

            // Update the diagonal and the lower subdiagonal of dest
            dest.coeffRef(i    , i    ) = c * m11 - s * m12;
            dest.coeffRef(i + 1, i    ) =         - s * m22;
            dest.coeffRef(i + 1, i + 1) =           c * m22;
        }

        // Copy the lower subdiagonal to upper subdiagonal
        dest.diagonal(1).noalias() = dest.diagonal(-1);

        // Add the shift to the diagonal
        dest.diagonal().array() += this->m_shift;
    }
};

///
/// @}
///


} // namespace Spectra

#endif // UPPER_HESSENBERG_QR_H
