// Written by Anna Araslanova
// Modified by Yixuan Qiu
// License: MIT

#ifndef SPECTRA_LOBPCG_SOLVER_H
#define SPECTRA_LOBPCG_SOLVER_H

#include <functional>
#include <map>

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/Eigenvalues>
#include <Eigen/SVD>
#include <Eigen/SparseCholesky>

#include "../SymGEigsSolver.h"

namespace Spectra {

///
/// \ingroup EigenSolver
///

///	*** METHOD
///	The class represent the LOBPCG algorithm, which was invented by Andrew Knyazev
///	Theoretical background of the procedure can be found in the articles below
///	- Knyazev, A.V., 2001. Toward the optimal preconditioned eigensolver : Locally optimal block preconditioned conjugate gradient method.SIAM journal on scientific computing, 23(2), pp.517 - 541.
///	- Knyazev, A.V., Argentati, M.E., Lashuk, I. and Ovtchinnikov, E.E., 2007. Block locally optimal preconditioned eigenvalue xolvers(BLOPEX) in HYPRE and PETSc.SIAM Journal on Scientific Computing, 29(5), pp.2224 - 2239.
///
///	*** CONDITIONS OF USE
///	Locally Optimal Block Preconditioned Conjugate Gradient(LOBPCG) is a method for finding the M smallest eigenvalues
/// and eigenvectors of a large symmetric positive definite generalized eigenvalue problem
///	\f$Ax=\lambda Bx,\f$
///	where \f$A_{NxN}\f$ is a symmetric matrix, \f$B\f$ is symmetric and positive - definite. \f$A and B\f$ are also assumed large and sparse
///	\f$\textit{M}\f$ should be \f$\<< textit{N}\f$ (at least \f$\textit{5M} < \textit{N} \f$)
///
///	*** ARGUMENTS
///	Eigen::SparseMatrix<long double> A; // N*N - Ax = lambda*Bx, lrage and sparse
///	Eigen::SparseMatrix<long double> X; // N*M - initial approximations to eigenvectors (random in general case)
///	Spectra::LOBPCGSolver<long double> solver(A, X);
///	*Eigen::SparseMatrix<long double> B; // N*N - Ax = lambda*Bx, sparse, positive definite
///	solver.setConstraints(B);
///	*Eigen::SparseMatrix<long double> Y; // N*K - constraints, already found eigenvectors
///	solver.setB(B);
///	*Eigen::SparseMatrix<long double> T; // N*N - preconditioner ~ A^-1
///	solver.setPreconditioner(T);
///
///	*** OUTCOMES
///	solver.solve(); // compute eigenpairs // void
///	solver.info(); // state of converjance // int
///	solver.residuals(); // get residuals to evaluate biases // Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>
///	solver.eigenvalues(); // get eigenvalues // Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>
///	solver.eigenvectors(); // get eigenvectors // Eigen::Matrix<Scalar, Eigen::Dynamic, 1>
///
///	*** EXAMPLE
/// \code{.cpp}
/// #include <Spectra/contrib/SymSparseEigsSolverLOBPCG.h>
///
///	// random A
///	Matrix a;
///	a = (Matrix::Random(10, 10).array() > 0.6).cast<long double>() * Matrix::Random(10, 10).array() * 5;
///	a = Matrix((a).triangularView<Eigen::Lower>()) + Matrix((a).triangularView<Eigen::Lower>()).transpose();
///	for (int i = 0; i < 10; i++)
///		a(i, i) = i + 0.5;
///	std::cout << a << "\n";
///	Eigen::SparseMatrix<long double> A(a.sparseView());
///	// random X
///	Eigen::Matrix<long double, 10, 2> x;
///	x = Matrix::Random(10, 2).array();
///	Eigen::SparseMatrix<long double> X(x.sparseView());
///	// solve Ax = lambda*x
///	Spectra::LOBPCGSolver<long double> solver(A, X);
///	solver.compute(10, 1e-4); // 10 iterations, L2_tolerance = 1e-4*N
/// std::cout << "info\n" << solver.info() << std::endl;
/// std::cout << "eigenvalues\n" << solver.eigenvalues() << std::endl;
/// std::cout << "eigenvectors\n" << solver.eigenvectors() << std::endl;
/// std::cout << "residuals\n" << solver.residuals() << std::endl;
/// \endcode
///

template <typename Scalar = long double>
class LOBPCGSolver
{
private:
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Matrix;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;

    typedef std::complex<Scalar> Complex;
    typedef Eigen::Matrix<Complex, Eigen::Dynamic, Eigen::Dynamic> ComplexMatrix;
    typedef Eigen::Matrix<Complex, Eigen::Dynamic, 1> ComplexVector;

    typedef Eigen::SparseMatrix<Scalar> SparseMatrix;
    typedef Eigen::SparseMatrix<Complex> SparseComplexMatrix;

    const int m_n;    // dimension of matrix A
    const int m_nev;  // number of eigenvalues requested
    SparseMatrix A, X;
    SparseMatrix m_Y, m_B, m_preconditioner;
    bool flag_with_constraints, flag_with_B, flag_with_preconditioner;

public:
    SparseMatrix m_residuals;
    Matrix m_evectors;
    Vector m_evalues;
    int m_info;

private:
    // B-orthonormalize matrix M
    int orthogonalizeInPlace(SparseMatrix& M, SparseMatrix& B,
                             SparseMatrix& true_BM, bool has_true_BM = false)
    {
        SparseMatrix BM;

        if (has_true_BM == false)
        {
            if (flag_with_B)
            {
                BM = B * M;
            }
            else
            {
                BM = M;
            }
        }
        else
        {
            BM = true_BM;
        }

        Eigen::SimplicialLDLT<SparseMatrix> chol_MBM(M.transpose() * BM);

        if (chol_MBM.info() != Eigen::Success)
        {
            // LDLT decomposition fail
            m_info = chol_MBM.info();
            return chol_MBM.info();
        }

        SparseComplexMatrix Upper_MBM = chol_MBM.matrixU().template cast<Complex>();
        ComplexVector D_MBM_vec = chol_MBM.vectorD().template cast<Complex>();

        D_MBM_vec = D_MBM_vec.cwiseSqrt();

        for (int i = 0; i < D_MBM_vec.rows(); i++)
        {
            D_MBM_vec(i) = Complex(1.0, 0.0) / D_MBM_vec(i);
        }

        SparseComplexMatrix D_MBM_mat(D_MBM_vec.asDiagonal());

        SparseComplexMatrix U_inv(Upper_MBM.rows(), Upper_MBM.cols());
        U_inv.setIdentity();
        Upper_MBM.template triangularView<Eigen::Upper>().solveInPlace(U_inv);

        SparseComplexMatrix right_product = U_inv * D_MBM_mat;
        M = M * right_product.real();
        if (flag_with_B)
        {
            true_BM = B * M;
        }
        else
        {
            true_BM = M;
        }

        return Eigen::Success;
    }

    void applyConstraintsInPlace(SparseMatrix& X, SparseMatrix& Y,
                                 SparseMatrix& B)
    {
        SparseMatrix BY;
        if (flag_with_B)
        {
            BY = B * Y;
        }
        else
        {
            BY = Y;
        }

        SparseMatrix YBY = Y.transpose() * BY;
        SparseMatrix BYX = BY.transpose() * X;

        SparseMatrix YBY_XYX = (Matrix(YBY).bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(Matrix(BYX))).sparseView();
        X = X - Y * YBY_XYX;
    }

    /*
                return
                'AB
                CD'
                */
    Matrix stack_4_matricies(Matrix A, Matrix B,
                             Matrix C, Matrix D)
    {
        Matrix result(A.rows() + C.rows(), A.cols() + B.cols());
        result.topLeftCorner(A.rows(), A.cols()) = A;
        result.topRightCorner(B.rows(), B.cols()) = B;
        result.bottomLeftCorner(C.rows(), C.cols()) = C;
        result.bottomRightCorner(D.rows(), D.cols()) = D;
        return result;
    }

    Matrix stack_9_matricies(Matrix A, Matrix B, Matrix C,
                             Matrix D, Matrix E, Matrix F,
                             Matrix G, Matrix H, Matrix I)
    {
        Matrix result(A.rows() + D.rows() + G.rows(), A.cols() + B.cols() + C.cols());
        result.block(0, 0, A.rows(), A.cols()) = A;
        result.block(0, A.cols(), B.rows(), B.cols()) = B;
        result.block(0, A.cols() + B.cols(), C.rows(), C.cols()) = C;
        result.block(A.rows(), 0, D.rows(), D.cols()) = D;
        result.block(A.rows(), A.cols(), E.rows(), E.cols()) = E;
        result.block(A.rows(), A.cols() + B.cols(), F.rows(), F.cols()) = F;
        result.block(A.rows() + D.rows(), 0, G.rows(), G.cols()) = G;
        result.block(A.rows() + D.rows(), A.cols(), H.rows(), H.cols()) = H;
        result.block(A.rows() + D.rows(), A.cols() + B.cols(), I.rows(), I.cols()) = I;

        return result;
    }

    void sort_epairs(Vector& evalues, Matrix& evectors, SortRule SelectionRule)
    {
        std::function<bool(Scalar, Scalar)> cmp;
        if (SelectionRule == SortRule::SmallestAlge)
            cmp = std::less<Scalar>{};
        else
            cmp = std::greater<Scalar>{};

        std::map<Scalar, Vector, decltype(cmp)> epairs(cmp);
        for (int i = 0; i < m_evectors.cols(); ++i)
            epairs.insert(std::make_pair(evalues(i), evectors.col(i)));

        int i = 0;
        for (auto& epair : epairs)
        {
            evectors.col(i) = epair.second;
            evalues(i) = epair.first;
            i++;
        }
    }

    void removeColumns(SparseMatrix& matrix, std::vector<int>& colToRemove)
    {
        // remove columns through matrix multiplication
        SparseMatrix new_matrix(matrix.cols(), matrix.cols() - int(colToRemove.size()));
        int iCol = 0;
        std::vector<Eigen::Triplet<Scalar>> tripletList;
        tripletList.reserve(matrix.cols() - int(colToRemove.size()));

        for (int iRow = 0; iRow < matrix.cols(); iRow++)
        {
            if (std::find(colToRemove.begin(), colToRemove.end(), iRow) == colToRemove.end())
            {
                tripletList.push_back(Eigen::Triplet<Scalar>(iRow, iCol, 1));
                iCol++;
            }
        }

        new_matrix.setFromTriplets(tripletList.begin(), tripletList.end());
        matrix = matrix * new_matrix;
    }

    int checkConvergence_getBlocksize(SparseMatrix& m_residuals, Scalar tolerance_L2, std::vector<int>& columnsToDelete)
    {
        // square roots from sum of squares by column
        int BlockSize = m_nev;
        Scalar sum, buffer;

        for (int iCol = 0; iCol < m_nev; iCol++)
        {
            sum = 0;
            for (int iRow = 0; iRow < m_n; iRow++)
            {
                buffer = m_residuals.coeff(iRow, iCol);
                sum += buffer * buffer;
            }

            if (sqrt(sum) < tolerance_L2)
            {
                BlockSize--;
                columnsToDelete.push_back(iCol);
            }
        }
        return BlockSize;
    }

public:
    LOBPCGSolver(const SparseMatrix& A, const SparseMatrix X) :
        m_n(A.rows()),
        m_nev(X.cols()),
        A(A),
        X(X),
        flag_with_constraints(false),
        flag_with_B(false),
        flag_with_preconditioner(false),
        m_info(Eigen::InvalidInput)
    {
        if (A.rows() != X.rows() || A.rows() != A.cols())
            throw std::invalid_argument("Wrong size");

        // if (m_n < 5* m_nev)
        //	throw std::invalid_argument("The problem size is small compared to the block size. Use standard eigensolver");
    }

    void setConstraints(const SparseMatrix& Y)
    {
        m_Y = Y;
        flag_with_constraints = true;
    }

    void setB(const SparseMatrix& B)
    {
        if (B.rows() != A.rows() || B.cols() != A.cols())
            throw std::invalid_argument("Wrong size");
        m_B = B;
        flag_with_B = true;
    }

    void setPreconditioner(const SparseMatrix& preconditioner)
    {
        m_preconditioner = preconditioner;
        flag_with_preconditioner = true;
    }

    void compute(int maxit = 10, Scalar tol_div_n = 1e-7)
    {
        Scalar tolerance_L2 = tol_div_n * m_n;
        int BlockSize;
        int max_iter = std::min(m_n, maxit);

        SparseMatrix directions, AX, AR, BX, AD, ADD, DD, BDD, BD, XAD, RAD, DAD, XBD, RBD, BR, sparse_eVecX, sparse_eVecR, sparse_eVecD, inverse_matrix;
        Matrix XAR, RAR, XBR, gramA, gramB, eVecX, eVecR, eVecD;
        std::vector<int> columnsToDelete;

        if (flag_with_constraints)
        {
            // Apply the constraints Y to X
            applyConstraintsInPlace(X, m_Y, m_B);
        }

        // Make initial vectors orthonormal
        // implicit BX declaration
        if (orthogonalizeInPlace(X, m_B, BX) != Eigen::Success)
        {
            max_iter = 0;
        }

        AX = A * X;
        // Solve the following NxN eigenvalue problem for all N eigenvalues and -vectors:
        // first approximation via a dense problem
        Eigen::EigenSolver<Matrix> eigs(Matrix(X.transpose() * AX));

        if (eigs.info() != Eigen::Success)
        {
            m_info = eigs.info();
            max_iter = 0;
        }
        else
        {
            m_evalues = eigs.eigenvalues().real();
            m_evectors = eigs.eigenvectors().real();
            sort_epairs(m_evalues, m_evectors, SortRule::SmallestAlge);
            sparse_eVecX = m_evectors.sparseView();

            X = X * sparse_eVecX;
            AX = AX * sparse_eVecX;
            BX = BX * sparse_eVecX;
        }

        for (int iter_num = 0; iter_num < max_iter; iter_num++)
        {
            m_residuals.resize(m_n, m_nev);
            for (int i = 0; i < m_nev; i++)
            {
                m_residuals.col(i) = AX.col(i) - m_evalues(i) * BX.col(i);
            }
            BlockSize = checkConvergence_getBlocksize(m_residuals, tolerance_L2, columnsToDelete);

            if (BlockSize == 0)
            {
                m_info = Eigen::Success;
                break;
            }

            // substitution of the original active mask
            if (columnsToDelete.size() > 0)
            {
                removeColumns(m_residuals, columnsToDelete);
                if (iter_num > 0)
                {
                    removeColumns(directions, columnsToDelete);
                    removeColumns(AD, columnsToDelete);
                    removeColumns(BD, columnsToDelete);
                }
                columnsToDelete.clear();  // for next iteration
            }

            if (flag_with_preconditioner)
            {
                // Apply the preconditioner to the residuals
                m_residuals = m_preconditioner * m_residuals;
            }

            if (flag_with_constraints)
            {
                // Apply the constraints Y to residuals
                applyConstraintsInPlace(m_residuals, m_Y, m_B);
            }

            if (orthogonalizeInPlace(m_residuals, m_B, BR) != Eigen::Success)
            {
                break;
            }
            AR = A * m_residuals;

            // Orthonormalize conjugate directions
            if (iter_num > 0)
            {
                if (orthogonalizeInPlace(directions, m_B, BD, true) != Eigen::Success)
                {
                    break;
                }
                AD = A * directions;
            }

            // Perform the Rayleigh Ritz Procedure
            XAR = Matrix(X.transpose() * AR);
            RAR = Matrix(m_residuals.transpose() * AR);
            XBR = Matrix(X.transpose() * BR);

            if (iter_num > 0)
            {
                XAD = X.transpose() * AD;
                RAD = m_residuals.transpose() * AD;
                DAD = directions.transpose() * AD;
                XBD = X.transpose() * BD;
                RBD = m_residuals.transpose() * BD;

                gramA = stack_9_matricies(m_evalues.asDiagonal(), XAR, XAD, XAR.transpose(), RAR, RAD, XAD.transpose(), RAD.transpose(), DAD.transpose());
                gramB = stack_9_matricies(Matrix::Identity(m_nev, m_nev), XBR, XBD, XBR.transpose(), Matrix::Identity(BlockSize, BlockSize), RBD, XBD.transpose(), RBD.transpose(), Matrix::Identity(BlockSize, BlockSize));
            }
            else
            {
                gramA = stack_4_matricies(m_evalues.asDiagonal(), XAR, XAR.transpose(), RAR);
                gramB = stack_4_matricies(Matrix::Identity(m_nev, m_nev), XBR, XBR.transpose(), Matrix::Identity(BlockSize, BlockSize));
            }

            // Calculate the lowest/largest m eigenpairs; Solve the generalized eigenvalue problem.
            DenseSymMatProd<Scalar> Aop(gramA);
            DenseCholesky<Scalar> Bop(gramB);

            SymGEigsSolver<DenseSymMatProd<Scalar>, DenseCholesky<Scalar>, GEigsMode::Cholesky>
                geigs(Aop, Bop, m_nev, (std::min)(10, int(gramA.rows()) - 1));

            geigs.init();
            geigs.compute(SortRule::SmallestAlge);

            // Mat evecs
            if (geigs.info() == CompInfo::Successful)
            {
                m_evalues = geigs.eigenvalues();
                m_evectors = geigs.eigenvectors();
                sort_epairs(m_evalues, m_evectors, SortRule::SmallestAlge);
            }
            else
            {
                // Problem With General EgenVec
                m_info = Eigen::NoConvergence;
                break;
            }

            // Compute Ritz vectors
            if (iter_num > 0)
            {
                eVecX = m_evectors.block(0, 0, m_nev, m_nev);
                eVecR = m_evectors.block(m_nev, 0, BlockSize, m_nev);
                eVecD = m_evectors.block(m_nev + BlockSize, 0, BlockSize, m_nev);

                sparse_eVecX = eVecX.sparseView();
                sparse_eVecR = eVecR.sparseView();
                sparse_eVecD = eVecD.sparseView();

                DD = m_residuals * sparse_eVecR;  // new conjugate directions
                ADD = AR * sparse_eVecR;
                BDD = BR * sparse_eVecR;

                DD = DD + directions * sparse_eVecD;
                ADD = ADD + AD * sparse_eVecD;
                BDD = BDD + BD * sparse_eVecD;
            }
            else
            {
                eVecX = m_evectors.block(0, 0, m_nev, m_nev);
                eVecR = m_evectors.block(m_nev, 0, BlockSize, m_nev);

                sparse_eVecX = eVecX.sparseView();
                sparse_eVecR = eVecR.sparseView();

                DD = m_residuals * sparse_eVecR;
                ADD = AR * sparse_eVecR;
                BDD = BR * sparse_eVecR;
            }

            X = X * sparse_eVecX + DD;
            AX = AX * sparse_eVecX + ADD;
            BX = BX * sparse_eVecX + BDD;

            directions = DD;
            AD = ADD;
            BD = BDD;

        }  // iteration loop

        // calculate last residuals
        m_residuals.resize(m_n, m_nev);
        for (int i = 0; i < m_nev; i++)
        {
            m_residuals.col(i) = AX.col(i) - m_evalues(i) * BX.col(i);
        }
        BlockSize = checkConvergence_getBlocksize(m_residuals, tolerance_L2, columnsToDelete);

        if (BlockSize == 0)
        {
            m_info = Eigen::Success;
        }
    }  // compute

    Vector eigenvalues()
    {
        return m_evalues;
    }

    Matrix eigenvectors()
    {
        return m_evectors;
    }

    Matrix residuals()
    {
        return Matrix(m_residuals);
    }

    int info() { return m_info; }
};

}  // namespace Spectra

#endif  // SPECTRA_LOBPCG_SOLVER_H
