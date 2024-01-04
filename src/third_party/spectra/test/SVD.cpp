#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/SVD>
#include <iostream>
#include <random>  // Requires C++ 11

#include <Spectra/contrib/PartialSVDSolver.h>

using namespace Spectra;

#include "catch.hpp"

using Matrix = Eigen::MatrixXd;
using Vector = Eigen::VectorXd;
using SpMatrix = Eigen::SparseMatrix<double>;

// Generate random sparse matrix
SpMatrix gen_sparse_data(int m, int n, double prob = 0.5)
{
    SpMatrix mat(m, n);
    std::default_random_engine gen;
    gen.seed(0);
    std::uniform_real_distribution<double> distr(0.0, 1.0);
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (distr(gen) < prob)
                mat.insert(i, j) = distr(gen) - 0.5;
        }
    }
    return mat;
}

template <typename MatType>
void run_test(const MatType& mat, int k, int m)
{
    PartialSVDSolver<MatType> svds(mat, k, m);
    int nconv = svds.compute();

    INFO("nconv = " << nconv);
    REQUIRE(nconv == k);

    Vector svals = svds.singular_values();
    Matrix U = svds.matrix_U(k);
    Matrix V = svds.matrix_V(k);

    // SVD solver from Eigen
    // Requires dense matrices
    Matrix mat_dense = Matrix(mat);
    Eigen::JacobiSVD<Matrix> svd(mat, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Vector svals_eigen = svd.singularValues();
    Matrix U_eigen = svd.matrixU();
    Matrix V_eigen = svd.matrixV();

    double err = (svals - svals_eigen.head(k)).array().abs().maxCoeff();
    INFO("Residual of singular values = " << err);
    REQUIRE(err == Approx(0.0).margin(1e-9));

    err = (U.array().abs() - U_eigen.leftCols(k).array().abs()).abs().maxCoeff();
    INFO("Residual of left singular vectors = " << err);
    REQUIRE(err == Approx(0.0).margin(1e-9));

    err = (V.array().abs() - V_eigen.leftCols(k).array().abs()).abs().maxCoeff();
    INFO("Residual of right singular vectors = " << err);
    REQUIRE(err == Approx(0.0).margin(1e-9));
}

TEST_CASE("Partial SVD of tall dense matrix [1000x100]", "[svds_dense_tall]")
{
    std::srand(123);

    const Matrix A = Matrix::Random(1000, 100);
    int k = 5;
    int m = 10;

    run_test<Matrix>(A, k, m);
}

TEST_CASE("Partial SVD of wide dense matrix [1000x100]", "[svds_dense_wide]")
{
    std::srand(123);

    const Matrix A = Matrix::Random(100, 1000);
    int k = 5;
    int m = 10;

    run_test<Matrix>(A, k, m);
}

TEST_CASE("Partial SVD of tall sparse matrix [1000x100]", "[svds_sparse_tall]")
{
    std::srand(123);

    const SpMatrix A = gen_sparse_data(1000, 100, 0.1);
    int k = 5;
    int m = 10;

    run_test<SpMatrix>(A, k, m);
}

TEST_CASE("Partial SVD of wide sparse matrix [1000x100]", "[svds_sparse_wide]")
{
    std::srand(123);

    const SpMatrix A = gen_sparse_data(100, 1000, 0.1);
    int k = 5;
    int m = 10;

    run_test<SpMatrix>(A, k, m);
}
