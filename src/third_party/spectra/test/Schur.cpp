// Test ../include/Spectra/LinAlg/UpperHessenbergSchur.h
#include <iostream>
#include <Eigen/Core>
#include <Spectra/LinAlg/UpperHessenbergSchur.h>

using namespace Spectra;

#include "catch.hpp"

using Matrix = Eigen::MatrixXd;
using Vector = Eigen::VectorXd;

// Schur decomposition on an upper Hessenberg matrix
void run_test(const Matrix& A)
{
    // Spectra implementation
    UpperHessenbergSchur<double> decomp(A);
    const Matrix& matT = decomp.matrix_T();
    const Matrix& matU = decomp.matrix_U();

    // Test whether T is upper Hessenberg
    Matrix lowerT = matT.triangularView<Eigen::StrictlyLower>();
    lowerT.diagonal(-1).setZero();
    INFO("||lower T||_inf = " << lowerT.cwiseAbs().maxCoeff());
    REQUIRE(lowerT.cwiseAbs().maxCoeff() == Approx(0.0).margin(1e-16));

    // Test whether Q is orthonormal
    constexpr double tol = 1e-12;
    const int n = A.rows();
    Matrix I = Matrix::Identity(n, n);
    Matrix resid1 = matU.transpose() * matU - I;
    Matrix resid2 = matU * matU.transpose() - I;
    INFO("||U'U - I||_inf = " << resid1.cwiseAbs().maxCoeff());
    REQUIRE(resid1.cwiseAbs().maxCoeff() == Approx(0.0).margin(tol));
    INFO("||UU' - I||_inf = " << resid2.cwiseAbs().maxCoeff());
    REQUIRE(resid2.cwiseAbs().maxCoeff() == Approx(0.0).margin(tol));

    // Test whether AU=UT
    Matrix resid3 = A * matU - matU * matT;
    INFO("||AU - UT||_inf = " << resid3.cwiseAbs().maxCoeff());
    REQUIRE(resid3.cwiseAbs().maxCoeff() == Approx(0.0).margin(tol));
}

Matrix gen_upper_hessenberg_mat(int n)
{
    std::srand(123);
    Matrix A = Matrix::Random(n, n);
    A.triangularView<Eigen::StrictlyLower>().setZero();
    A.diagonal(-1).array() = Vector::Random(n - 1).array();
    return A;
}

TEST_CASE("Schur decomposition of upper Hessenberg matrix [10x10]", "[Schur]")
{
    Matrix A = gen_upper_hessenberg_mat(10);
    run_test(A);
}

TEST_CASE("Schur decomposition of upper Hessenberg matrix [100x100]", "[Schur]")
{
    Matrix A = gen_upper_hessenberg_mat(100);
    run_test(A);
}

TEST_CASE("Schur decomposition of upper Hessenberg matrix [500x500]", "[Schur]")
{
    Matrix A = gen_upper_hessenberg_mat(500);
    run_test(A);
}
