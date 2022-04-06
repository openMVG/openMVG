// Test ../include/Spectra/LinAlg/BKLDLT.h
#include <Eigen/Core>
#include <Spectra/LinAlg/BKLDLT.h>

using namespace Spectra;

#include "catch.hpp"

using Matrix = Eigen::MatrixXd;
using Vector = Eigen::VectorXd;

// Solve (A - s * I)x = b
void run_test(const Matrix& A, const Vector& b, double s)
{
    // Test decomposition using only the lower triangular part
    BKLDLT<double> decompL(A, Eigen::Lower, s);
    REQUIRE(decompL.info() == CompInfo::Successful);

    // Test decomposition using only the upper triangular part
    BKLDLT<double> decompU(A, Eigen::Upper, s);
    REQUIRE(decompU.info() == CompInfo::Successful);

    // Test whether the solutions are identical
    Vector solL = decompL.solve(b);
    Vector solU = decompU.solve(b);
    REQUIRE((solL - solU).cwiseAbs().maxCoeff() == 0.0);

    // Test the accuracy of the solution
    constexpr double tol = 1e-9;
    Vector resid = A * solL - s * solL - b;
    INFO("||(A - s * I)x - b||_inf = " << resid.cwiseAbs().maxCoeff());
    REQUIRE(resid.cwiseAbs().maxCoeff() == Approx(0.0).margin(tol));
}

TEST_CASE("BKLDLT decomposition of symmetric real matrix [10x10]", "[BKLDLT]")
{
    std::srand(123);
    const int n = 10;
    Matrix A = Matrix::Random(n, n);
    A = (A + A.transpose()).eval();
    Vector b = Vector::Random(n);
    const double shift = 1.0;

    run_test(A, b, shift);
}

TEST_CASE("BKLDLT decomposition of symmetric real matrix [100x100]", "[BKLDLT]")
{
    std::srand(123);
    const int n = 100;
    Matrix A = Matrix::Random(n, n);
    A = (A + A.transpose()).eval();
    Vector b = Vector::Random(n);
    const double shift = 1.0;

    run_test(A, b, shift);
}

TEST_CASE("BKLDLT decomposition of symmetric real matrix [1000x1000]", "[BKLDLT]")
{
    std::srand(123);
    const int n = 1000;
    Matrix A = Matrix::Random(n, n);
    A = (A + A.transpose()).eval();
    Vector b = Vector::Random(n);
    const double shift = 1.0;

    run_test(A, b, shift);
}
