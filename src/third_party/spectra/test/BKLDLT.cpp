// Test ../include/Spectra/LinAlg/BKLDLT.h
#include <Eigen/Core>
#include <Spectra/LinAlg/BKLDLT.h>

using namespace Spectra;

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

using Eigen::MatrixXd;
using Eigen::VectorXd;

// Solve (A - s * I)x = b
void run_test(const MatrixXd& A, const VectorXd& b, double s)
{
    BKLDLT<double> decompL(A, Eigen::Lower, s);
    REQUIRE(decompL.info() == SUCCESSFUL);

    BKLDLT<double> decompU(A, Eigen::Upper, s);
    REQUIRE( decompU.info() == SUCCESSFUL );

    VectorXd solL = decompL.solve(b);
    VectorXd solU = decompU.solve(b);
    REQUIRE( (solL - solU).cwiseAbs().maxCoeff() == 0.0 );

    const double tol = 1e-9;
    VectorXd resid = A * solL - s * solL - b;
    INFO( "||(A - s * I)x - b||_inf = " << resid.cwiseAbs().maxCoeff() );
    REQUIRE( resid.cwiseAbs().maxCoeff() == Approx(0.0).margin(tol) );
}

TEST_CASE("BKLDLT decomposition of symmetric real matrix [10x10]", "[BKLDLT]")
{
    std::srand(123);
    const int n = 10;
    MatrixXd A = MatrixXd::Random(n, n);
    A = (A + A.transpose()).eval();
    VectorXd b = VectorXd::Random(n);
    const double shift = 1.0;

    run_test(A, b, shift);
}

TEST_CASE("BKLDLT decomposition of symmetric real matrix [100x100]", "[BKLDLT]")
{
    std::srand(123);
    const int n = 100;
    MatrixXd A = MatrixXd::Random(n, n);
    A = (A + A.transpose()).eval();
    VectorXd b = VectorXd::Random(n);
    const double shift = 1.0;

    run_test(A, b, shift);
}

TEST_CASE("BKLDLT decomposition of symmetric real matrix [1000x1000]", "[BKLDLT]")
{
    std::srand(3);
    const int n = 1000;
    MatrixXd A = MatrixXd::Random(n, n);
    A = (A + A.transpose()).eval();
    VectorXd b = VectorXd::Random(n);
    const double shift = 1.0;

    run_test(A, b, shift);
}
