// Test ../include/Spectra/LinAlg/UpperHessenbergEigen.h and
//      ../include/Spectra/LinAlg/TridiagEigen.h
#include <ctime>
#include <iostream>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <Spectra/LinAlg/UpperHessenbergEigen.h>
#include <Spectra/LinAlg/TridiagEigen.h>

using namespace Spectra;

#include "catch.hpp"

using Matrix = Eigen::MatrixXd;
using Vector = Eigen::VectorXd;
using ComplexMatrix = Eigen::MatrixXcd;
using ComplexVector = Eigen::VectorXcd;

TEST_CASE("Eigen decomposition of upper Hessenberg matrix", "[Eigen]")
{
    std::srand(123);
    const int n = 100;
    Matrix M = Matrix::Random(n, n);
    M.array() -= 0.5;
    // H is upper Hessenberg
    Matrix H = M.triangularView<Eigen::Upper>();
    H.diagonal(-1) = M.diagonal(-1);

    UpperHessenbergEigen<double> decomp(H);
    ComplexVector evals = decomp.eigenvalues();
    ComplexMatrix evecs = decomp.eigenvectors();

    // Test accuracy
    ComplexMatrix err = H * evecs - evecs * evals.asDiagonal();
    INFO("||HU - UD||_inf = " << err.cwiseAbs().maxCoeff());
    REQUIRE(err.cwiseAbs().maxCoeff() == Approx(0.0).margin(1e-12));

    clock_t t1, t2;
    t1 = clock();
    for (int i = 0; i < 100; i++)
    {
        UpperHessenbergEigen<double> decomp(H);
        ComplexVector evals = decomp.eigenvalues();
        ComplexMatrix evecs = decomp.eigenvectors();
    }
    t2 = clock();
    std::cout << "elapsed time for UpperHessenbergEigen: "
              << double(t2 - t1) / CLOCKS_PER_SEC << " secs\n";

    t1 = clock();
    for (int i = 0; i < 100; i++)
    {
        Eigen::EigenSolver<Matrix> decomp(H);
        ComplexVector evals = decomp.eigenvalues();
        ComplexMatrix evecs = decomp.eigenvectors();
    }
    t2 = clock();
    std::cout << "elapsed time for Eigen::EigenSolver: "
              << double(t2 - t1) / CLOCKS_PER_SEC << " secs\n";
}

TEST_CASE("Eigen decomposition of symmetric tridiagonal matrix", "[Eigen]")
{
    std::srand(123);
    const int n = 100;
    Matrix M = Matrix::Random(n, n);
    M.array() -= 0.5;
    // H is symmetric tridiagonal
    Matrix H = Matrix::Zero(n, n);
    H.diagonal() = M.diagonal();
    H.diagonal(-1) = M.diagonal(-1);
    H.diagonal(1) = M.diagonal(-1);

    TridiagEigen<double> decomp(H);
    Vector evals = decomp.eigenvalues();
    Matrix evecs = decomp.eigenvectors();

    // Test accuracy
    Matrix err = H * evecs - evecs * evals.asDiagonal();
    INFO("||HU - UD||_inf = " << err.cwiseAbs().maxCoeff());
    REQUIRE(err.cwiseAbs().maxCoeff() == Approx(0.0).margin(1e-12));

    clock_t t1, t2;
    t1 = clock();
    for (int i = 0; i < 100; i++)
    {
        TridiagEigen<double> decomp(H);
        Vector evals = decomp.eigenvalues();
        Matrix evecs = decomp.eigenvectors();
    }
    t2 = clock();
    std::cout << "elapsed time for TridiagEigen: "
              << double(t2 - t1) / CLOCKS_PER_SEC << " secs\n";

    t1 = clock();
    for (int i = 0; i < 100; i++)
    {
        Eigen::SelfAdjointEigenSolver<Matrix> decomp(H);
        Vector evals = decomp.eigenvalues();
        Matrix evecs = decomp.eigenvectors();
    }
    t2 = clock();
    std::cout << "elapsed time for Eigen::SelfAdjointEigenSolver: "
              << double(t2 - t1) / CLOCKS_PER_SEC << " secs\n";
}
