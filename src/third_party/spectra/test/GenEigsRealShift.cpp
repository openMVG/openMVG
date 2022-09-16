#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <iostream>
#include <type_traits>
#include <random>  // Requires C++ 11

#include <Spectra/GenEigsRealShiftSolver.h>
#include <Spectra/MatOp/DenseGenRealShiftSolve.h>
#include <Spectra/MatOp/SparseGenRealShiftSolve.h>

using namespace Spectra;

#include "catch.hpp"

using Matrix = Eigen::MatrixXd;
using Vector = Eigen::VectorXd;
using ComplexMatrix = Eigen::MatrixXcd;
using ComplexVector = Eigen::VectorXcd;
using SpMatrix = Eigen::SparseMatrix<double>;

// Generate random sparse matrix
SpMatrix gen_sparse_data(int n, double prob = 0.5)
{
    SpMatrix mat(n, n);
    std::default_random_engine gen;
    gen.seed(0);
    std::uniform_real_distribution<double> distr(0.0, 1.0);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (distr(gen) < prob)
                mat.insert(i, j) = distr(gen) - 0.5;
        }
    }
    return mat;
}

template <typename MatType, typename Solver>
void run_test(const MatType& mat, Solver& eigs, SortRule selection, bool allow_fail = false)
{
    eigs.init();
    // maxit = 200 to reduce running time for failed cases
    int nconv = eigs.compute(selection, 200);
    int niter = eigs.num_iterations();
    int nops = eigs.num_operations();

    if (allow_fail && eigs.info() != CompInfo::Successful)
    {
        WARN("FAILED on this test");
        std::cout << "nconv = " << nconv << std::endl;
        std::cout << "niter = " << niter << std::endl;
        std::cout << "nops  = " << nops << std::endl;
        return;
    }
    else
    {
        INFO("nconv = " << nconv);
        INFO("niter = " << niter);
        INFO("nops  = " << nops);
        REQUIRE(eigs.info() == CompInfo::Successful);
    }

    ComplexVector evals = eigs.eigenvalues();
    ComplexMatrix evecs = eigs.eigenvectors();

    ComplexMatrix resid = mat * evecs - evecs * evals.asDiagonal();
    const double err = resid.array().abs().maxCoeff();

    INFO("||AU - UD||_inf = " << err);
    REQUIRE(err == Approx(0.0).margin(1e-8));
}

template <typename MatType>
void run_test_sets(const MatType& A, int k, int m, double sigma)
{
    constexpr bool is_dense = std::is_same<MatType, Matrix>::value;
    using DenseOp = DenseGenRealShiftSolve<double>;
    using SparseOp = SparseGenRealShiftSolve<double>;
    using OpType = typename std::conditional<is_dense, DenseOp, SparseOp>::type;

    OpType op(A);
    GenEigsRealShiftSolver<OpType> eigs(op, k, m, sigma);

    SECTION("Largest Magnitude")
    {
        run_test(A, eigs, SortRule::LargestMagn);
    }
    SECTION("Largest Real Part")
    {
        run_test(A, eigs, SortRule::LargestReal);
    }
    SECTION("Largest Imaginary Part")
    {
        run_test(A, eigs, SortRule::LargestImag);
    }
    SECTION("Smallest Magnitude")
    {
        run_test(A, eigs, SortRule::SmallestMagn, true);
    }
    SECTION("Smallest Real Part")
    {
        run_test(A, eigs, SortRule::SmallestReal);
    }
    SECTION("Smallest Imaginary Part")
    {
        run_test(A, eigs, SortRule::SmallestImag, true);
    }
}

TEST_CASE("Eigensolver of general real matrix [10x10]", "[eigs_gen]")
{
    std::srand(123);

    const Matrix A = Matrix::Random(10, 10);
    int k = 3;
    int m = 8;
    double sigma = 1.0;

    run_test_sets(A, k, m, sigma);
}

TEST_CASE("Eigensolver of general real matrix [100x100]", "[eigs_gen]")
{
    std::srand(123);

    const Matrix A = Matrix::Random(100, 100);
    int k = 10;
    int m = 20;
    double sigma = 10.0;

    run_test_sets(A, k, m, sigma);
}

TEST_CASE("Eigensolver of general real matrix [1000x1000]", "[eigs_gen]")
{
    std::srand(123);

    const Matrix A = Matrix::Random(1000, 1000);
    int k = 20;
    int m = 50;
    double sigma = 100.0;

    run_test_sets(A, k, m, sigma);
}

TEST_CASE("Eigensolver of sparse real matrix [10x10]", "[eigs_gen]")
{
    std::srand(123);

    const SpMatrix A = gen_sparse_data(10, 0.5);
    int k = 3;
    int m = 6;
    double sigma = 1.0;

    run_test_sets(A, k, m, sigma);
}

TEST_CASE("Eigensolver of sparse real matrix [100x100]", "[eigs_gen]")
{
    std::srand(123);

    const SpMatrix A = gen_sparse_data(100, 0.1);
    int k = 10;
    int m = 30;
    double sigma = 10.0;

    run_test_sets(A, k, m, sigma);
}

TEST_CASE("Eigensolver of sparse real matrix [1000x1000]", "[eigs_gen]")
{
    std::srand(123);

    const SpMatrix A = gen_sparse_data(1000, 0.01);
    int k = 20;
    int m = 50;
    double sigma = 100.0;

    run_test_sets(A, k, m, sigma);
}
