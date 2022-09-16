#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <iostream>
#include <random>  // Requires C++ 11

#include <Spectra/SymGEigsSolver.h>
#include <Spectra/MatOp/SparseSymMatProd.h>
#include <Spectra/MatOp/SparseRegularInverse.h>

using namespace Spectra;

#include "catch.hpp"

using Matrix = Eigen::MatrixXd;
using Vector = Eigen::VectorXd;
using SpMatrix = Eigen::SparseMatrix<double>;

// Generate random sparse matrix
SpMatrix sprand(int size, double prob = 0.5)
{
    SpMatrix mat(size, size);
    std::default_random_engine gen;
    gen.seed(0);
    std::uniform_real_distribution<double> distr(0.0, 1.0);
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            if (distr(gen) < prob)
                mat.insert(i, j) = distr(gen) - 0.5;
        }
    }
    return mat;
}

void gen_sparse_data(int n, SpMatrix& A, SpMatrix& B, double prob = 0.1)
{
    // Eigen solver only uses the lower triangle of A,
    // so we don't need to make A symmetric here.
    A = sprand(n, prob);
    B = A.transpose() * A;
    // To make sure B is positive definite
    for (int i = 0; i < n; i++)
        B.coeffRef(i, i) += 0.1;
}

template <typename Solver>
void run_test(const SpMatrix& A, const SpMatrix& B, Solver& eigs, SortRule selection, bool allow_fail = false)
{
    eigs.init();
    // maxit = 100 to reduce running time for failed cases
    int nconv = eigs.compute(selection, 100);
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

    Vector evals = eigs.eigenvalues();
    Matrix evecs = eigs.eigenvectors();

    Matrix resid = A.template selfadjointView<Eigen::Lower>() * evecs -
        B.template selfadjointView<Eigen::Lower>() * evecs * evals.asDiagonal();
    const double err = resid.array().abs().maxCoeff();

    INFO("||AU - BUD||_inf = " << err);
    REQUIRE(err == Approx(0.0).margin(1e-9));
}

void run_test_sets(const SpMatrix& A, const SpMatrix& B, int k, int m)
{
    using OpType = SparseSymMatProd<double>;
    using BOpType = SparseRegularInverse<double>;
    OpType op(A);
    BOpType Bop(B);
    SymGEigsSolver<OpType, BOpType, GEigsMode::RegularInverse> eigs(op, Bop, k, m);

    SECTION("Largest Magnitude")
    {
        run_test(A, B, eigs, SortRule::LargestMagn);
    }
    SECTION("Largest Value")
    {
        run_test(A, B, eigs, SortRule::LargestAlge);
    }
    SECTION("Smallest Magnitude")
    {
        run_test(A, B, eigs, SortRule::SmallestMagn, true);
    }
    SECTION("Smallest Value")
    {
        run_test(A, B, eigs, SortRule::SmallestAlge);
    }
    SECTION("Both Ends")
    {
        run_test(A, B, eigs, SortRule::BothEnds);
    }
}

TEST_CASE("Generalized eigensolver of sparse symmetric real matrix [10x10]", "[geigs_sym]")
{
    std::srand(123);

    // Eigen solver only uses the lower triangle
    SpMatrix A, B;
    gen_sparse_data(10, A, B, 0.5);
    int k = 3;
    int m = 6;

    run_test_sets(A, B, k, m);
}

TEST_CASE("Generalized eigensolver of sparse symmetric real matrix [100x100]", "[geigs_sym]")
{
    std::srand(123);

    // Eigen solver only uses the lower triangle
    SpMatrix A, B;
    gen_sparse_data(100, A, B, 0.1);
    int k = 10;
    int m = 20;

    run_test_sets(A, B, k, m);
}

TEST_CASE("Generalized eigensolver of sparse symmetric real matrix [1000x1000]", "[geigs_sym]")
{
    std::srand(123);

    // Eigen solver only uses the lower triangle
    SpMatrix A, B;
    gen_sparse_data(1000, A, B, 0.01);
    int k = 20;
    int m = 50;

    run_test_sets(A, B, k, m);
}
