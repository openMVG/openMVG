#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <iostream>
#include <random>  // Requires C++ 11

#include <Spectra/SymGEigsShiftSolver.h>
#include <Spectra/MatOp/SymShiftInvert.h>
#include <Spectra/MatOp/DenseSymMatProd.h>
#include <Spectra/MatOp/SparseSymMatProd.h>

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

// Generate data for testing
void gen_dense_data(int n, Matrix& A, Matrix& B)
{
    // Eigen solver only uses the lower triangle of A,
    // so we don't need to make A symmetric here
    A = Matrix::Random(n, n);
    B = A.transpose() * A;
    // To make sure B is positive definite
    B.diagonal() += Eigen::VectorXd::Random(n).cwiseAbs();
}

void gen_sparse_data(int n, SpMatrix& A, SpMatrix& B, double prob = 0.1)
{
    // Eigen solver only uses the lower triangle of A,
    // so we don't need to make A symmetric here
    A = sprand(n, prob);
    B = A.transpose() * A;
    // To make sure B is positive definite
    for (int i = 0; i < n; i++)
        B.coeffRef(i, i) += 0.1;
}

template <typename TypeA, typename TypeB, typename Solver>
void run_test(const TypeA& A, const TypeB& B, Solver& eigs, SortRule selection, bool allow_fail = false)
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

    INFO("||AU - UD||_inf = " << err);
    REQUIRE(err == Approx(0.0).margin(1e-9));
}

template <typename TypeA, typename TypeB, typename Solver>
void run_test_sets(const TypeA& A, const TypeB& B, Solver& eigs)
{
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

// ======================== Shift-and-invert mode ======================== //
// A is sparse, B is sparse
TEST_CASE("Generalized eigensolver of symmetric real matrices [shift-invert, ss, 100x100]", "[geigs_sym]")
{
    std::srand(123);

    // Eigen solver only uses the lower triangle
    SpMatrix A, B;
    gen_sparse_data(100, A, B, 0.1);
    int k = 10;
    int m = 20;
    double sigma = 1.2345;

    using OpType = SymShiftInvert<double, Eigen::Sparse, Eigen::Sparse>;
    using BOpType = SparseSymMatProd<double>;

    OpType op(A, B);
    BOpType Bop(B);
    SymGEigsShiftSolver<OpType, BOpType, GEigsMode::ShiftInvert> eigs(op, Bop, k, m, sigma);

    run_test_sets(A, B, eigs);
}

// A is sparse, B is dense
TEST_CASE("Generalized eigensolver of symmetric real matrices [shift-invert, sd, 100x100]", "[geigs_sym]")
{
    std::srand(123);

    // Eigen solver only uses the lower triangle
    SpMatrix As, Bs;
    gen_sparse_data(100, As, Bs, 0.1);
    Matrix Ad, Bd;
    gen_dense_data(100, Ad, Bd);
    int k = 10;
    int m = 20;
    double sigma = 1.2345;

    using OpType = SymShiftInvert<double, Eigen::Sparse, Eigen::Dense>;
    using BOpType = DenseSymMatProd<double>;

    OpType op(As, Bd);
    BOpType Bop(Bd);
    SymGEigsShiftSolver<OpType, BOpType, GEigsMode::ShiftInvert> eigs(op, Bop, k, m, sigma);

    run_test_sets(As, Bd, eigs);
}

// A is dense, B is sparse
TEST_CASE("Generalized eigensolver of symmetric real matrices [shift-invert, ds, 100x100]", "[geigs_sym]")
{
    std::srand(123);

    // Eigen solver only uses the lower triangle
    SpMatrix As, Bs;
    gen_sparse_data(100, As, Bs, 0.1);
    Matrix Ad, Bd;
    gen_dense_data(100, Ad, Bd);
    int k = 10;
    int m = 20;
    double sigma = 1.2345;

    using OpType = SymShiftInvert<double, Eigen::Dense, Eigen::Sparse>;
    using BOpType = SparseSymMatProd<double>;

    OpType op(Ad, Bs);
    BOpType Bop(Bs);
    SymGEigsShiftSolver<OpType, BOpType, GEigsMode::ShiftInvert> eigs(op, Bop, k, m, sigma);

    run_test_sets(Ad, Bs, eigs);
}

// A is dense, B is dense
TEST_CASE("Generalized eigensolver of symmetric real matrices [shift-invert, dd, 100x100]", "[geigs_sym]")
{
    std::srand(123);

    // Eigen solver only uses the lower triangle
    Matrix A, B;
    gen_dense_data(100, A, B);
    int k = 10;
    int m = 20;
    double sigma = 1.2345;

    using OpType = SymShiftInvert<double, Eigen::Dense, Eigen::Dense>;
    using BOpType = DenseSymMatProd<double>;

    OpType op(A, B);
    BOpType Bop(B);
    SymGEigsShiftSolver<OpType, BOpType, GEigsMode::ShiftInvert> eigs(op, Bop, k, m, sigma);

    run_test_sets(A, B, eigs);
}

// ======================== Buckling mode ======================== //
// K is sparse, KG is sparse
TEST_CASE("Generalized eigensolver of symmetric real matrices [buckling, ss, 100x100]", "[geigs_sym]")
{
    std::srand(123);

    // Eigen solver only uses the lower triangle
    SpMatrix K, KG;
    gen_sparse_data(100, KG, K, 0.1);
    int k = 10;
    int m = 20;
    double sigma = 1.2345;

    using OpType = SymShiftInvert<double, Eigen::Sparse, Eigen::Sparse>;
    using BOpType = SparseSymMatProd<double>;

    OpType op(K, KG);
    BOpType Bop(K);
    SymGEigsShiftSolver<OpType, BOpType, GEigsMode::Buckling> eigs(op, Bop, k, m, sigma);

    run_test_sets(K, KG, eigs);
}

// K is sparse, KG is dense
TEST_CASE("Generalized eigensolver of symmetric real matrices [buckling, sd, 100x100]", "[geigs_sym]")
{
    std::srand(123);

    // Eigen solver only uses the lower triangle
    SpMatrix Ks, KGs;
    gen_sparse_data(100, KGs, Ks, 0.1);
    Matrix Kd, KGd;
    gen_dense_data(100, KGd, Kd);
    int k = 10;
    int m = 20;
    double sigma = 1.2345;

    using OpType = SymShiftInvert<double, Eigen::Sparse, Eigen::Dense>;
    using BOpType = SparseSymMatProd<double>;

    OpType op(Ks, KGd);
    BOpType Bop(Ks);
    SymGEigsShiftSolver<OpType, BOpType, GEigsMode::Buckling> eigs(op, Bop, k, m, sigma);

    run_test_sets(Ks, KGd, eigs);
}

// K is dense, KG is sparse
TEST_CASE("Generalized eigensolver of symmetric real matrices [buckling, ds, 100x100]", "[geigs_sym]")
{
    std::srand(123);

    // Eigen solver only uses the lower triangle
    SpMatrix Ks, KGs;
    gen_sparse_data(100, KGs, Ks, 0.1);
    Matrix Kd, KGd;
    gen_dense_data(100, KGd, Kd);
    int k = 10;
    int m = 20;
    double sigma = 1.2345;

    using OpType = SymShiftInvert<double, Eigen::Dense, Eigen::Sparse>;
    using BOpType = DenseSymMatProd<double>;

    OpType op(Kd, KGs);
    BOpType Bop(Kd);
    SymGEigsShiftSolver<OpType, BOpType, GEigsMode::Buckling> eigs(op, Bop, k, m, sigma);

    run_test_sets(Kd, KGs, eigs);
}

// K is dense, KG is dense
TEST_CASE("Generalized eigensolver of symmetric real matrices [buckling, dd, 100x100]", "[geigs_sym]")
{
    std::srand(123);

    // Eigen solver only uses the lower triangle
    Matrix K, KG;
    gen_dense_data(100, KG, K);
    int k = 10;
    int m = 20;
    double sigma = 1.2345;

    using OpType = SymShiftInvert<double, Eigen::Dense, Eigen::Dense>;
    using BOpType = DenseSymMatProd<double>;

    OpType op(K, KG);
    BOpType Bop(K);
    SymGEigsShiftSolver<OpType, BOpType, GEigsMode::Buckling> eigs(op, Bop, k, m, sigma);

    run_test_sets(K, KG, eigs);
}

// ======================== Cayley mode ======================== //
// A is sparse, B is sparse
TEST_CASE("Generalized eigensolver of symmetric real matrices [cayley, ss, 100x100]", "[geigs_sym]")
{
    std::srand(123);

    // Eigen solver only uses the lower triangle
    SpMatrix A, B;
    gen_sparse_data(100, A, B, 0.1);
    int k = 10;
    int m = 20;
    double sigma = 1.2345;

    using OpType = SymShiftInvert<double, Eigen::Sparse, Eigen::Sparse>;
    using BOpType = SparseSymMatProd<double>;

    OpType op(A, B);
    BOpType Bop(B);
    SymGEigsShiftSolver<OpType, BOpType, GEigsMode::Cayley> eigs(op, Bop, k, m, sigma);

    run_test_sets(A, B, eigs);
}

// A is sparse, B is dense
TEST_CASE("Generalized eigensolver of symmetric real matrices [cayley, sd, 100x100]", "[geigs_sym]")
{
    std::srand(123);

    // Eigen solver only uses the lower triangle
    SpMatrix As, Bs;
    gen_sparse_data(100, As, Bs, 0.1);
    Matrix Ad, Bd;
    gen_dense_data(100, Ad, Bd);
    int k = 10;
    int m = 20;
    double sigma = 1.2345;

    using OpType = SymShiftInvert<double, Eigen::Sparse, Eigen::Dense>;
    using BOpType = DenseSymMatProd<double>;

    OpType op(As, Bd);
    BOpType Bop(Bd);
    SymGEigsShiftSolver<OpType, BOpType, GEigsMode::Cayley> eigs(op, Bop, k, m, sigma);

    run_test_sets(As, Bd, eigs);
}

// A is dense, B is sparse
TEST_CASE("Generalized eigensolver of symmetric real matrices [cayley, ds, 100x100]", "[geigs_sym]")
{
    std::srand(123);

    // Eigen solver only uses the lower triangle
    SpMatrix As, Bs;
    gen_sparse_data(100, As, Bs, 0.1);
    Matrix Ad, Bd;
    gen_dense_data(100, Ad, Bd);
    int k = 10;
    int m = 20;
    double sigma = 1.2345;

    using OpType = SymShiftInvert<double, Eigen::Dense, Eigen::Sparse>;
    using BOpType = SparseSymMatProd<double>;

    OpType op(Ad, Bs);
    BOpType Bop(Bs);
    SymGEigsShiftSolver<OpType, BOpType, GEigsMode::Cayley> eigs(op, Bop, k, m, sigma);

    run_test_sets(Ad, Bs, eigs);
}

// A is dense, B is dense
TEST_CASE("Generalized eigensolver of symmetric real matrices [cayley, dd, 100x100]", "[geigs_sym]")
{
    std::srand(123);

    // Eigen solver only uses the lower triangle
    Matrix A, B;
    gen_dense_data(100, A, B);
    int k = 10;
    int m = 20;
    double sigma = 1.2345;

    using OpType = SymShiftInvert<double, Eigen::Dense, Eigen::Dense>;
    using BOpType = DenseSymMatProd<double>;

    OpType op(A, B);
    BOpType Bop(B);
    SymGEigsShiftSolver<OpType, BOpType, GEigsMode::Cayley> eigs(op, Bop, k, m, sigma);

    run_test_sets(A, B, eigs);
}
