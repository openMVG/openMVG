#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <iostream>
#include <random> // Requires C++ 11

#include <Spectra/GenEigsSolver.h>
#include <Spectra/MatOp/DenseGenMatProd.h>
#include <Spectra/MatOp/SparseGenMatProd.h>

using namespace Spectra;

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

typedef Eigen::MatrixXd Matrix;
typedef Eigen::VectorXd Vector;
typedef Eigen::MatrixXcd ComplexMatrix;
typedef Eigen::VectorXcd ComplexVector;
typedef Eigen::SparseMatrix<double> SpMatrix;

// Traits to obtain operation type from matrix type
template <typename MatType>
struct OpTypeTrait
{
    typedef DenseGenMatProd<double> OpType;
};

template <>
struct OpTypeTrait<SpMatrix>
{
    typedef SparseGenMatProd<double> OpType;
};

// Generate random sparse matrix
SpMatrix gen_sparse_data(int n, double prob = 0.5)
{
    SpMatrix mat(n, n);
    std::default_random_engine gen;
    gen.seed(0);
    std::uniform_real_distribution<double> distr(0.0, 1.0);
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
        {
            if(distr(gen) < prob)
                mat.insert(i, j) = distr(gen) - 0.5;
        }
    }
    return mat;
}



template <typename MatType, int SelectionRule>
void run_test(const MatType& mat, int k, int m, bool allow_fail = false)
{
    typename OpTypeTrait<MatType>::OpType op(mat);
    GenEigsSolver<double, SelectionRule, typename OpTypeTrait<MatType>::OpType>
        eigs(&op, k, m);
    eigs.init();
    int nconv = eigs.compute(500); // maxit = 500 to reduce running time for failed cases
    int niter = eigs.num_iterations();
    int nops  = eigs.num_operations();

    if(allow_fail)
    {
        if( eigs.info() != SUCCESSFUL )
        {
            WARN( "FAILED on this test" );
            std::cout << "nconv = " << nconv << std::endl;
            std::cout << "niter = " << niter << std::endl;
            std::cout << "nops  = " << nops  << std::endl;
            return;
        }
    } else {
        INFO( "nconv = " << nconv );
        INFO( "niter = " << niter );
        INFO( "nops  = " << nops );
        REQUIRE( eigs.info() == SUCCESSFUL );
    }

    ComplexVector evals = eigs.eigenvalues();
    ComplexMatrix evecs = eigs.eigenvectors();

    ComplexMatrix resid = mat * evecs - evecs * evals.asDiagonal();
    const double err = resid.array().abs().maxCoeff();

    INFO( "||AU - UD||_inf = " << err );
    REQUIRE( err == Approx(0.0).margin(1e-9) );
}

template <typename MatType>
void run_test_sets(const MatType& A, int k, int m)
{
    SECTION( "Largest Magnitude" )
    {
        run_test<MatType, LARGEST_MAGN>(A, k, m);
    }
    SECTION( "Largest Real Part" )
    {
        run_test<MatType, LARGEST_REAL>(A, k, m);
    }
    SECTION( "Largest Imaginary Part" )
    {
        run_test<MatType, LARGEST_IMAG>(A, k, m);
    }
    SECTION( "Smallest Magnitude" )
    {
        run_test<MatType, SMALLEST_MAGN>(A, k, m, true);
    }
    SECTION( "Smallest Real Part" )
    {
        run_test<MatType, SMALLEST_REAL>(A, k, m);
    }
    SECTION( "Smallest Imaginary Part" )
    {
        run_test<MatType, SMALLEST_IMAG>(A, k, m, true);
    }
}

TEST_CASE("Eigensolver of general real matrix [10x10]", "[eigs_gen]")
{
    std::srand(123);

    const Matrix A = Eigen::MatrixXd::Random(10, 10);
    int k = 3;
    int m = 6;

    run_test_sets(A, k, m);
}

TEST_CASE("Eigensolver of general real matrix [100x100]", "[eigs_gen]")
{
    std::srand(123);

    const Matrix A = Eigen::MatrixXd::Random(100, 100);
    int k = 10;
    int m = 20;

    run_test_sets(A, k, m);
}

TEST_CASE("Eigensolver of general real matrix [1000x1000]", "[eigs_gen]")
{
    std::srand(123);

    const Matrix A = Eigen::MatrixXd::Random(1000, 1000);
    int k = 20;
    int m = 50;

    run_test_sets(A, k, m);
}

TEST_CASE("Eigensolver of sparse real matrix [10x10]", "[eigs_gen]")
{
    std::srand(123);

    const SpMatrix A = gen_sparse_data(10, 0.5);
    int k = 3;
    int m = 6;

    run_test_sets(A, k, m);
}

TEST_CASE("Eigensolver of sparse real matrix [100x100]", "[eigs_gen]")
{
    std::srand(123);

    const SpMatrix A = gen_sparse_data(100, 0.5);
    int k = 10;
    int m = 20;

    run_test_sets(A, k, m);
}

TEST_CASE("Eigensolver of sparse real matrix [1000x1000]", "[eigs_gen]")
{
    std::srand(123);

    const SpMatrix A = gen_sparse_data(1000, 0.5);
    int k = 20;
    int m = 50;

    run_test_sets(A, k, m);
}
