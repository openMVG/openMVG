#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <iostream>
#include <random> // Requires C++ 11

#include <Spectra/SymEigsSolver.h>
#include <Spectra/MatOp/DenseSymMatProd.h>
#include <Spectra/MatOp/SparseSymMatProd.h>

using namespace Spectra;

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

typedef Eigen::MatrixXd Matrix;
typedef Eigen::VectorXd Vector;
typedef Eigen::SparseMatrix<double> SpMatrix;

// Traits to obtain operation type from matrix type
template <typename MatType>
struct OpTypeTrait
{
    typedef DenseSymMatProd<double> OpType;
};

template <>
struct OpTypeTrait<SpMatrix>
{
    typedef SparseSymMatProd<double> OpType;
};

// Generate data for testing
Matrix gen_dense_data(int n)
{
    const Matrix mat = Eigen::MatrixXd::Random(n, n);
    return mat + mat.transpose();
}

SpMatrix gen_sparse_data(int n, double prob = 0.5)
{
    // Eigen solver only uses the lower triangle of mat,
    // so we don't need to make mat symmetric here.
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
void run_test(const MatType& mat, int k, int m)
{
    typename OpTypeTrait<MatType>::OpType op(mat);
    SymEigsSolver<double, SelectionRule, typename OpTypeTrait<MatType>::OpType>
        eigs(&op, k, m);
    eigs.init();
    int nconv = eigs.compute();
    int niter = eigs.num_iterations();
    int nops  = eigs.num_operations();

    INFO( "nconv = " << nconv );
    INFO( "niter = " << niter );
    INFO( "nops  = " << nops  );
    REQUIRE( eigs.info() == SUCCESSFUL );

    Vector evals = eigs.eigenvalues();
    Matrix evecs = eigs.eigenvectors();

    Matrix resid = mat.template selfadjointView<Eigen::Lower>() * evecs - evecs * evals.asDiagonal();
    const double err = resid.array().abs().maxCoeff();

    INFO( "||AU - UD||_inf = " << err );
    REQUIRE( err == Approx(0.0).margin(1e-9) );
}

template <typename MatType>
void run_test_sets(const MatType& mat, int k, int m)
{
    SECTION( "Largest Magnitude" )
    {
        run_test<MatType, LARGEST_MAGN>(mat, k, m);
    }
    SECTION( "Largest Value" )
    {
        run_test<MatType, LARGEST_ALGE>(mat, k, m);
    }
    SECTION( "Smallest Magnitude" )
    {
        run_test<MatType, SMALLEST_MAGN>(mat, k, m);
    }
    SECTION( "Smallest Value" )
    {
        run_test<MatType, SMALLEST_ALGE>(mat, k, m);
    }
    SECTION( "Both Ends" )
    {
        run_test<MatType, BOTH_ENDS>(mat, k, m);
    }
}

TEST_CASE("Eigensolver of symmetric real matrix [10x10]", "[eigs_sym]")
{
    std::srand(123);

    const Matrix A = gen_dense_data(10);
    int k = 3;
    int m = 6;

    run_test_sets(A, k, m);
}

TEST_CASE("Eigensolver of symmetric real matrix [100x100]", "[eigs_sym]")
{
    std::srand(123);

    const Matrix A = gen_dense_data(100);
    int k = 10;
    int m = 20;

    run_test_sets(A, k, m);
}

TEST_CASE("Eigensolver of symmetric real matrix [1000x1000]", "[eigs_sym]")
{
    std::srand(123);

    const Matrix A = gen_dense_data(1000);
    int k = 20;
    int m = 50;

    run_test_sets(A, k, m);
}

TEST_CASE("Eigensolver of sparse symmetric real matrix [10x10]", "[eigs_sym]")
{
    std::srand(123);

    // Eigen solver only uses the lower triangle
    const SpMatrix A = gen_sparse_data(10, 0.5);
    int k = 3;
    int m = 6;

    run_test_sets(A, k, m);
}

TEST_CASE("Eigensolver of sparse symmetric real matrix [100x100]", "[eigs_sym]")
{
    std::srand(123);

    // Eigen solver only uses the lower triangle
    const SpMatrix A = gen_sparse_data(100, 0.5);
    int k = 10;
    int m = 20;

    run_test_sets(A, k, m);
}

TEST_CASE("Eigensolver of sparse symmetric real matrix [1000x1000]", "[eigs_sym]")
{
    std::srand(123);

    // Eigen solver only uses the lower triangle
    const SpMatrix A = gen_sparse_data(1000, 0.5);
    int k = 20;
    int m = 50;

    run_test_sets(A, k, m);
}
