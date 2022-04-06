#include <Eigen/Core>
#include <Eigen/SparseCore>

#include <Spectra/DavidsonSymEigsSolver.h>
#include <Spectra/MatOp/DenseSymMatProd.h>
#include <Spectra/MatOp/SparseSymMatProd.h>

using namespace Spectra;

#include "catch.hpp"

template <typename T>
using Matrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

template <typename T>
using Vector = Eigen::Matrix<T, Eigen::Dynamic, 1>;

template <typename T>
using SpMatrix = Eigen::SparseMatrix<T>;

// Traits to obtain operation type from matrix type
template <typename MatType>
struct OpTypeTrait
{
    using Scalar = typename MatType::Scalar;
    using OpType = DenseSymMatProd<Scalar>;
};
template <typename T>
struct OpTypeTrait<SpMatrix<T>>
{
    using OpType = SparseSymMatProd<T>;
};

// Generate data for testing
template <typename Matrix>
Matrix gen_sym_data_dense(int n)
{
    Matrix mat = 0.03 * Matrix::Random(n, n);
    Matrix mat1 = mat + mat.transpose();
    for (Eigen::Index i = 0; i < n; i++)
    {
        mat1(i, i) += i + 1;
    }
    return mat1;
}

template <typename SpMatrix>
SpMatrix gen_sym_data_sparse(int n)
{
    // Eigen solver only uses the lower triangle of mat,
    // so we don't need to make mat symmetric here.
    double prob = 0.5;
    SpMatrix mat(n, n);
    std::default_random_engine gen;
    gen.seed(0);
    std::uniform_real_distribution<double> distr(0.0, 1.0);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (distr(gen) < prob)
                mat.insert(i, j) = 0.1 * (distr(gen) - 0.5);
            if (i == j)
            {
                mat.coeffRef(i, j) = i + 1;
            }
        }
    }
    return mat;
}

template <typename MatType>
void run_test(const MatType& mat, int nev, SortRule selection)
{
    using OpType = typename OpTypeTrait<MatType>::OpType;
    OpType op(mat);
    DavidsonSymEigsSolver<OpType> eigs(op, nev);
    int nconv = eigs.compute(selection);

    int niter = eigs.num_iterations();
    REQUIRE(nconv == nev);
    INFO("nconv = " << nconv);
    INFO("niter = " << niter);
    REQUIRE(eigs.info() == CompInfo::Successful);
    using T = typename OpType::Scalar;
    Vector<T> evals = eigs.eigenvalues();
    Matrix<T> evecs = eigs.eigenvectors();

    Matrix<T> resid = op * evecs - evecs * evals.asDiagonal();
    const T err = resid.array().abs().maxCoeff();

    INFO("||AU - UD||_inf = " << err);
    REQUIRE(err < 100 * Eigen::NumTraits<T>::dummy_precision());
}

template <typename MatType>
void run_test_set(const MatType& mat, int k)
{
    SECTION("Largest Value")
    {
        run_test<MatType>(mat, k, SortRule::LargestAlge);
    }

    SECTION("Smallest Value")
    {
        run_test<MatType>(mat, k, SortRule::SmallestAlge);
    }
}

TEMPLATE_TEST_CASE("Davidson Solver of dense symmetric real matrix [1000x1000]", "", double)
{
    std::srand(123);
    const Matrix<TestType> A = gen_sym_data_dense<Matrix<TestType>>(1000);
    int k = 10;
    run_test_set<Matrix<TestType>>(A, k);
}

TEMPLATE_TEST_CASE("Davidson Solver of sparse symmetric real matrix [1000x1000]", "", double)
{
    std::srand(123);
    int k = 10;
    const SpMatrix<TestType> A = gen_sym_data_sparse<SpMatrix<TestType>>(1000);
    run_test_set<SpMatrix<TestType>>(A, k);
}
