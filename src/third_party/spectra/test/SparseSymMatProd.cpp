
#include <Spectra/MatOp/SparseSymMatProd.h>
#include <Eigen/Core>
#include <Eigen/SparseCore>

using namespace Spectra;

#include "catch.hpp"

using Eigen::Index;
template <typename TestType>
using Matrix = Eigen::Matrix<TestType, Eigen::Dynamic, Eigen::Dynamic>;

template <typename T>
Eigen::SparseMatrix<T> generate_random_sparse(Index rows, Index cols)
{
    std::default_random_engine gen;
    std::uniform_real_distribution<T> dist(0.0, 1.0);

    std::vector<Eigen::Triplet<T>> tripletVector;
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
        {
            auto v_ij = dist(gen);
            if (v_ij < 0.5)
            {
                // if larger than treshold, insert it
                tripletVector.push_back(Eigen::Triplet<T>(i, j, v_ij));
            }
        }
    Eigen::SparseMatrix<T> mat(rows, cols);
    // create the matrix
    mat.setFromTriplets(tripletVector.begin(), tripletVector.end());

    return mat;
}

TEMPLATE_TEST_CASE("matrix operations", "[SparseSymMatProd]", float, double)
{
    std::srand(123);
    constexpr Index n = 100;

    Eigen::SparseMatrix<TestType> mat = generate_random_sparse<TestType>(n, n);
    Eigen::SparseMatrix<TestType> mat1 = mat + Eigen::SparseMatrix<TestType>(mat.transpose());  // It needs to be symetric
    Matrix<TestType> mat2 = Matrix<TestType>::Random(n, n);

    SparseSymMatProd<TestType> sparse1(mat1);
    Matrix<TestType> xs = sparse1 * mat2;
    Matrix<TestType> ys = mat1 * mat2;

    INFO("The matrix-matrix product must be the same as in eigen.")
    REQUIRE(xs.isApprox(ys));
    INFO("The accesor operator must produce the same element as in eigen")
    REQUIRE(mat1.coeff(45, 22) == sparse1(45, 22));
}
