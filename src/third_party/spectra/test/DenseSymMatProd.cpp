
#include <Spectra/MatOp/DenseSymMatProd.h>
#include <Eigen/Dense>
#include <Eigen/Core>

using namespace Spectra;

#include "catch.hpp"
#include <complex>

template <typename TestType>
using Matrix = Eigen::Matrix<TestType, Eigen::Dynamic, Eigen::Dynamic>;
using Eigen::Index;

TEMPLATE_TEST_CASE("matrix operations", "[DenseSymMatProd]", float, double)
{
    std::srand(123);
    constexpr Index n = 100;

    Matrix<TestType> mat = Matrix<TestType>::Random(n, n);
    Matrix<TestType> mat1 = mat + mat.transpose();  // It needs to be symetric
    Matrix<TestType> mat2 = Matrix<TestType>::Random(n, n);

    DenseSymMatProd<TestType> dense1(mat1);
    Matrix<TestType> xs = dense1 * mat2;
    Matrix<TestType> ys = mat1 * mat2;

    INFO("The matrix-matrix product must be the same as in eigen.")
    REQUIRE(xs.isApprox(ys));
    INFO("The accesor operator must produce the same element as in eigen")
    REQUIRE(mat1(15, 23) == dense1(23, 15));
}
