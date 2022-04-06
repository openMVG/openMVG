
#include <Spectra/MatOp/DenseGenMatProd.h>
#include <Eigen/Dense>

using namespace Spectra;

#include "catch.hpp"

using Eigen::Matrix;
using Eigen::Index;

TEMPLATE_TEST_CASE("matrix operations", "[DenseGenMatProd]", float, double)
{
    std::srand(123);
    constexpr Index n = 100;

    Matrix<TestType, -1, -1> mat1 = Matrix<TestType, -1, -1>::Random(n, n);
    Matrix<TestType, -1, -1> mat2 = Matrix<TestType, -1, -1>::Random(n, n);

    DenseGenMatProd<TestType> dense1(mat1);
    Matrix<TestType, -1, -1> xs = dense1 * mat2;
    Matrix<TestType, -1, -1> ys = mat1 * mat2;

    INFO("The matrix-matrix product must be the same as in eigen.")
    REQUIRE(xs.isApprox(ys));
    INFO("The accesor operator must produce the same element as in eigen")
    REQUIRE(mat1(23, 87) == dense1(23, 87));
}
