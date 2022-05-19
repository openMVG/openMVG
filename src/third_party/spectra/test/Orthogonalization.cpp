#include <Eigen/Core>
#include <Spectra/LinAlg/Orthogonalization.h>
#include <iostream>
using namespace Spectra;

#include "catch.hpp"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Index;

template <typename Matrix>
void check_orthogonality(const Matrix& basis)
{
    const double tol = 1e-12;
    Matrix xs = basis.transpose() * basis;
    INFO("The orthonormalized basis must fulfill that basis.T * basis = I");
    INFO("Matrix is\n " << basis);
    INFO("Overlap is\n " << xs);
    CHECK(xs.isIdentity(tol));
}

TEST_CASE("complete orthonormalization", "[orthogonalisation]")
{
    std::srand(123);
    const Index n = 20;

    MatrixXd mat = MatrixXd::Random(n, n);

    SECTION("MGS")
    {
        MGS_orthogonalisation(mat);
        check_orthogonality(mat);
    }

    SECTION("GS")
    {
        GS_orthogonalisation(mat);
        check_orthogonality(mat);
    }

    SECTION("QR")
    {
        QR_orthogonalisation(mat);
        check_orthogonality(mat);
    }

    SECTION("twice_is_enough")
    {
        twice_is_enough_orthogonalisation(mat);
        check_orthogonality(mat);
    }

    SECTION("JensWehner")
    {
        JensWehner_orthogonalisation(mat);
        check_orthogonality(mat);
    }
}

TEST_CASE("Partial orthonormalization", "[orthogonalisation]")
{
    std::srand(123);
    const Index n = 20;
    const Index sub = 5;
    Index start = n - sub;

    // Create a n x 20 orthonormal basis
    MatrixXd mat = MatrixXd::Random(n, start);
    QR_orthogonalisation(mat);

    mat.conservativeResize(Eigen::NoChange, n);
    mat.rightCols(sub) = MatrixXd::Random(n, sub);

    SECTION("MGS")
    {
        MGS_orthogonalisation(mat, start);
        check_orthogonality(mat);
    }

    SECTION("GS")
    {
        GS_orthogonalisation(mat, start);
        check_orthogonality(mat);
    }

    SECTION("twice_is_enough")
    {
        twice_is_enough_orthogonalisation(mat, start);
        check_orthogonality(mat);
    }

    SECTION("JensWehner")
    {
        JensWehner_orthogonalisation(mat, start);
        check_orthogonality(mat);
    }
}