#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <iostream>
#include <Spectra/MatOp/DenseGenMatProd.h>
#include <Spectra/LinAlg/SearchSpace.h>
#include <Spectra/LinAlg/Orthogonalization.h>

using namespace Spectra;

#include "catch.hpp"

using Matrix = Eigen::MatrixXd;
using Vector = Eigen::VectorXd;
using ComplexMatrix = Eigen::MatrixXcd;
using ComplexVector = Eigen::VectorXcd;
using SpMatrix = Eigen::SparseMatrix<double>;
using Index = Eigen::Index;

TEST_CASE("CompleteSearchSpace", "[SearchSpace]")
{
    SearchSpace<double> space;
    Matrix initial_space = Matrix::Random(10, 3);
    Spectra::twice_is_enough_orthogonalisation(initial_space);
    space.initialize_search_space(initial_space);
    REQUIRE(space.basis_vectors().cols() == 3);
    REQUIRE(space.operator_basis_product().cols() == 0);

    Matrix A = Eigen::MatrixXd::Random(10, 10);
    Matrix B = A + A.transpose();
    DenseGenMatProd<double> op(B);

    space.update_operator_basis_product(op);
    REQUIRE(space.basis_vectors().cols() == 3);
    REQUIRE(space.operator_basis_product().cols() == 3);
    REQUIRE(space.operator_basis_product().isApprox(B * initial_space));

    Matrix append_space = Matrix::Random(10, 3);
    space.extend_basis(append_space);
    REQUIRE((space.basis_vectors().transpose() * space.basis_vectors()).isIdentity(1e-12));
    REQUIRE(space.basis_vectors().cols() == 6);
    REQUIRE(space.operator_basis_product().cols() == 3);
    space.update_operator_basis_product(op);
    REQUIRE(space.operator_basis_product().cols() == 6);

    RitzPairs<double> ritzpair;
    ritzpair.compute_eigen_pairs(space);
    REQUIRE(ritzpair.size() == 6);
    space.restart(ritzpair, 2);

    REQUIRE(space.basis_vectors().cols() == 2);
    REQUIRE(space.operator_basis_product().cols() == 2);
}
