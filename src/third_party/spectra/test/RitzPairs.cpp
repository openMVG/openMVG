#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <iostream>
#include <Spectra/MatOp/DenseGenMatProd.h>
#include <Spectra/Util/SelectionRule.h>
#include <Spectra/LinAlg/SearchSpace.h>
#include <Spectra/LinAlg/RitzPairs.h>
#include <Spectra/LinAlg/Orthogonalization.h>

using namespace Spectra;

#include "catch.hpp"

using Matrix = Eigen::MatrixXd;
using Vector = Eigen::VectorXd;
using ComplexMatrix = Eigen::MatrixXcd;
using ComplexVector = Eigen::VectorXcd;
using SpMatrix = Eigen::SparseMatrix<double>;
using Index = Eigen::Index;

TEST_CASE("compute_eigen_pairs", "[RitzPairs]")
{
    Matrix A = Eigen::MatrixXd::Random(10, 10);
    Matrix B = A + A.transpose();
    DenseGenMatProd<double> op(B);

    SearchSpace<double> space;
    Matrix initial_space = Matrix::Random(10, 3);
    Spectra::twice_is_enough_orthogonalisation(initial_space);
    space.initialize_search_space(initial_space);
    space.update_operator_basis_product(op);

    RitzPairs<double> ritzpair;
    ritzpair.compute_eigen_pairs(space);

    SECTION("Largest Magnitude")
    {
    }
}

TEST_CASE("Convergence", "[RitzPairs]")
{
    Matrix A = Eigen::MatrixXd::Random(10, 10);
    Matrix B = A + A.transpose();
    DenseGenMatProd<double> op(B);
}
