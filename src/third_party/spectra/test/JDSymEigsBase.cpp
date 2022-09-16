#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <iostream>
#include <random>  // Requires C++ 11

#include <Spectra/JDSymEigsBase.h>
#include <Spectra/MatOp/DenseGenMatProd.h>

using namespace Spectra;

#include "catch.hpp"

using Matrix = Eigen::MatrixXd;
using Vector = Eigen::VectorXd;
using ComplexMatrix = Eigen::MatrixXcd;
using ComplexVector = Eigen::VectorXcd;
using SpMatrix = Eigen::SparseMatrix<double>;
using Index = Eigen::Index;
template <typename OpType>
class JDMock : public JDSymEigsBase<JDMock<OpType>, OpType>
{
public:
    JDMock(OpType& op, Index nev) :
        JDSymEigsBase<JDMock<OpType>, OpType>(op, nev) {}
    Matrix setup_initial_search_space(SortRule) const
    {
        return Matrix::Zero(0, 0);
    }

    Matrix calculate_correction_vector() const
    {
        return Matrix::Zero(0, 0);
    }
};

TEST_CASE("Constructing JDSymObject", "[eigs_gen]")
{
    const Matrix A = Eigen::MatrixXd::Random(10, 10);
    DenseGenMatProd<double> op(A);
    JDMock<DenseGenMatProd<double>> eigs(op, 5);
    REQUIRE(eigs.num_iterations() == 0);
    REQUIRE(eigs.info() == Spectra::CompInfo::NotComputed);
}