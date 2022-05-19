
#include <Eigen/Dense>
#include <Spectra/DavidsonSymEigsSolver.h>
#include <Spectra/MatOp/DenseSymMatProd.h>
#include <iostream>

using namespace Spectra;

int main()
{
    Eigen::Index n = 1000;
    Eigen::MatrixXd mat = 0.03 * Eigen::MatrixXd::Random(n, n);
    Eigen::MatrixXd mat1 = mat + mat.transpose();
    for (Eigen::Index i = 0; i < n; i++)
    {
        mat1(i, i) += i + 1;
    }

    DenseSymMatProd<double> op_dense(mat1);  // Create the Matrix Product operation

    Eigen::Index num_of_eigenvalues = 5;
    DavidsonSymEigsSolver<DenseSymMatProd<double>> solver(op_dense, num_of_eigenvalues);  // Create Solver
    Eigen::Index iterations = 100;
    double tolerance = 1e-3;
    int nconv = solver.compute(SortRule::LargestAlge, iterations, tolerance);

    // Retrieve results
    Eigen::VectorXd evalues;
    if (solver.info() == CompInfo::Successful)
    {
        evalues = solver.eigenvalues();

        std::cout << nconv << " Eigenvalues found:\n"
                  << evalues << std::endl;
    }
    else
    {
        std::cout << "Calculation failed" << std::endl;
    }
    return 0;
}