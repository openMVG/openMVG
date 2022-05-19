#include <Eigen/Core>
#include <iostream>
#include "timer.h"

#include <Spectra/SymEigsSolver.h>
#include <Spectra/GenEigsSolver.h>

using namespace Spectra;

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::MatrixXcd;
using Eigen::VectorXcd;

void eigs_sym_Cpp(MatrixXd &M, VectorXd &init_resid, int k, int m,
                  double &time_used, double &prec_err, int &nops)
{
    double start, end;
    start = get_wall_time();

    DenseSymMatProd<double> op(M);
    SymEigsSolver<DenseSymMatProd<double>> eigs(op, k, m);
    eigs.init(init_resid.data());

    int nconv = eigs.compute(SortRule::LargestMagn);
    int niter = eigs.num_iterations();
    nops = eigs.num_operations();

    VectorXd evals = eigs.eigenvalues();
    MatrixXd evecs = eigs.eigenvectors();

    /* std::cout << "computed eigenvalues D = \n" << evals.transpose() << std::endl;
    std::cout << "first 5 rows of computed eigenvectors U = \n" << evecs.topRows<5>() << std::endl;
    std::cout << "nconv = " << nconv << std::endl;
    std::cout << "niter = " << niter << std::endl;
    std::cout << "nops = " << nops << std::endl; */

    end = get_wall_time();
    time_used = (end - start) * 1000;

    MatrixXd err = M * evecs - evecs * evals.asDiagonal();
    prec_err = err.cwiseAbs().maxCoeff();
}

void eigs_gen_Cpp(MatrixXd &M, VectorXd &init_resid, int k, int m,
                  double &time_used, double &prec_err, int &nops)
{
    double start, end;
    start = get_wall_time();

    DenseGenMatProd<double> op(M);
    GenEigsSolver<DenseGenMatProd<double>> eigs(op, k, m);
    eigs.init(init_resid.data());

    int nconv = eigs.compute(SortRule::LargestMagn);
    int niter = eigs.num_iterations();
    nops = eigs.num_operations();

    VectorXcd evals = eigs.eigenvalues();
    MatrixXcd evecs = eigs.eigenvectors();

    /* std::cout << "computed eigenvalues D = \n" << evals.transpose() << std::endl;
    std::cout << "first 5 rows of computed eigenvectors U = \n" << evecs.topRows<5>() << std::endl;
    std::cout << "nconv = " << nconv << std::endl;
    std::cout << "niter = " << niter << std::endl;
    std::cout << "nops = " << nops << std::endl;

    MatrixXcd err = M * evecs - evecs * evals.asDiagonal();
    std::cout << "||AU - UD||_inf = " << err.array().abs().maxCoeff() << std::endl; */

    end = get_wall_time();
    time_used = (end - start) * 1000;

    MatrixXcd err = M * evecs - evecs * evals.asDiagonal();
    prec_err = err.cwiseAbs().maxCoeff();
}
