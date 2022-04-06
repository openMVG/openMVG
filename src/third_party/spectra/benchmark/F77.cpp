#include <Eigen/Core>
#include <iostream>
#include "timer.h"
#include "ArpackFun.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::MatrixXcd;
using Eigen::VectorXcd;
using Eigen::Lower;
typedef Eigen::Map<VectorXd> MapVec;

void eigs_sym_F77(MatrixXd &M, VectorXd &init_resid, int k, int m,
                  double &time_used, double &prec_err, int &nops)
{
    double start, end;
    prec_err = -1.0;
    start = get_wall_time();

    // Begin ARPACK
    //
    // Initial value of ido
    int ido = 0;
    // 'I' means standard eigen value problem, A * x = lambda * x
    char bmat = 'I';
    // dimension of A (n by n)
    int n = M.rows();
    // Specify selection criteria
    // "LM": largest magnitude
    char which[3] = {'L', 'M', '\0'};
    // Number of eigenvalues requested
    int nev = k;
    // Precision
    double tol = 1e-10;
    // Residual vector
    double *resid = new double[n]();
    std::copy(init_resid.data(), init_resid.data() + n, resid);
    // Number of Ritz values used
    int ncv = m;
    // Vector of eigenvalues
    VectorXd evals(nev);
    // Matrix of eigenvectors
    MatrixXd evecs(n, ncv);

    // Store final results of eigenvectors
    // double *V = new double[n * ncv]();
    double *V = evecs.data();
    // Leading dimension of V, required by FORTRAN
    int ldv = n;
    // Control parameters
    int *iparam = new int[11]();
    iparam[1 - 1] = 1;     // ishfts
    iparam[3 - 1] = 1000;  // maxitr
    iparam[7 - 1] = 1;     // mode
    // Some pointers
    int *ipntr = new int[11]();
    /* workd has 3 columns.
     * ipntr[2] - 1 ==> first column to store B * X,
     * ipntr[1] - 1 ==> second to store Y,
     * ipntr[0] - 1 ==> third to store X. */
    double *workd = new double[3 * n]();
    int lworkl = ncv * (ncv + 8);
    double *workl = new double[lworkl]();
    // Error flag. 0 means random initialization,
    // otherwise using resid as initial value
    int info = 1;

    saupd(ido, bmat, n, which,
          nev, tol, resid,
          ncv, V, ldv,
          iparam, ipntr, workd,
          workl, lworkl, info);
    // ido == -1 or ido == 1 means more iterations needed
    while (ido == -1 || ido == 1)
    {
        MapVec vec_in(&workd[ipntr[0] - 1], n);
        MapVec vec_out(&workd[ipntr[1] - 1], n);
        vec_out.noalias() = M.selfadjointView<Lower>() * vec_in;

        saupd(ido, bmat, n, which,
              nev, tol, resid,
              ncv, V, ldv,
              iparam, ipntr, workd,
              workl, lworkl, info);
    }

    // info > 0 means warning, < 0 means error
    if (info > 0)
        std::cout << "warnings occured" << std::endl;
    if (info < 0)
    {
        delete[] workl;
        delete[] workd;
        delete[] ipntr;
        delete[] iparam;
        delete[] resid;

        std::cout << "errors occured" << std::endl;
        end = get_wall_time();
        time_used = (end - start) * 1000;

        return;
    }

    // Retrieve results
    //
    // Whether to calculate eigenvectors or not.
    bool rvec = true;
    // 'A' means to calculate Ritz vectors
    // 'P' to calculate Schur vectors
    char howmny = 'A';
    // Vector of eigenvalues
    double *d = evals.data();
    // Used to store results, will use V instead.
    double *Z = V;
    // Leading dimension of Z, required by FORTRAN
    int ldz = n;
    // Shift
    double sigma = 0;
    // Error information
    int ierr = 0;

    // Number of converged eigenvalues
    int nconv = 0;
    // Number of iterations
    int niter = 0;

    // Use seupd() to retrieve results
    seupd(rvec, howmny, d,
          Z, ldz, sigma, bmat,
          n, which, nev, tol,
          resid, ncv, V, ldv,
          iparam, ipntr, workd, workl,
          lworkl, ierr);

    // Obtain 'nconv' converged eigenvalues
    nconv = iparam[5 - 1];
    // 'niter' number of iterations
    niter = iparam[9 - 1];

    // Free memory of temp arrays
    delete[] workl;
    delete[] workd;
    delete[] ipntr;
    delete[] iparam;
    delete[] resid;

    // ierr < 0 means error
    if (ierr < 0)
    {
        std::cout << "errors occured" << std::endl;
        end = get_wall_time();
        time_used = (end - start) * 1000;

        return;
    }

    /* std::cout << "computed eigenvalues D = \n" << evals.transpose() << std::endl;
    std::cout << "first 5 rows of computed eigenvectors U = \n" <<
    evecs.topLeftCorner(5, nconv) << std::endl;
    std::cout << "nconv = " << nconv << std::endl;
    std::cout << "nops = " << niter << std::endl; */

    end = get_wall_time();
    time_used = (end - start) * 1000;
    MatrixXd err = M * evecs.leftCols(nev) - evecs.leftCols(nev) * evals.asDiagonal();
    prec_err = err.cwiseAbs().maxCoeff();
    nops = niter;
}

void eigs_gen_F77(MatrixXd &M, VectorXd &init_resid, int k, int m,
                  double &time_used, double &prec_err, int &nops)
{
    double start, end;
    prec_err = -1.0;
    start = get_wall_time();

    // Begin ARPACK
    //
    // Initial value of ido
    int ido = 0;
    // 'I' means standard eigen value problem, A * x = lambda * x
    char bmat = 'I';
    // dimension of A (n by n)
    int n = M.rows();
    // Specify selection criteria
    // "LM": largest magnitude
    char which[3] = {'L', 'M', '\0'};
    // Number of eigenvalues requested
    int nev = k;
    // Precision
    double tol = 1e-10;
    // Residual vector
    double *resid = new double[n]();
    std::copy(init_resid.data(), init_resid.data() + n, resid);
    // Number of Ritz values used
    int ncv = m;
    // Vector of eigenvalues
    VectorXd evals_re(nev + 1);
    VectorXd evals_im(nev + 1);
    // Matrix of eigenvectors
    MatrixXd evecs(n, ncv);

    // Store final results of eigenvectors
    // double *V = new double[n * ncv]();
    double *V = evecs.data();
    // Leading dimension of V, required by FORTRAN
    int ldv = n;
    // Control parameters
    int *iparam = new int[11]();
    iparam[1 - 1] = 1;     // ishfts
    iparam[3 - 1] = 1000;  // maxitr
    iparam[7 - 1] = 1;     // mode
    // Some pointers
    int *ipntr = new int[14]();
    /* workd has 3 columns.
     * ipntr[2] - 1 ==> first column to store B * X,
     * ipntr[1] - 1 ==> second to store Y,
     * ipntr[0] - 1 ==> third to store X. */
    double *workd = new double[3 * n]();
    int lworkl = 3 * ncv * ncv + 6 * ncv;
    double *workl = new double[lworkl]();
    // Error flag. 0 means random initialization,
    // otherwise using resid as initial value
    int info = 1;

    naupd(ido, bmat, n, which,
          nev, tol, resid,
          ncv, V, ldv,
          iparam, ipntr, workd,
          workl, lworkl, info);
    // ido == -1 or ido == 1 means more iterations needed
    while (ido == -1 || ido == 1)
    {
        MapVec vec_in(&workd[ipntr[0] - 1], n);
        MapVec vec_out(&workd[ipntr[1] - 1], n);
        vec_out.noalias() = M * vec_in;

        naupd(ido, bmat, n, which,
              nev, tol, resid,
              ncv, V, ldv,
              iparam, ipntr, workd,
              workl, lworkl, info);
    }

    // info > 0 means warning, < 0 means error
    if (info > 0)
        std::cout << "warnings occured" << std::endl;
    if (info < 0)
    {
        delete[] workl;
        delete[] workd;
        delete[] ipntr;
        delete[] iparam;
        delete[] resid;

        std::cout << "errors occured" << std::endl;
        end = get_wall_time();
        time_used = (end - start) * 1000;

        return;
    }

    // Retrieve results
    //
    // Whether to calculate eigenvectors or not.
    bool rvec = true;
    // 'A' means to calculate Ritz vectors
    // 'P' to calculate Schur vectors
    char howmny = 'A';
    // Vector of eigenvalues
    double *dr = evals_re.data();
    double *di = evals_im.data();
    // Used to store results, will use V instead.
    double *Z = V;
    // Leading dimension of Z, required by FORTRAN
    int ldz = n;
    // Shift
    double sigmar = 0;
    double sigmai = 0;
    // Working space
    double *workv = new double[3 * ncv]();
    // Error information
    int ierr = 0;

    // Number of converged eigenvalues
    int nconv = 0;
    // Number of iterations
    int niter = 0;

    // Use seupd() to retrieve results
    neupd(rvec, howmny, dr, di,
          Z, ldz, sigmar, sigmai, workv,
          bmat, n, which, nev, tol,
          resid, ncv, V, ldv, iparam,
          ipntr, workd, workl, lworkl, ierr);

    // Obtain 'nconv' converged eigenvalues
    nconv = iparam[5 - 1];
    // 'niter' number of iterations
    niter = iparam[9 - 1];

    // Free memory of temp arrays
    delete[] workv;
    delete[] workl;
    delete[] workd;
    delete[] ipntr;
    delete[] iparam;
    delete[] resid;

    // ierr < 0 means error
    if (ierr < 0)
    {
        std::cout << "errors occured" << std::endl;
        end = get_wall_time();
        time_used = (end - start) * 1000;

        return;
    }

    /* VectorXcd evals(evals_re.size());
    evals.real() = evals_re;
    evals.imag() = evals_im;
    std::cout << "computed eigenvalues D = \n" << evals << std::endl;
    std::cout << "first 5 rows of computed eigenvectors U = \n" <<
        evecs.topLeftCorner(5, nconv) << std::endl;
    std::cout << "nconv = " << nconv << std::endl;
    std::cout << "nops = " << niter << std::endl; */

    end = get_wall_time();
    time_used = (end - start) * 1000;
    nops = niter;
}
