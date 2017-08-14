// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Chris Sweeney (cmsweeney@cs.ucsb.edu)
// Copyright (c) 2016 Pierre Moulon

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_NUMERIC_L1_SOLVER_ADMM_HPP
#define OPENMVG_NUMERIC_L1_SOLVER_ADMM_HPP

#include <Eigen/Core>
#ifdef EIGEN_MPL2_ONLY
#include <Eigen/SparseLU>
#else
#include <Eigen/Cholesky>
#include <Eigen/SparseCholesky>
#endif


#include <algorithm>
#include <iostream>
#include <string>

namespace openMVG {

// These are template overrides that allow the sparse linear solvers to work
// with sparse or dense matrices. The sparseView() method is not implemented for
// Eigen::SparseMatrix.
namespace l1_solver_internal {

template <typename Linear_SolverT>
inline void Compute
(
  const Eigen::SparseMatrix<double>& spd_mat,
  Linear_SolverT * linear_solver
)
{
  linear_solver->compute(spd_mat);
}

template <typename Linear_SolverT>
inline void Compute
(
  const Eigen::MatrixXd& spd_mat,
  Linear_SolverT * linear_solver
)
{
  linear_solver->compute(spd_mat.sparseView());
}

}  // namespace l1_solver_internal

// A L1 norm approximation solver. This class will attempt to solve the
// problem: || A * x - b || under L1-norm (as opposed to L2 i.e. "least-squares"
// norm). This problem can be solved with the alternating direction method of
// multipliers (ADMM) as a least unsquared deviations minimizer. A full
// description of the method, including how to use ADMM for L1 minimization can
// be found in "Distributed Optimization and Statistical Learning via the
// Alternating Direction Method of Multipliers" by Boyd et al, Foundations and
// Trends in Machine Learning (2012). The paper can be found at:
//   https://web.stanford.edu/~boyd/papers/pdf/admm_distr_stats.pdf
//
// ADMM can be much faster than interior point methods but convergence may be
// slower. Generally speaking, ADMM solvers converge to good solutions in only a
// few number of iterations, but can spend many iterations subsequently refining
// the solution to obtain the global optimum. The speed improvements are because
// the matrix A only needs to be factorized (by Cholesky decomposition) once, as
// opposed to every iteration.
//
// This implementation is based off of the code found at:
//   https://web.stanford.edu/~boyd/papers/admm/least_abs_deviations/lad.html
template <class MatrixType>
class L1Solver {
 public:
  struct Options {
    int max_num_iterations = 1000;
    // Rho is the augmented Lagrangian parameter.
    double rho = 1.0;
    // Alpha is the over-relaxation parameter (typically between 1.0 and 1.8).
    double alpha = 1.0;

    double absolute_tolerance = 1e-4;
    double relative_tolerance = 1e-2;
  };

  L1Solver
  (
    const Options& options,
    const MatrixType& mat
  )
  : options_(options), a_(mat)
  {
    // Analyze the sparsity pattern once. Only the values of the entries will be
    // changed with each iteration.
    const MatrixType spd_mat = a_.transpose() * a_;
    l1_solver_internal::Compute(spd_mat, &linear_solver_);
  }

  void SetMaxIterations
  (
    const int max_iterations
  )
  {
    options_.max_num_iterations = max_iterations;
  }

  bool Status() const
  {
    return linear_solver_.info() == Eigen::Success;
  }

  // Solves ||Ax - b||_1 for the optimal L1 solution given an initial guess for
  // x. To solve this we introduce an auxiliary variable y such that the
  // solution to:
  //        min   1 * y
  //   s.t. [  A   -I ] [ x ] < [  b ]
  //        [ -A   -I ] [ y ]   [ -b ]
  // which is an equivalent linear program.
  bool Solve
  (
    const Eigen::VectorXd& rhs,
    Eigen::VectorXd* solution
  )
  {
    // Since constructor was called before we check Compute status
    if (linear_solver_.info() != Eigen::Success)
    {
      std::cerr << "Cannot compute the matrix factorization" << std::endl;
      return false;
    }

    Eigen::VectorXd& x = *solution;
    Eigen::VectorXd z(a_.rows()), u(a_.rows());
    z.setZero();
    u.setZero();

    Eigen::VectorXd a_times_x(a_.rows()), z_old(z.size()), ax_hat(a_.rows());
    // Precompute some convergence terms.
    const double rhs_norm = rhs.norm();
    const double primal_abs_tolerance_eps =
      std::sqrt(a_.rows()) * options_.absolute_tolerance;
    const double dual_abs_tolerance_eps =
      std::sqrt(a_.cols()) * options_.absolute_tolerance;

    for (int i = 0; i < options_.max_num_iterations; ++i)
    {
      // Update x.
      x.noalias() = linear_solver_.solve(a_.transpose() * (rhs + z - u));
      a_times_x.noalias() = a_ * x;
      ax_hat.noalias() = options_.alpha * a_times_x;
      ax_hat.noalias() += (1.0 - options_.alpha) * (z + rhs);

      // Update z and set z_old.
      std::swap(z, z_old);
      z.noalias() = Shrinkage(ax_hat - rhs + u, 1.0 / options_.rho);

      // Update u.
      u.noalias() += ax_hat - z - rhs;

      // Compute the convergence terms.
      const double r_norm = (a_times_x - z - rhs).norm();
      const double s_norm =
        (-options_.rho * a_.transpose() * (z - z_old)).norm();
      const double max_norm =
        std::max({a_times_x.norm(), z.norm(), rhs_norm});
      const double primal_eps =
        primal_abs_tolerance_eps + options_.relative_tolerance * max_norm;
      const double dual_eps =
        dual_abs_tolerance_eps +
        options_.relative_tolerance *
          (options_.rho * a_.transpose() * u).norm();

      // Log the result to the screen.
      // std::ostringstream os;
      // os << "Iteration: " << i << "\n"
      //   << "R norm: " << r_norm << "\n"
      //   << "S norm: " << s_norm << "\n"
      //   << "Primal eps: " << primal_eps << "\n"
      //   << "Dual eps: " << dual_eps << std::endl;
      // std::cout << os.str() << std::endl;

      // Determine if the minimizer has converged.
      if (r_norm < primal_eps && s_norm < dual_eps)
      {
        return true;
      }
    }
    return false;
  }

 private:
  Options options_;

  // Matrix A where || Ax - b ||_1 is the problem we are solving.
  MatrixType a_;

  // Cholesky linear solver.
#ifdef EIGEN_MPL2_ONLY
  using Linear_Solver_T = Eigen::SparseLU<Eigen::SparseMatrix<double> >;
#else
  // Since our linear system will be a SPD matrix we can
  // utilize the Cholesky factorization.
  using Linear_Solver_T = Eigen::SimplicialLLT<Eigen::SparseMatrix<double> >;
#endif
  Linear_Solver_T linear_solver_;

  Eigen::VectorXd Shrinkage
  (
    const Eigen::VectorXd& vec, const double kappa
  ) const
  {
    Eigen::ArrayXd zero_vec(vec.size());
    zero_vec.setZero();
    return zero_vec.max( vec.array() - kappa) -
           zero_vec.max(-vec.array() - kappa);
  }
};

}  // namespace openMVG

#endif  // OPENMVG_NUMERIC_L1_SOLVER_ADMM_HPP
