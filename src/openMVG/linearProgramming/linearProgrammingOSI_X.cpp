// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/linearProgramming/linearProgrammingOSI_X.hpp"
#include <assert.h>
#include <cstddef>
#include "CoinPackedVector.hpp"
#include "OsiClpSolverInterface.hpp"
#include "openMVG/numeric/eigen_alias_definition.hpp"

namespace openMVG   {
namespace linearProgramming  {

OSI_X_SolverWrapper::OSI_X_SolverWrapper(int nbParams) : LP_Solver(nbParams)
{
  si.reset(new OsiClpSolverInterface);
  si->setLogLevel(0);
}

bool OSI_X_SolverWrapper::setup(const LP_Constraints & cstraints) //cstraints <-> constraints
{
  if (!si)
  {
    return false;
  }
  assert(nbParams_ == cstraints.nbParams_);

  const unsigned int NUMVAR = cstraints.constraint_mat_.cols();
  std::vector<double>
    col_lb(NUMVAR), // the column lower bounds
    col_ub(NUMVAR); // the column upper bounds

  this->nbParams_ = NUMVAR;

  si->setObjSense( ((cstraints.bminimize_) ? 1 : -1) );

  const Mat & A = cstraints.constraint_mat_;

  //Equality constraint will be done by two constraints due to the API limitation ( >= & <=).
  const size_t nbLine = A.rows() +
    std::count(cstraints.vec_sign_.begin(), cstraints.vec_sign_.end(), LP_Constraints::LP_EQUAL);

  // Define default lower and upper bound to [-inf, inf]
  std::vector<double>
    row_lb(nbLine, -si->getInfinity()), // the row lower bounds
    row_ub(nbLine, si->getInfinity()); // the row upper bounds

  std::unique_ptr<CoinPackedMatrix> matrix(new CoinPackedMatrix(false,0,0));
  matrix->setDimensions(0, NUMVAR);

  //-- Add row-wise constraint
  size_t indexRow = 0;
  for (int i=0; i < A.rows(); ++i)
  {
    const Vec temp = A.row(i);
    CoinPackedVector row;

    if (cstraints.vec_sign_[i] == LP_Constraints::LP_EQUAL ||
        cstraints.vec_sign_[i] == LP_Constraints::LP_LESS_OR_EQUAL)
    {
      const int coef = 1;
      for ( int j = 0; j < A.cols(); ++j )
      {
        row.insert(j, coef * temp.data()[j]);
      }
      row_ub[indexRow] = coef * cstraints.constraint_objective_(i);
      matrix->appendRow(row);
      ++indexRow;
    }

    if (cstraints.vec_sign_[i] == LP_Constraints::LP_EQUAL ||
        cstraints.vec_sign_[i] == LP_Constraints::LP_GREATER_OR_EQUAL)
    {
      const int coef = -1;
      for ( int j = 0; j < A.cols(); ++j )
      {
        row.insert(j, coef * temp.data()[j]);
      }
      row_ub[indexRow] = coef * cstraints.constraint_objective_(i);
      matrix->appendRow(row);
      ++indexRow;
    }
  }

  //-- Setup bounds for all the parameters
  if (cstraints.vec_bounds_.size() == 1)
  {
    // Setup the same bound for all the parameters
    std::fill(col_lb.begin(), col_lb.end(), cstraints.vec_bounds_[0].first);
    std::fill(col_ub.begin(), col_ub.end(), cstraints.vec_bounds_[0].second);
  }
  else // each parameter have its own bounds
  {
    for (int i=0; i < this->nbParams_; ++i)
    {
      col_lb[i] = cstraints.vec_bounds_[i].first;
      col_ub[i] = cstraints.vec_bounds_[i].second;
    }
  }

  si->loadProblem(*matrix, &col_lb[0], &col_ub[0], cstraints.vec_cost_.empty() ? nullptr : &cstraints.vec_cost_[0], &row_lb[0], &row_ub[0] );

  return true;
}

bool OSI_X_SolverWrapper::setup(const LP_Constraints_Sparse & cstraints) //cstraints <-> constraints
{
  if (!si)
  {
    return false;
  }
  assert(nbParams_ == cstraints.nbParams_);

  const int NUMVAR = cstraints.constraint_mat_.cols();
  std::vector<double>
    col_lb(NUMVAR), // the column lower bounds
    col_ub(NUMVAR); // the column upper bounds

  this->nbParams_ = NUMVAR;

  si->setObjSense( ((cstraints.bminimize_) ? 1 : -1) );

  const sRMat & A = cstraints.constraint_mat_;

  //Equality constraint will be done by two constraints due to the API limitation (>= & <=)
  const size_t nbLine = A.rows() +
    std::count(cstraints.vec_sign_.begin(), cstraints.vec_sign_.end(), LP_Constraints::LP_EQUAL);

  // Define default lower and upper bound to [-inf, inf]
  std::vector<double>
    row_lb(nbLine, -si->getInfinity()), // the row lower bounds
    row_ub(nbLine, si->getInfinity()); // the row upper bounds

  std::unique_ptr<CoinPackedMatrix> matrix(new CoinPackedMatrix(false,0,0));
  matrix->setDimensions(0, NUMVAR);

  //-- Add row-wise constraint
  size_t rowindex = 0;
  for (int i=0; i < A.rows(); ++i)
  {
    std::vector<int> vec_colno;
    std::vector<double> vec_value;
    for (sRMat::InnerIterator it(A,i); it; ++it)
    {
      vec_colno.push_back(it.col());
      vec_value.push_back(it.value());
    }

    if ( cstraints.vec_sign_[i] == LP_Constraints::LP_EQUAL ||
         cstraints.vec_sign_[i] == LP_Constraints::LP_LESS_OR_EQUAL )
    {
      const int coef = 1;
      row_ub[rowindex] = coef * cstraints.constraint_objective_(i);
      matrix->appendRow( vec_colno.size(),
                   &vec_colno[0],
                   &vec_value[0] );
      ++rowindex;
    }

    if ( cstraints.vec_sign_[i] == LP_Constraints::LP_EQUAL ||
         cstraints.vec_sign_[i] == LP_Constraints::LP_GREATER_OR_EQUAL )
    {
      const int coef = -1;
      for (auto &  iter_val : vec_value)
      {
        iter_val *= coef;
      }
      row_ub[rowindex] = coef * cstraints.constraint_objective_(i);
      matrix->appendRow( vec_colno.size(),
                   &vec_colno[0],
                   &vec_value[0] );
      ++rowindex;
    }
  }

  //-- Setup bounds for all the parameters
  if (cstraints.vec_bounds_.size() == 1)
  {
    // Setup the same bound for all the parameters
    std::fill(col_lb.begin(), col_lb.end(), cstraints.vec_bounds_[0].first);
    std::fill(col_ub.begin(), col_ub.end(), cstraints.vec_bounds_[0].second);
  }
  else  // each parameter have its own bounds
  {
    for (int i=0; i < this->nbParams_; ++i)
    {
      col_lb[i] = cstraints.vec_bounds_[i].first;
      col_ub[i] = cstraints.vec_bounds_[i].second;
    }
  }

  si->loadProblem(
    *matrix,
    &col_lb[0],
    &col_ub[0],
    cstraints.vec_cost_.empty() ? nullptr : &cstraints.vec_cost_[0],
    &row_lb[0],
    &row_ub[0]);

  return true;
}


bool OSI_X_SolverWrapper::solve()
{
  //-- Compute solution
  if ( si )
  {
    si->getModelPtr()->setPerturbation(50);
    si->initialSolve();
    return si->isProvenOptimal();
  }
  return false;
}


bool OSI_X_SolverWrapper::getSolution(std::vector<double> & estimatedParams)
{
  if ( si )
  {
    const int n = si->getNumCols();
    std::memcpy(&estimatedParams[0], si->getColSolution(), n * sizeof(double));
    return true;
  }
  return false;
}

} // namespace linearProgramming
} // namespace openMVG
