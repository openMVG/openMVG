// Copyright (c) 2012 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef MIMATTE_LINEAR_PROGRAMMING_INTERFACE_OSICLP_H_
#define MIMATTE_LINEAR_PROGRAMMING_INTERFACE_OSICLP_H_

#include "OsiClpSolverInterface.hpp"
#ifdef OPENMVG_HAVE_MOSEK
#include "OsiMskSolverInterface.hpp"
#endif

#include "openMVG/numeric/numeric.h"
#include "openMVG/linearProgramming/linearProgrammingInterface.hpp"

#include "CoinPackedMatrix.hpp"
#include "CoinPackedVector.hpp"

#include <vector>

namespace openMVG   {
namespace linearProgramming  {

/// OSI_X wrapper for the LP_Solver
template<typename SOLVERINTERFACE>
class OSI_X_SolverWrapper : public LP_Solver
{
public :
  OSI_X_SolverWrapper(int nbParams);

  ~OSI_X_SolverWrapper();

  //--
  // Inherited functions :
  //--

  bool setup(const LP_Constraints & constraints);
  bool setup(const LP_Constraints_Sparse & constraints);

  bool solve();

  bool getSolution(std::vector<double> & estimatedParams);

private :
  SOLVERINTERFACE *si;
};


typedef OSI_X_SolverWrapper<OsiClpSolverInterface> OSI_CLP_SolverWrapper;
#ifdef OPENMVG_HAVE_MOSEK
typedef OSI_X_SolverWrapper<OsiMskSolverInterface> OSI_MOSEK_SolverWrapper;
#endif // OPENMVG_HAVE_MOSEK



template<typename SOLVERINTERFACE>
OSI_X_SolverWrapper<SOLVERINTERFACE>::OSI_X_SolverWrapper(int nbParams) : LP_Solver(nbParams)
{
  si = new SOLVERINTERFACE;
  si->setLogLevel(0);
}

#ifdef OPENMVG_HAVE_MOSEK
template<>
OSI_X_SolverWrapper<OsiMskSolverInterface>::OSI_X_SolverWrapper(int nbParams) : LP_Solver(nbParams)
{
  si = new OsiMskSolverInterface();
  //si->setLogLevel(0);
}
#endif // OPENMVG_HAVE_MOSEK

template<typename SOLVERINTERFACE>
OSI_X_SolverWrapper<SOLVERINTERFACE>::~OSI_X_SolverWrapper()
{
  // Memory cleaning.
  if ( si != NULL )
  {
    delete si;
    si = NULL;
  }
}

template<typename SOLVERINTERFACE>
bool OSI_X_SolverWrapper<SOLVERINTERFACE>::setup(const LP_Constraints & cstraints) //cstraints <-> constraints
{
  bool bOk = true;
  if ( si == NULL )
  {
    return false;
  }
  assert(_nbParams == cstraints._nbParams);

  const unsigned int NUMVAR = cstraints._constraintMat.cols();
  std::vector<double> col_lb(NUMVAR);//the column lower bounds
  std::vector<double> col_ub(NUMVAR);//the column upper bounds

  this->_nbParams = NUMVAR;

  si->setObjSense( ((cstraints._bminimize) ? 1 : -1) );

  const Mat & A = cstraints._constraintMat;

  //Equality constraint will be done by two constraints due to the API limitation ( >= & <=).
  const size_t nbLine = A.rows() +
    std::count(cstraints._vec_sign.begin(), cstraints._vec_sign.end(), LP_Constraints::LP_EQUAL);

  std::vector<double> row_lb(nbLine);//the row lower bounds
  std::vector<double> row_ub(nbLine);//the row upper bounds

  CoinPackedMatrix * matrix = new CoinPackedMatrix(false,0,0);
  matrix->setDimensions(0, NUMVAR);

  //-- Add row-wise constraint
  size_t indexRow = 0;
  for (int i=0; i < A.rows(); ++i)
  {
    Vec temp = A.row(i);

    CoinPackedVector row;
    if ( cstraints._vec_sign[i] == LP_Constraints::LP_EQUAL || cstraints._vec_sign[i] == LP_Constraints::LP_LESS_OR_EQUAL )
    {
      int coef = 1;
      for ( int j = 0; j < A.cols() ; j++ )
      {
        row.insert(j, coef * temp.data()[j]);
      }
      row_lb[indexRow] = -1.0 * si->getInfinity();
      row_ub[indexRow] = coef * cstraints._Cst_objective(i);
      matrix->appendRow(row);
      indexRow++;
    }
    if ( cstraints._vec_sign[i] == LP_Constraints::LP_EQUAL || cstraints._vec_sign[i] == LP_Constraints::LP_GREATER_OR_EQUAL )
    {
      int coef = -1;
      for ( int j = 0; j < A.cols() ; j++ )
      {
	      row.insert(j, coef * temp.data()[j]);
      }
      row_lb[indexRow] = -1.0 * si->getInfinity();
      row_ub[indexRow] = coef * cstraints._Cst_objective(i);
      matrix->appendRow(row);
      indexRow++;
    }
  }

  //-- Setup bounds for all the parameters
  if (cstraints._vec_bounds.size() == 1)
  {
    // Setup the same bound for all the parameters
    for (int i=0; i < this->_nbParams; ++i)
    {
      col_lb[i] = cstraints._vec_bounds[0].first;
      col_ub[i] = cstraints._vec_bounds[0].second;
    }
  }
  else // each parameter have it's own bounds
  {
    for (int i=0; i < this->_nbParams; ++i)
    {
      col_lb[i] = cstraints._vec_bounds[i].first;
      col_ub[i] = cstraints._vec_bounds[i].second;
    }
  }

  si->loadProblem(*matrix, &col_lb[0], &col_ub[0], cstraints._vec_cost.empty() ? NULL : &cstraints._vec_cost[0], &row_lb[0], &row_ub[0] );

  delete matrix;

  return bOk;
}

template<typename SOLVERINTERFACE>
bool OSI_X_SolverWrapper<SOLVERINTERFACE>::setup(const LP_Constraints_Sparse & cstraints) //cstraints <-> constraints
{
  bool bOk = true;
  if ( si == NULL )
  {
    return false;
  }
  assert(_nbParams == cstraints._nbParams);

  const int NUMVAR = cstraints._constraintMat.cols();
  std::vector<double> col_lb(NUMVAR);//the column lower bounds
  std::vector<double> col_ub(NUMVAR);//the column upper bounds

  this->_nbParams = NUMVAR;

  si->setObjSense( ((cstraints._bminimize) ? 1 : -1) );

  const sRMat & A = cstraints._constraintMat;

  //Equality constraint will be done by two constraints due to the API limitation (>= & <=)
  const size_t nbLine = A.rows() +
    std::count(cstraints._vec_sign.begin(), cstraints._vec_sign.end(), LP_Constraints::LP_EQUAL);

  std::vector<double> row_lb(nbLine);//the row lower bounds
  std::vector<double> row_ub(nbLine);//the row upper bounds

  CoinPackedMatrix * matrix = new CoinPackedMatrix(false,0,0);
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


    if ( cstraints._vec_sign[i] == LP_Constraints::LP_EQUAL || cstraints._vec_sign[i] == LP_Constraints::LP_LESS_OR_EQUAL )
    {
      int coef = 1;
      row_lb[rowindex] = -1.0 * si->getInfinity();
      row_ub[rowindex] = coef * cstraints._Cst_objective(i);
      matrix->appendRow( vec_colno.size(),
	                 &vec_colno[0],
	                 &vec_value[0] );
      rowindex++;
    }

    if ( cstraints._vec_sign[i] == LP_Constraints::LP_EQUAL || cstraints._vec_sign[i] == LP_Constraints::LP_GREATER_OR_EQUAL )
    {
      int coef = -1;
      for ( std::vector<double>::iterator iter_val = vec_value.begin();
	      iter_val != vec_value.end();
	      iter_val++)
      {
  	    *iter_val *= coef;
      }
      row_lb[rowindex] = -1.0 * si->getInfinity();
      row_ub[rowindex] = coef * cstraints._Cst_objective(i);
      matrix->appendRow( vec_colno.size(),
	                 &vec_colno[0],
	                 &vec_value[0] );
      rowindex++;
    }
  }

  //-- Setup bounds for all the parameters
  if (cstraints._vec_bounds.size() == 1)
  {
    // Setup the same bound for all the parameters
    for (int i=0; i < this->_nbParams; ++i)
    {
      col_lb[i] = cstraints._vec_bounds[0].first;
      col_ub[i] = cstraints._vec_bounds[0].second;
    }
  }
  else  // each parameter have it's own bounds
  {
    for (int i=0; i < this->_nbParams; ++i)
    {
      col_lb[i] = cstraints._vec_bounds[i].first;
      col_ub[i] = cstraints._vec_bounds[i].second;
    }
  }

  si->loadProblem(
    *matrix,
    &col_lb[0],
    &col_ub[0],
    cstraints._vec_cost.empty() ? NULL : &cstraints._vec_cost[0],
    &row_lb[0],
    &row_ub[0]);

  delete matrix;

  return bOk;
}

template<typename SOLVERINTERFACE>
bool OSI_X_SolverWrapper<SOLVERINTERFACE>::solve()
{
  //-- Compute solution
  if ( si != NULL )
  {
    si->getModelPtr()->setPerturbation(50);
    si->initialSolve();
    return si->isProvenOptimal();
  }
  return false;
}

template<typename SOLVERINTERFACE>
bool OSI_X_SolverWrapper<SOLVERINTERFACE>::getSolution(std::vector<double> & estimatedParams)
{
  if ( si != NULL )
  {
    const int n = si->getNumCols();
    memcpy(&estimatedParams[0], si->getColSolution(), n * sizeof(double));
    return true;
  }
  return false;
}

} // namespace linearProgramming
} // namespace openMVG


#endif // MIMATTE_LINEAR_PROGRAMMING_INTERFACE_OSICLP_H_

