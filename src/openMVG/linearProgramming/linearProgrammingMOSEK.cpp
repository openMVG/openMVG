// Copyright (c) 2012 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifdef OPENMVG_HAVE_MOSEK

#include "openMVG/linearProgramming/linearProgrammingMOSEK.hpp"
#include <iostream>

namespace openMVG {
namespace linearProgramming {

using namespace std;

// This function prints log output from MOSEK to the terminal.
static void MSKAPI printstr(void *handle,
                            char str[])
{
  std::cout << str;
}

MSKenv_t     env = NULL;

MOSEK_SolveWrapper::MOSEK_SolveWrapper(int nbParams) : LP_Solver(nbParams)
{
  task = NULL;
  //env  = NULL;

  //-- Initialize MOSEK framework

  // Create the mosek environment.
  /*MSKrescodee r = MSK_makeenv(&env,NULL,NULL,NULL,NULL);
  if ( r!=MSK_RES_OK )  {
    std::cerr << "Cannot create the MOSEK environment" << std::endl;
  }*/
  if (env == NULL)
  {
     MSKrescodee r = MSK_makeenv(&env,NULL,NULL,NULL,NULL);

    // Directs the env log stream to the 'printstr' function.
    if ( r==MSK_RES_OK )
      MSK_linkfunctoenvstream(env,MSK_STREAM_LOG,NULL,printstr);

    // Initialize the environment.
    if ( r==MSK_RES_OK )
      r = MSK_initenv(env);
    else  {
      std::cerr << "Cannot create the MOSEK environment" << std::endl;
    }
  }

}

MOSEK_SolveWrapper::~MOSEK_SolveWrapper()
{
  // Memory cleaning.
  MSK_deletetask(&task);
  //MSK_deleteenv(&env);
}


inline MSKboundkey_enum convertSign(LP_Constraints::eLP_SIGN sign) {

  switch(sign) {
    case LP_Constraints::LP_LESS_OR_EQUAL:    // = 1,  // cst (<=) VAL
      return MSK_BK_UP;
    case LP_Constraints::LP_GREATER_OR_EQUAL: // = 2,  // cst (>=) VAL
      return MSK_BK_LO;
    case LP_Constraints::LP_EQUAL:            // = 3,  // cst (=) VAL
      return MSK_BK_FX;
    case LP_Constraints::LP_FREE:
      return MSK_BK_FR;
    default:
      std::cerr << "Error unknow constraint sign : " << sign << "\n";
  }
}

bool MOSEK_SolveWrapper::setup(const LP_Constraints & cstraints) //cstraints <-> constraints
{
  assert(nbParams_ == cstraints.nbParams_);

  MSK_deletetask(&task);

  int NUMVAR = cstraints.constraint_mat_.cols();
  int NUMCON = cstraints.constraint_mat_.rows();
  int NUMANZ = cstraints.constraint_mat_.cols() * cstraints.constraint_mat_.rows(); //DENSE MATRIX

  MSKrescodee r = MSK_RES_OK;
  if ( r==MSK_RES_OK )
  {
    // Create the optimization task.
    r = MSK_maketask(env,NUMCON,NUMVAR,&task);

    // Directs the log task stream to the 'printstr' function.
    if ( r==MSK_RES_OK )
      MSK_linkfunctotaskstream(task,MSK_STREAM_LOG,NULL,printstr);

    // Give MOSEK an estimate of the size of the input data.
    //This is done to increase the speed of inputting data.
    // However, it is optional.
    if (r == MSK_RES_OK)
      r = MSK_putmaxnumvar(task,NUMVAR);

    if (r == MSK_RES_OK)
      r = MSK_putmaxnumcon(task,NUMCON);

    if (r == MSK_RES_OK)
      r = MSK_putmaxnumanz(task,NUMANZ);

    // Append 'NUMCON' empty constraints. The constraints will initially have no bounds.
    if ( r == MSK_RES_OK )
      r = MSK_append(task,MSK_ACC_CON,NUMCON);

    // Append 'NUMVAR' variables. The variables will initially be fixed at zero (x=0).
    if ( r == MSK_RES_OK )
      r = MSK_append(task,MSK_ACC_VAR,NUMVAR);
  }

  this->nbParams_ = NUMVAR;

  if (cstraints.bminimize_) {
    r = MSK_putobjsense(task, MSK_OBJECTIVE_SENSE_MINIMIZE);
  }
  else  {
    r = MSK_putobjsense(task, MSK_OBJECTIVE_SENSE_MAXIMIZE);
  }

  // Optionally add a constant term to the objective.
  if ( r ==MSK_RES_OK )
    r = MSK_putcfix(task,0.0);

  // Add the objective function if any
  if (!cstraints.vec_cost_.empty())  {
    // Set objective
    for (size_t i = 0; i < cstraints.vec_cost_.size(); ++i)
      MSK_putcj(task, i, cstraints.vec_cost_[i]);
  }

  //Add constraint row by row
  const Mat & A = cstraints.constraint_mat_;

  for (int i=0; i<A.rows() && r == MSK_RES_OK; ++i)
  {
    std::vector<int> vec_colno(A.cols(), 0);
    for (int ii=0; ii<A.cols(); ++ii){
      vec_colno[ii] = ii;
    }

    // Insert a row
    Vec temp = A.row(i);
    r = MSK_putavec(task,
                    MSK_ACC_CON,       // Input row of A.
                    MSKidxt (i),       // Variable row index.
                    vec_colno.size(),  // Number of non-zeros in row(i) A.
                    &vec_colno[0],     // Pointer to row indexes of row i.
                    (double*)temp.data());    // Pointer to Values of row i.
  }

  //Add bound on variables:
  if (cstraints.vec_bounds_.size() == 1) {
    for (size_t i = 0; i < NUMVAR; ++i) {
      if (r == MSK_RES_OK)
        r = MSK_putbound(task,
                         MSK_ACC_VAR, // Put bounds on variables.
                         i,           // Index of variable.
                         convertSign(cstraints.vec_sign_[0]),      // Bound key.
                         cstraints.vec_bounds_[0].first,  // Numerical value of lower bound.
                         cstraints.vec_bounds_[0].second); // Numerical value of upper bound.
    }
  }
  else{
    for (size_t i = 0; i < NUMVAR; ++i) {
    // Set the bounds on variable j.
    //    lowerbound <= x_j <= upper bound
      if (r == MSK_RES_OK)
        r = MSK_putbound(task,
                         MSK_ACC_VAR, // Put bounds on variables.
                         i,           // Index of variable.
                         convertSign(cstraints.vec_sign_[i]),      // Bound key.
                         cstraints.vec_bounds_[i].first,  // Numerical value of lower bound.
                         cstraints.vec_bounds_[i].second); // Numerical value of upper bound.
    // in order to add sparse bounds use: MSK_putboundlist
    }
  }

  // Add bounds on constraint
  for (size_t i=0; i<NUMCON && r==MSK_RES_OK; ++i)
  {
    r = MSK_putbound(task,
                     MSK_ACC_CON, // Put bounds on constraints.
                     i,           // Index of constraint.
                     convertSign(cstraints.vec_sign_[i]),      // Bound key.
                     -MSK_INFINITY,  // Numerical value of lower bound.
                     cstraints.constraint_objective_(i)); // Numerical value of upper bound.
  }

  return r == MSK_RES_OK;
}

bool MOSEK_SolveWrapper::setup(const LP_Constraints_Sparse & cstraints) //cstraints <-> constraints
{
  assert(nbParams_ == cstraints.nbParams_);

  MSK_deletetask(&task);

  int NUMVAR = this->nbParams_;
  int NUMCON = cstraints.constraint_mat_.rows();
  int NUMANZ = cstraints.constraint_mat_.nonZeros();

  MSKrescodee r = MSK_RES_OK;
  if ( r==MSK_RES_OK )
  {
    // Directs the log task stream to the 'printstr' function.
    if ( r==MSK_RES_OK )
      MSK_linkfunctotaskstream(task,MSK_STREAM_LOG,NULL,printstr);

    // Create the optimization task.
    r = MSK_maketask(env,NUMCON,NUMVAR,&task);

    // Give MOSEK an estimate of the size of the input data.
    //This is done to increase the speed of inputting data.
    // However, it is optional.
    if (r == MSK_RES_OK)
      r = MSK_putmaxnumvar(task,NUMVAR);

    if (r == MSK_RES_OK)
      r = MSK_putmaxnumcon(task,NUMCON);

    if (r == MSK_RES_OK)
      r = MSK_putmaxnumanz(task,NUMANZ);

    // Append 'NUMCON' empty constraints. The constraints will initially have no bounds.
    if ( r == MSK_RES_OK )
      r = MSK_append(task,MSK_ACC_CON,NUMCON);

    // Append 'NUMVAR' variables. The variables will initially be fixed at zero (x=0).
    if ( r == MSK_RES_OK )
      r = MSK_append(task,MSK_ACC_VAR,NUMVAR);
  }

  this->nbParams_ = NUMVAR;

  if (cstraints.bminimize_) {
    r = MSK_putobjsense(task, MSK_OBJECTIVE_SENSE_MINIMIZE);
  }
  else  {
    r = MSK_putobjsense(task, MSK_OBJECTIVE_SENSE_MAXIMIZE);
  }

  // Optionally add a constant term to the objective.
  if ( r ==MSK_RES_OK )
    r = MSK_putcfix(task,0.0);

  // Add the objective function if any
  if (!cstraints.vec_cost_.empty())  {
    // Set objective
    for (size_t i = 0; i < cstraints.vec_cost_.size(); ++i)
      r = MSK_putcj(task, i, cstraints.vec_cost_[i]);
  }

  //Add constraint row by row
  const sRMat & A = cstraints.constraint_mat_;
  std::vector<int> vec_colno;
  std::vector<double> vec_value;
  for (int i=0; i<A.rows(); ++i)
  {
    vec_colno.resize(0);
    vec_value.resize(0);
    for (sRMat::InnerIterator it(A,i); it; ++it)
    {
      vec_colno.push_back(it.col());
      vec_value.push_back(it.value());
    }
    // Insert a row
    r = MSK_putavec(task,
                    MSK_ACC_CON,      // Input row of A.
                    MSKidxt(i),       // Variable row index.
                    vec_colno.size(), // Number of non-zeros in row(i) A.
                    &vec_colno[0],    // Pointer to row indexes of row i.
                    &vec_value[0]);   // Pointer to Values of row i.
  }

  //Add bound on variables:
  if (cstraints.vec_bounds_.size() == 1) {

    for( size_t i = 0; i < NUMVAR; ++i) {
      if(r == MSK_RES_OK)
        r = MSK_putbound(task,
                         MSK_ACC_VAR,                         // Put bounds on variables.
                         MSKidxt(i),                          // Index of variable.
                         MSK_BK_RA, // Bound key.
                         cstraints.vec_bounds_[0].first,      // Numerical value of lower bound.
                         cstraints.vec_bounds_[0].second);    // Numerical value of upper bound.
    }
  }
  else{

    for( size_t i = 0; i < NUMVAR; ++i) {
    // Set the bounds on variable j.
    //    lower bound <= x_j <= upper bound

      if(r == MSK_RES_OK)
        r = MSK_putbound(task,
                         MSK_ACC_VAR,                         // Put bounds on variables.
                         MSKidxt(i),                          // Index of variable.
                         MSK_BK_RA, // Bound key.
                         cstraints.vec_bounds_[i].first,      // Numerical value of lower bound.
                         cstraints.vec_bounds_[i].second);    // Numerical value of upper bound.
    // in order to add sparse bounds use: MSK_putboundlist
    }
  }

  // Add bounds on constraint
  for(size_t i=0; i<NUMCON && r==MSK_RES_OK; ++i)
  {
    MSKboundkey_enum cst = convertSign(cstraints.vec_sign_[i]);
    if (cst == MSK_BK_UP || cst == MSK_BK_RA)
      r = MSK_putbound(task,
                     MSK_ACC_CON,   // Put bounds on constraints.
                     MSKidxt(i),    // Index of constraint.
                     cst,           // Bound key.
                     -MSK_INFINITY, // Numerical value of lower bound.
                     cstraints.constraint_objective_(i));      //upper bound.
    else if (cst == MSK_BK_LO)
      r = MSK_putbound(task,
                     MSK_ACC_CON,   // Put bounds on constraints.
                     MSKidxt(i),    // Index of constraint.
                     cst,           // Bound key.
                     cstraints.constraint_objective_(i),
                    +MSK_INFINITY);
    else if (cst == MSK_BK_FX)
      r = MSK_putbound(task,
                     MSK_ACC_CON,   // Put bounds on constraints.
                     MSKidxt(i),    // Index of constraint.
                     cst,           // Bound key.
                     cstraints.constraint_objective_(i),
                     cstraints.constraint_objective_(i));
  }

  return r == MSK_RES_OK;
}

bool MOSEK_SolveWrapper::solve()
{
  MSKrescodee trmcode;

  // Run optimizer
  //MSK_putintparam(task,MSK_IPAR_OPTIMIZER,MSK_OPTIMIZER_DUAL_SIMPLEX);
  MSKrescodee r = MSK_optimizetrm(task,&trmcode);

  // Print a summary containing information
  // about the solution for debugging purposes.
  MSK_solutionsummary (task,MSK_STREAM_LOG);

  if ( r== MSK_RES_OK )
  {
    MSKsolstae solsta;
    MSK_getsolutionstatus (task,
                           MSK_SOL_BAS,
                           NULL,
                           &solsta);
    switch(solsta)
    {
      case MSK_SOL_STA_OPTIMAL:
      case MSK_SOL_STA_NEAR_OPTIMAL:{
        return true;
        break;
      }
      case MSK_SOL_STA_DUAL_INFEAS_CER:
      case MSK_SOL_STA_PRIM_INFEAS_CER:
      case MSK_SOL_STA_NEAR_DUAL_INFEAS_CER:
      case MSK_SOL_STA_NEAR_PRIM_INFEAS_CER:
        //printf("Primal or dual infeasibility certificate found.\n");
        break;
      case MSK_SOL_STA_UNKNOWN:
        //printf("The status of the solution could not be determined.\n");
        break;
      default:
        //printf("Other solution status.");
        break;
    }
  }
  else
  {
    printf("Error while optimizing.\n");
  }
  return false;
}

bool MOSEK_SolveWrapper::getSolution(std::vector<double> & estimatedParams)
{
  bool bRet = false;
  MSKsolstae solsta;
  MSK_getsolutionstatus (task,
                           MSK_SOL_BAS,
                           NULL,
                           &solsta);
  switch(solsta)
  {
    case MSK_SOL_STA_OPTIMAL:
    case MSK_SOL_STA_NEAR_OPTIMAL:
      MSK_getsolutionslice(task,
                           MSK_SOL_BAS,    // Request the basic solution.
                           MSK_SOL_ITEM_XX,// Which part of solution.
                           0,              // Index of first variable.
                           estimatedParams.size(), // Index of last variable+1.
                           &estimatedParams[0]);
      bRet = true;
    //MSK_writedata(task,"taskdump.opf");

      /*printf("Optimal primal solution\n");
      for(size_t j=0; j<estimatedParams.size(); ++j)
        std::cout << estimatedParams[j] << " ";
      std::cout << std::endl;*/

      break;
    case MSK_SOL_STA_DUAL_INFEAS_CER:
    case MSK_SOL_STA_PRIM_INFEAS_CER:
    case MSK_SOL_STA_NEAR_DUAL_INFEAS_CER:
    case MSK_SOL_STA_NEAR_PRIM_INFEAS_CER:
      //printf("Primal or dual infeasibility certificate found.\n");
      break;

    case MSK_SOL_STA_UNKNOWN:
      //printf("The status of the solution could not be determined.\n");
      break;
    default:
      //printf("Other solution status.");
      break;
  }
  return bRet;
}

} //namespace linearProgramming
} //namespace openMVG

#endif // OPENMVG_HAVE_MOSEK
