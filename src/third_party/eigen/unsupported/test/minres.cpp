// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2011 Gael Guennebaud <g.gael@free.fr>
// Copyright (C) 2012 Giacomo Po <gpo@ucla.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
#include <cmath>

#include "../../test/sparse_solver.h"
#include <Eigen/IterativeSolvers>

template<typename T> void test_minres_T()
{
  MINRES<SparseMatrix<T>, Lower, DiagonalPreconditioner<T> > minres_colmajor_diag;
  MINRES<SparseMatrix<T>, Lower, IdentityPreconditioner    > minres_colmajor_I;
//  MINRES<SparseMatrix<T>, Lower, IncompleteLUT<T> >           minres_colmajor_ilut;
  //minres<SparseMatrix<T>, SSORPreconditioner<T> >     minres_colmajor_ssor;

  CALL_SUBTEST( check_sparse_square_solving(minres_colmajor_diag)  );
  CALL_SUBTEST( check_sparse_spd_solving(minres_colmajor_I) );
 // CALL_SUBTEST( check_sparse_square_solving(minres_colmajor_ilut)     );
  //CALL_SUBTEST( check_sparse_square_solving(minres_colmajor_ssor)     );
}

void test_minres()
{
  CALL_SUBTEST_1(test_minres_T<double>());
//  CALL_SUBTEST_2(test_minres_T<std::complex<double> >());
}
