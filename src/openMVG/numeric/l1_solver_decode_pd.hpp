// Copyright (c) 2014 cDc and Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_NUMERIC_L1_SOLVER_DECODE_PD_H_
#define OPENMVG_NUMERIC_L1_SOLVER_DECODE_PD_H_

#include <Eigen/Dense>

namespace openMVG {

// Minimum l1 error approximation:
//
// Let A be a M x N matrix with full rank. Given y of R^M, the problem
// (PA) minx ||y - Ax||1
// finds the vector x of R^N such that the error y - Ax has minimum l1 norm
// (i.e. we are asking that the difference between Ax and y be sparse).
// This problem arises in the context of channel coding
// (see E. J. Candes and T. Tao. "Decoding by linear programming". in IEEE Trans.
// Inform. Theory, December 2005).
//
// Suppose we have a channel code that produces a codeword c = Ax for a message x. The
// message travels over the channel, and has an unknown number of its entries corrupted.
// The decoder observes y = c + e, where e is the corruption. If e is sparse enough, then
// the decoder can use (PA) to recover x exactly. When x, A, y have real-valued entries,
// (PA) can be recast as an LP.
//
typedef double REAL;

template<typename MATRIX_TYPE>
inline bool TRobustRegressionL1PD(
  const MATRIX_TYPE& A,
  const Eigen::Matrix<REAL, Eigen::Dynamic, 1>& y,
  Eigen::Matrix<REAL, Eigen::Dynamic, 1>& xp,
  REAL pdtol=1e-3, unsigned pdmaxiter=50)
{
  typedef Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;
  typedef Eigen::Matrix<REAL, Eigen::Dynamic, 1> Vector;
  const unsigned M = (unsigned)y.size();
  const unsigned N = (unsigned)xp.size();
  assert(A.rows() == M && A.cols() == N);

  const REAL alpha(0.01);
  const REAL beta(0.5);
  const REAL mu(10);

  Vector x(xp);
  Vector Ax(A*x);
  Vector tmpM1(y-Ax);
  Vector tmpM2(-tmpM1);
  Vector tmpM3(tmpM1.cwiseAbs()), tmpM4(M);
  Vector u = (tmpM3*REAL(0.95)).array() + tmpM3.maxCoeff()*REAL(0.10);
  Vector fu1 = tmpM2-u;
  Vector fu2 = tmpM1-u;

  Vector lamu1(M), lamu2(M);
  for (unsigned i=0; i<M; ++i) {
    lamu1(i) = -1.0/fu1(i);
    lamu2(i) = -1.0/fu2(i);
  }
  const MATRIX_TYPE At(A.transpose());
  Vector Atv(At*(lamu1-lamu2));
  REAL AtvNormSq = Atv.squaredNorm();
  Vector rdual((-lamu1-lamu2).array() + REAL(1));
  REAL rdualNormSq = rdual.squaredNorm();

  Vector w2(M), sig1(M), sig2(M), sigx(M), dx(N), up(N), Atdv(N);
  Vector Axp(M), Atvp(M);
  Vector &Adx(sigx), &du(w2), &w1p(dx);
  Matrix H11p(N,N);
  Vector &dlamu1(tmpM3), &dlamu2(tmpM4);
  for (unsigned pditer=0; pditer<pdmaxiter; ++pditer) {
    // surrogate duality gap
    const REAL sdg(-(fu1.dot(lamu1) + fu2.dot(lamu2)));
    if (sdg < pdtol)
      break;
    const REAL tau(mu*2*M/sdg);
    const REAL inv_tau = REAL(-1)/tau;
    tmpM1 = (-lamu1.cwiseProduct(fu1)).array() + inv_tau;
    tmpM2 = (-lamu2.cwiseProduct(fu2)).array() + inv_tau;
    const REAL resnorm = sqrt(AtvNormSq + rdualNormSq + tmpM1.squaredNorm() + tmpM2.squaredNorm());

    for (unsigned i=0; i<M; ++i) {
      REAL& tmpM3i = tmpM3(i);
      tmpM3i = inv_tau/fu1(i);
      REAL& tmpM4i = tmpM4(i);
      tmpM4i = inv_tau/fu2(i);
      w2(i) = tmpM3i + tmpM4i - REAL(1);
    }

    tmpM1 = lamu1.cwiseQuotient(fu1);
    tmpM2 = lamu2.cwiseQuotient(fu2);
    sig1 = -tmpM1 - tmpM2;
    sig2 = tmpM1 - tmpM2;
    sigx = sig1 - sig2.cwiseAbs2().cwiseQuotient(sig1);

    H11p = At*(Eigen::DiagonalMatrix<REAL,Eigen::Dynamic>(sigx)*A);
    w1p = At*(tmpM4 - tmpM3 - (sig2.cwiseQuotient(sig1).cwiseProduct(w2)));

    // optimized solver as A is positive definite and symmetric
    dx = H11p.ldlt().solve(w1p);

    Adx = A*dx;

    du = (w2 - sig2.cwiseProduct(Adx)).cwiseQuotient(sig1);

    dlamu1 = -tmpM1.cwiseProduct(Adx-du) - lamu1 + tmpM3;
    dlamu2 =  tmpM2.cwiseProduct(Adx+du) - lamu2 + tmpM4;
    Atdv = At*(dlamu1-dlamu2);

    // make sure that the step is feasible: keeps lamu1,lamu2 > 0, fu1,fu2 < 0
    REAL s(1);
    for (unsigned i=0; i<M; ++i) {
      REAL& dlamu1i = dlamu1(i);
      if (dlamu1i < 0) {
        const REAL tmp = -lamu1(i)/dlamu1i;
        if (s > tmp)
          s = tmp;
      }
      REAL& dlamu2i = dlamu2(i);
      if (dlamu2i < 0) {
        const REAL tmp = -lamu2(i)/dlamu2i;
        if (s > tmp)
          s = tmp;
      }
    }
    for (unsigned i=0; i<M; ++i) {
      REAL& Adxi = Adx(i);
      REAL& dui = du(i);
      REAL Adx_du = Adxi-dui;
      if (Adx_du > 0) {
        const REAL tmp = -fu1(i)/Adx_du;
        if (s > tmp)
          s = tmp;
      }
      Adx_du = -Adxi-dui;
      if (Adx_du > 0) {
        const REAL tmp = -fu2(i)/Adx_du;
        if (s > tmp)
          s = tmp;
      }
    }
    s *= REAL(0.99);

    // backtrack
    lamu1 += s*dlamu1;  lamu2 += s*dlamu2;
    rdual = (-lamu1-lamu2).array() + REAL(1);
    rdualNormSq = rdual.squaredNorm();
    bool suffdec = false;
    unsigned backiter = 0;
    do {
      xp = x + s*dx;  up = u + s*du;
      Axp = Ax + s*Adx;  Atvp = Atv + s*Atdv;
      fu1 = Axp - y - up;  fu2 = -Axp + y - up;
      AtvNormSq = Atvp.squaredNorm();
      tmpM1 = (-lamu1.cwiseProduct(fu1)).array() + inv_tau;
      tmpM2 = (-lamu2.cwiseProduct(fu2)).array() + inv_tau;
      const REAL newresnorm = sqrt(AtvNormSq + rdualNormSq + tmpM1.squaredNorm() + tmpM2.squaredNorm());
      suffdec = (newresnorm <= (REAL(1)-alpha*s)*resnorm);
      s = beta*s;
      if (++backiter > 32) {
        //("error: stuck backtracking, returning last iterate"); // see Section 4 of notes for more information
        xp.swap(x);
        return false;
      }
    } while (!suffdec);

    // next iteration
    x.swap(xp);  u.swap(up);
    Ax.swap(Axp);  Atv.swap(Atvp);
  }
  return true;
}

}  // namespace openMVG

#endif  // OPENMVG_NUMERIC_L1_SOLVER_DECODE_PD_H_

