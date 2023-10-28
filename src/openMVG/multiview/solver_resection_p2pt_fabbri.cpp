// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2020 Ricardo Fabbri, Ariel Kovaljski
//
// Author: Ariel Kovaljski and Ricardo Fabbri
// Rio de Janeiro State University

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
//
// This works OK with --ffast_math
//
// TODO(optimization)
//  - variabele alignment and order
//    - specially for fn_t and pose_from_point_tangents_2
//  - make pose_poly member to p2pt class, eliminate pose poly class
//  - eliminate pointers to root_ids
//  - root_ids also as member of pose_poly
//  - when we mark as 1, store index and only redo that part.
//  - many times we have ints multiplying doubles
//

#include <array>
#include <cmath>
#include <iomanip>

#include "openMVG/multiview/solver_resection_p2pt_fabbri.hpp"
#include "openMVG/multiview/projection.hpp"


// Visual studio bug work-around
#ifndef M_PI_2
    #define M_PI_2 1.57079632679489661923
#endif


#if defined(_MSC_VER)
#define __attribute__(x) /* blank - should simply ignore thanks to C preprocessor */
#endif

namespace openMVG
{
namespace euclidean_resection
{
  
// At most 8 solutions with positive depth, 16 total
static constexpr unsigned char TS_MAX_LEN = 32;
static constexpr unsigned char RT_MAX_LEN = 4*TS_MAX_LEN;

template <typename T=double>
class p2pt { // fully static, not to be instantiated - just used for templating
	public:
  static bool pose_from_point_tangents(
    const T gama1[3], const T tgt1[3],
    const T gama2[3], const T tgt2[3],
    const T Gama1[3], const T Tgt1[3],
    const T Gama2[3], const T Tgt2[3],
    T (*output_RT)[RT_MAX_LEN][4][3],
    unsigned char *output_RT_len,
    T *output_degen
  );
};

// polynomial in t
// solution is pose(t)
template<typename T>
struct pose_poly {
  // t-independent terms
	T A0, A1, A2, B0, B1, B2, B3, C0, C1, C2, C3, C4,
		E0, E1, E2, F0, F1, F2, F3, G0, G1, G2, G3, G4,
		H0, H1, H2, H3, H4, J0, J1, J2, J3, K0, K1, K2, K3,
		L0, L1, L2, alpha, beta, theta, sth, cth;

  static constexpr unsigned short T_LEN = 2001 /* fine-tune based on data, but has to be about >= 1500 */, 
                                  ROOT_IDS_LEN = T_LEN - 1;
  static constexpr double T_LEN_2 = 2./(T_LEN-1);
  inline T t_vec(unsigned short i) { return T_LEN_2*i -1.; }

	void pose_from_point_tangents_2(
		const T gama1[3], const T tgt1[3], const T gama2[3], const T tgt2[3],
		const T Gama1[3], const T Tgt1[3], const T Gama2[3], const T Tgt2[3]);
  
	inline T __attribute__((always_inline)) fn_t(const T t, T p[10]) { // function of t part
    const T t2 = t*t, t3 = t2*t, t4 = t3*t, t5 = t4*t, t6 = t5*t, t7 = t6*t, t8 = t7*t;
    const T t2p12 = (t2 + 1.) * (t2 + 1.), t2p13 = t2p12 * (t2 + 1.), t2p14 = t2p13 * (t2 + 1.);

    T &A = p[0]; T &B = p[1]; T &C = p[2]; T &E = p[3]; T &F = p[4]; 
    T &G = p[5]; T &H = p[6]; T &J = p[7]; T &K = p[8]; T &L = p[9];
    A = (A0+A1*t+A2*t2-A1*t3+A0*t4)/t2p12;
    B = (B0+B1*t+B2*t2+B3*t3-B2*t4+B1*t5-B0*t6)/t2p13;
    C = (C0+C1*t+C2*t2+C3*t3+C4*t4-C3*t5+C2*t6-C1*t7+C0*t8)/t2p14;
    E = (E0+E1*t+E2*t2-E1*t3+E0*t4)/t2p12;
    F = (F0+F1*t+F2*t2+F3*t3-F2*t4+F1*t5-F0*t6)/t2p13;
    G = (G0+G1*t+G2*t2+G3*t3+G4*t4-G3*t5+G2*t6-G1*t7+G0*t8)/t2p14;
    H = (H0+H1*t+H2*t2+H3*t3+H4*t4-H3*t5+H2*t6-H1*t7+H0*t8)/t2p14;
    J = (J0+J1*t+J2*t2+J3*t3-J2*t4+J1*t5-J0*t6)/t2p13;
    K = (K0+K1*t+K2*t2+K3*t3-K2*t4+K1*t5-K0*t6)/t2p13;
    L = (L0+L1*t+L2*t2-L1*t3+L0*t4)/t2p12;

    const T AA=A*A, BB=B*B, CC=C*C, EE=E*E, FF=F*F, GG=G*G, HH=H*H, 
            H3=HH*H, H4=H3*H, JJ=J*J, J3=JJ*J, KK=K*K, K3=KK*K, LL=L*L, L3=LL*L;
    
    return EE*BB*HH*JJ +GG*CC*L3*L +GG*AA*K3*K +EE*AA*H4 +EE*CC*J3*J
    -2.*E*A*HH*G*C*LL +2.*EE*A*HH*C*JJ -2.*EE*C*J3*B*H +2.*E*CC*JJ*G*LL
    +2.*E*AA*HH*G*KK -2.*EE*A*H3*B*J -2.*E*A*HH*G*B*K*L -2.*E*C*JJ*G*B*K*L
    -2.*E*C*JJ*G*A*KK -2.*E*B*H*J*G*C*LL -2.*E*B*H*J*G*A*KK +GG*BB*KK*LL
    -2.*GG*B*K*L3*C -2.*GG*B*K3*L*A +2.*GG*C*LL*A*KK -2.*F*E*AA*H3*K
    -2.*F*E*A*H*K*C*JJ +3.*F*E*A*HH*K*B*J +3.*F*A*H*KK*G*B*L
    -2.*F*A*H*K*G*C*LL -2.*F*AA*H*K3*G +F*E*B*H3*L*A +3.*F*E*B*H*L*C*JJ
    -F*E*BB*HH*L*J -F*BB*H*LL*G*K +F*B*H*L3*G*C +F*E*B*K*J3*C
    -F*E*BB*K*JJ*H -F*BB*KK*J*G*L +3.*F*B*K*J*G*C*LL +F*B*K3*J*G*A
    -2.*F*E*C*J*L*A*HH -2.*F*E*CC*J3*L -2.*F*CC*J*L3*G -2.*F*C*J*L*G*A*KK
    +FF*AA*KK*HH +FF*A*KK*C*JJ -FF*A*KK*B*H*J -FF*B*K*L*A*HH
    -FF*B*K*L*C*JJ +FF*BB*K*L*H*J +FF*C*LL*A*HH +FF*CC*LL*JJ
    -FF*C*LL*B*H*J +G*E*BB*HH*LL +G*E*BB*KK*JJ +8.*G*E*A*H*K*C*J*L;
  }

	inline void find_bounded_root_intervals(bool (*root_ids_out)[ROOT_IDS_LEN]) {
	  T p[10];
    T curr_val = fn_t(t_vec(0), p), next_val;
    // std::cout << "fn_t [";
    for (unsigned short i = 0; i < ROOT_IDS_LEN; i++) {
      next_val = fn_t(t_vec(i+1), p);
      static constexpr T eps = std::numeric_limits<T>::epsilon();
      (*root_ids_out)[i] = curr_val > +eps && next_val < -eps || 
                           curr_val < -eps && next_val > +eps;
      // std::cout << curr_val << " ";
      curr_val = next_val;
    }
    // std::cout << "]\n";
  }
  
	// inline T operator()(T t) { return fn_t(t); }
  
	inline void rhos_from_root_ids(const bool (&root_ids)[ROOT_IDS_LEN], 
                                 T (*out)[3][TS_MAX_LEN], 
                                 unsigned char *out_ts_len) {

    T (&ts)[TS_MAX_LEN] = (*out)[0];
    T (&rhos1)[TS_MAX_LEN] = (*out)[1]; T (&rhos2)[TS_MAX_LEN] = (*out)[2];
    T p[10];
    unsigned char &ts_end = *out_ts_len; ts_end = 0;
    for (unsigned short i = 0; i < ROOT_IDS_LEN; i++) {
      if (!root_ids[i]) continue;
      T t0 = t_vec(i), t1 = t_vec(i+1), t;
      T f0 = fn_t(t_vec(i), p), f1 = fn_t(t_vec(i+1), p);
      for (unsigned char k = 0; k < 4; ++k) {
        t = t1 - f1*(t1-t0)/(f1-f0); t0 = t1; t1 = t;
        f0 = f1; if (k + 1 < 4) f1 = fn_t(t, p);
      }
      // Root is t, plus minus t_stddev. Now get rho1(t):

      const T tt = t*t, alpha_times_2 = 2.*alpha,
      alpha_ts_new2 = alpha_times_2 * t, beta_1_minus_tt = beta * (1. - tt);
      const T r1 =  alpha_ts_new2 * cth + beta_1_minus_tt * sth; if (r1 <= 1e-12) continue;
      const T r2 = -alpha_ts_new2 * sth + beta_1_minus_tt * cth; if (r2 <= 1e-12) continue;
      const T ts_den = 1. + tt; 
      rhos1[ts_end] = r1 / ts_den; rhos2[ts_end] = r2 / ts_den; ts[ts_end++] = t;
      assert(ts_end <= TS_MAX_LEN);
    }
    // std::cout << "ts_end: " << ts_end << std::endl;
  }
  
	void get_sigmas(const unsigned char ts_len, const T (&ts)[TS_MAX_LEN], 
      T (*out)[2][TS_MAX_LEN][4], unsigned char out_len[TS_MAX_LEN]);
  
	void get_r_t_from_rhos(
		const unsigned char ts_len,
		const T sigmas1[TS_MAX_LEN][4], const unsigned char sigmas_len[TS_MAX_LEN],
		const T sigmas2[TS_MAX_LEN][4],
		const T rhos1[TS_MAX_LEN], const T rhos2[TS_MAX_LEN],
		const T gama1[3], const T tgt1[3], const T gama2[3], const T tgt2[3],
		const T Gama1[3], const T Tgt1[3], const T Gama2[3], const T Tgt2[3],
		T (*out)[RT_MAX_LEN][4][3], unsigned char *out_len
	);
};
  
// This is the main routine ----------------------------------------------------
template <typename T>
bool p2pt<T>::
pose_from_point_tangents(
	const T gama1[3], const T tgt1[3], const T gama2[3], const T tgt2[3],
	const T Gama1[3], const T Tgt1[3], const T Gama2[3], const T Tgt2[3],
	T (*output_RT)[RT_MAX_LEN][4][3], unsigned char *output_RT_len, T *output_degen
)
{
  { // test for geometric degeneracy -------------------------------
    T DGama[3] = { Gama1[0] - Gama2[0], Gama1[1] - Gama2[1], Gama1[2] - Gama2[2] };
    const T norm = sqrt(DGama[0]*DGama[0] + DGama[1]*DGama[1] + DGama[2]*DGama[2]);
    const T d[3][3] = { // Matrix for degeneracy calculation
      DGama[0]/norm, Tgt1[0], Tgt2[0],
      DGama[1]/norm, Tgt1[1], Tgt2[1],
      DGama[2]/norm, Tgt1[2], Tgt2[2]
    };
    T &degen = *output_degen;
    degen = (d[0][0]*d[1][1]*d[2][2]+d[0][1]*d[1][2]*d[2][0]+d[0][2]*d[1][0]*d[2][1]) // det(d)
           -(d[2][0]*d[1][1]*d[0][2]+d[2][1]*d[1][2]*d[0][0]+d[2][2]*d[1][0]*d[0][1]);

    // std::cout << "degen: " << degen << std::endl;
    if (std::fabs(degen) < 0.09) {
      *output_RT_len = 0;
      return false;  // can still solve this in many cases, but lets not fool around
    }
  }

	// compute roots -------------------------------
	pose_poly<T> p;
	p.pose_from_point_tangents_2(gama1, tgt1, gama2, tgt2, Gama1, Tgt1, Gama2, Tgt2);

	bool root_ids[pose_poly<T>::ROOT_IDS_LEN];
	p.find_bounded_root_intervals(&root_ids);

	// compute rhos, r, t --------------------------
	T rhos[3][TS_MAX_LEN];
	unsigned char ts_len;
	p.rhos_from_root_ids(root_ids, &rhos, &ts_len);

	const T (&ts)[TS_MAX_LEN]    = rhos[0]; 
  const T (&rhos1)[TS_MAX_LEN] = rhos[1]; 
  const T (&rhos2)[TS_MAX_LEN] = rhos[2];
	T sigmas[2][TS_MAX_LEN][4]; unsigned char sigmas_len[TS_MAX_LEN];

 	p.get_sigmas(ts_len, ts, &sigmas, sigmas_len);

	const T (&sigmas1)[TS_MAX_LEN][4] = sigmas[0];
	const T (&sigmas2)[TS_MAX_LEN][4] = sigmas[1];

	T (&RT)[RT_MAX_LEN][4][3] = *output_RT;
	unsigned char &RT_len     = *output_RT_len;

	p.get_r_t_from_rhos(ts_len, sigmas1, sigmas_len, sigmas2,
		rhos1, rhos2, gama1, tgt1, gama2, tgt2, Gama1, Tgt1, Gama2, Tgt2, 
    &RT, &RT_len);
  return true;
}

// From problem input, build poly coefficients (independent of t).
// Behind the scenes, there is a polynomial in t
// The polynomial itself can be valuated using fn_t()
template<typename T>
void pose_poly<T>::
pose_from_point_tangents_2(
	const T gama1[3], const T tgt1[3],
	const T gama2[3], const T tgt2[3],
	const T Gama1[3], const T Tgt1[3],
	const T Gama2[3], const T Tgt2[3]
)
{
	const T g11 = gama1[0], g12 = gama1[1], g21 = gama2[0], g22 = gama2[1],
          g11_2 = g11*g11, g11_3 = g11_2*g11, g11_4 = g11_3*g11,
          g12_2 = g12*g12, g12_3 = g12_2*g12, g12_4 = g12_3*g12,
          g21_2 = g21*g21, g21_3 = g21_2*g21, g21_4 = g21_3*g21,
          g22_2 = g22*g22, g22_3 = g22_2*g22, g22_4 = g22_3*g22,
          h11 = tgt1[0],  h12 = tgt1[1], h21 = tgt2[0],  h22 = tgt2[1];

	T *V = &A0; // reusing memory from poly
  V[0] = Gama1[0] - Gama2[0]; V[1] = Gama1[1] - Gama2[1]; V[2] = Gama1[2] - Gama2[2];
  
	const T 
  a1 = V[0]*V[0]+V[1]*V[1]+V[2]*V[2],
  a2 = Tgt1[0]*Tgt1[0]+Tgt1[1]*Tgt1[1]+Tgt1[2]*Tgt1[2],
  a3 = Tgt2[0]*Tgt2[0]+Tgt2[1]*Tgt2[1]+Tgt2[2]*Tgt2[2],
  a4 = V[0]*Tgt1[0]+V[1]*Tgt1[1]+V[2]*Tgt1[2],
  a5 = Tgt1[0]*Tgt2[0]+Tgt1[1]*Tgt2[1]+Tgt1[2]*Tgt2[2],
  a6 = V[0]*Tgt2[0]+V[1]*Tgt2[1]+V[2]*Tgt2[2];

	theta = 0.5 * atan( 2.*(1.+g11*g21+g12*g22)/(g11_2+g12_2-g21_2-g22_2) );
	if (theta < 0) theta += M_PI_2;
	sth = sin(theta); const T s2 = sth*sth; cth = cos(theta);
  const T c2 = cth*cth, c2th = 2.*c2-1., s2th = 2.*sth*cth;
  
	const T den1 = 2.*s2th*(g11*g21+g12*g22+1.)+c2th*(g11_2+g12_2-g21_2-g22_2),
	        den2 = g11_2 + g12_2 + g21_2 + g22_2 + 2.;
	beta  = sqrt(-2.*a1 / (den1 - den2)); // sqrt(t25)
	alpha = sqrt(2.*a1 / (den1 + den2));  // sqrt(t24)

	// Coefficient code adapted from Maple ::: can be further cleaned up but works
  // TODO: remove parenthesis. Further reduce to known gates
  // Perhaps use gate optimization.
	A0 = a4*a4*g12_2
	+a4*a4*g11_2
	+a4*a4
	+2.0*a2*g11_3*g21*beta*beta*sth*cth
	+2.0*a2*g21*g11*g12_2*beta*beta*sth*cth
	-2.0*a2*g11_2*g12_2*beta*beta*s2
	-a2*g12_4*beta*beta*s2
	-a2*g21_2*g11_2*beta*beta*c2
	+2.0*a2*g12_2*beta*beta*sth*cth
	+2.0*a2*g11_2*beta*beta*sth*cth
	+2.0*a2*g11_2*g22*g12*beta*beta*sth*cth
	-a2*beta*beta*c2
	+2.0*a2*g12_3*g22*beta*beta*sth*cth
	-a2*g11_4*beta*beta*s2
	-2.0*a2*g11_2*beta*beta*s2
	-2.0*a2*g12_2*beta*beta*s2
	+2.0*a2*beta*beta*sth*cth
-2.0*a2*g21*g11*g22*g12*beta*beta*c2
	-a2*beta*beta*s2
	+2.0*a2*g21*g11*beta*beta*sth*cth
	-a2*g22_2*g12_2*beta*beta*c2
	-2.0*a2*g22*g12*beta*beta*c2
	-2.0*a2*g21*g11*beta*beta*c2
	+2.0*a2*g22*g12*beta*beta*sth*cth;

	A1 = 4.*a2*alpha*c2*beta
	-4.*a2*beta*s2*alpha
	+4.*a2*g21_2*g11_2*alpha*sth*beta*cth
	+4.*a2*g12_2*alpha*c2*beta
	+8.*a2*g21*g11*alpha*sth*beta*cth
	-4.*a2*g12_2*beta*s2*alpha
	+4.*a2*g22*g12*alpha*c2*beta
	-4.*a2*g22*g12*beta*s2*alpha
	+4.*a2*g22_2*g12_2*alpha*sth*beta*cth
	-4.*a2*g21*g11*g12_2*beta*s2*alpha
	+4.*a2*g21*g11*alpha*c2*beta
	+4.*a2*g11_2*alpha*c2*beta
	-4.*a2*g11_2*beta*s2*alpha
	+4.*a2*g21*g11*g12_2*alpha*c2*beta
	-4.*a2*g21*g11*beta*s2*alpha
	-8.*a2*g11_2*g12_2*alpha*sth*beta*cth
	-4.*a2*g11_4*alpha*sth*beta*cth
	-8.*a2*g11_2*alpha*sth*beta*cth
	+8.*a2*g21*g11*g22*g12*alpha*sth*beta*cth
	+4.*a2*g12_3*g22*alpha*c2*beta
	-4.*a2*g12_3*g22*beta*s2*alpha
	-4.*a2*g12_4*alpha*sth*beta*cth
	-8.*a2*g12_2*alpha*sth*beta*cth
	+8.*a2*g22*g12*alpha*sth*beta*cth
	+4.*a2*g11_3*g21*alpha*c2*beta
	-4.*a2*g11_3*g21*beta*s2*alpha
	+4.*a2*g11_2*g22*g12*alpha*c2*beta
	-4.*a2*g11_2*g22*g12*beta*s2*alpha;

	A2 = (2*a4*a4*g12_2)
	+(2*a4*a4*g11_2)
	+(2*a4*a4)
	+2.*a2*g12_4*beta*beta*s2
	+2.*a2*g11_4*beta*beta*s2
	+4.*a2*(g11_2)*beta*beta*s2
	+4.*a2*(g12_2)*beta*beta*s2
	-4.*a2*beta*beta*sth*cth
	+2.*a2*beta*beta*s2
	+2.*a2*beta*beta*c2
	-4.*a2*g21*g11*(g12_2)*beta*beta*sth*cth
	+2.*a2*g21_2*(g11_2)*beta*beta*c2
	-4.*a2*(g12_2)*beta*beta*sth*cth
	+4.*a2*g21*g11*beta*beta*c2
	+2.*a2*g22_2*(g12_2)*beta*beta*c2
	-4.*a2*g22*g12*beta*beta*sth*cth
	-4.*a2*(g11_2)*beta*beta*sth*cth
	-4.*a2*g21*g11*beta*beta*sth*cth
	+4.*a2*(g11_2)*(g12_2)*beta*beta*s2
	+4.*a2*g22*g12*beta*beta*c2
	-4.*a2*g11_3*g21*beta*beta*sth*cth
	-4.*a2*(g11_2)*g22*g12*beta*beta*sth*cth
	+4.*a2*g21*g11*g22*g12*beta*beta*c2
	-4.*a2*g12_3*g22*beta*beta*sth*cth
	-4.*a2*g11_4*alpha*alpha*c2
	-8.*a2*(g11_2)*alpha*alpha*c2
	-4.*a2*g12_4*alpha*alpha*c2
	-8.*a2*(g12_2)*alpha*alpha*c2
	-8.*a2*alpha*alpha*cth*sth
	-4.*a2*alpha*alpha*c2
	-4.*a2*alpha*alpha*s2
	-8.*a2*g22*g12*alpha*alpha*cth*sth
	-4.*a2*g21_2*(g11_2)*alpha*alpha*s2
	-8.*a2*(g12_2)*alpha*alpha*cth*sth
	-8.*a2*g21*g11*alpha*alpha*s2
	-4.*a2*g22_2*(g12_2)*alpha*alpha*s2
	-8.*a2*g21*g11*alpha*alpha*cth*sth
	-8.*a2*(g11_2)*(g12_2)*alpha*alpha*c2
	-8.*a2*(g11_2)*alpha*alpha*cth*sth
	-8.*a2*g21*g11*(g12_2)*alpha*alpha*cth*sth
	-8.*a2*g21*g11*g22*g12*alpha*alpha*s2
	-8.*a2*g12_3*g22*alpha*alpha*cth*sth
	-8.*a2*g22*g12*alpha*alpha*s2
	-8.*a2*g11_3*g21*alpha*alpha*cth*sth
	-8.*a2*(g11_2)*g22*g12*alpha*alpha*cth*sth;

	B0 = -2.*beta*sth*(a2*g21*g11*g22*h12*beta*beta*c2
	+a2*g12_3*h12*beta*beta*s2
	+a2*g21*h11*g22*g12*beta*beta*c2
	+a2*g11*h11*g12_2*beta*beta*s2
	-a2*g11*h11*beta*beta*sth*cth
	-a4*a4*h11*g11
	+a2*g11*h11*beta*beta*s2
	+a2*g22_2*h12*g12*beta*beta*c2
	+a2*g22*h12*beta*beta*c2
	+a2*g12*h12*beta*beta*s2
	-a2*g11*h11*g22*g12*beta*beta*sth*cth
	-a2*g12*h12*beta*beta*sth*cth
	-a2*g21*h11*g12_2*beta*beta*sth*cth
	-a2*g21*h11*beta*beta*sth*cth
	-a2*g11_2*g22*h12*beta*beta*sth*cth
	-2.*a2*g12_2*h12*g22*beta*beta*sth*cth
	-2.*a2*g11_2*h11*g21*beta*beta*sth*cth
	-a2*g22*h12*beta*beta*sth*cth
	+a2*g11_3*h11*beta*beta*s2
	+a2*g21_2*h11*g11*beta*beta*c2
	+a2*g21*h11*beta*beta*c2
	-a2*g21*g11*g12*h12*beta*beta*sth*cth
	+a2*g11_2*g12*h12*beta*beta*s2
	-a4*a4*h12*g12);

	B1 = -2.*beta*sth*(2.*a2*g11_2*g22*h12*beta*s2*alpha
	-2.*a2*g11_2*g22*h12*alpha*c2*beta
	-2.*a2*g21*h11*g12_2*alpha*c2*beta
	+2.*a2*g21*g11*g12*h12*beta*s2*alpha
	-4.*a2*g12_2*h12*g22*alpha*c2*beta
	+4.*a2*g12_2*h12*g22*beta*s2*alpha
	+2.*a2*g21*h11*g12_2*beta*s2*alpha
	-2.*a2*g11*h11*g22*g12*alpha*c2*beta
	+2.*a2*g11*h11*g22*g12*beta*s2*alpha
	+4.*a2*g11*h11*g12_2*alpha*sth*beta*cth
	-2.*a2*g22*h12*alpha*c2*beta
	+2.*a2*g22*h12*beta*s2*alpha
	-4.*a2*g22_2*h12*g12*alpha*sth*beta*cth
	-2.*a2*g12*h12*alpha*c2*beta
	+2.*a2*g12*h12*beta*s2*alpha
	-4.*a2*g22*h12*alpha*sth*beta*cth
	-2.*a2*g11*h11*alpha*c2*beta
	-4.*a2*g11_2*h11*g21*alpha*c2*beta
	+4.*a2*g11_2*h11*g21*beta*s2*alpha
	-4.*a2*g21*h11*g22*g12*alpha*sth*beta*cth
	+4.*a2*g12*h12*alpha*sth*beta*cth
	+4.*a2*g12_3*h12*alpha*sth*beta*cth
	+4.*a2*g11_3*h11*alpha*sth*beta*cth
	-2.*a2*g21*h11*alpha*c2*beta
	+2.*a2*g21*h11*beta*s2*alpha
	-4.*a2*g21*h11*alpha*sth*beta*cth
	-2.*a2*g21*g11*g12*h12*alpha*c2*beta
	+4.*a2*g11_2*g12*h12*alpha*sth*beta*cth
	+4.*a2*g11*h11*alpha*sth*beta*cth
	-4.*a2*g21*g11*g22*h12*alpha*sth*beta*cth
	+2.*a2*g11*h11*beta*s2*alpha
	-4.*a2*g21_2*h11*g11*alpha*sth*beta*cth)
	-4.*alpha*cth*(a2*g21*g11*g22*h12*beta*beta*c2
	+a2*g12_3*h12*beta*beta*s2
	+a2*g21*h11*g22*g12*beta*beta*c2
	+a2*g11*h11*g12_2*beta*beta*s2
	-a2*g11*h11*beta*beta*sth*cth
	-a4*a4*h11*g11
	+a2*g11*h11*beta*beta*s2
	+a2*g22_2*h12*g12*beta*beta*c2
	+a2*g22*h12*beta*beta*c2
	+a2*g12*h12*beta*beta*s2
	-a2*g11*h11*g22*g12*beta*beta*sth*cth
	-a2*g12*h12*beta*beta*sth*cth
	-a2*g21*h11*g12_2*beta*beta*sth*cth
	-a2*g21*h11*beta*beta*sth*cth
	-a2*g11_2*g22*h12*beta*beta*sth*cth
	-2.*a2*g12_2*h12*g22*beta*beta*sth*cth
	-2.*a2*g11_2*h11*g21*beta*beta*sth*cth
	-a2*g22*h12*beta*beta*sth*cth
	+a2*g11_3*h11*beta*beta*s2
	+a2*g21_2*h11*g11*beta*beta*c2
	+a2*g21*h11*beta*beta*c2
	-a2*g21*g11*g12*h12*beta*beta*sth*cth
	+a2*g11_2*g12*h12*beta*beta*s2
	-a4*a4*h12*g12);

	B2 = -2.*beta*sth*(4.*a2*g21_2*h11*g11*alpha*alpha*s2
	+4.*a2*g11_2*g12*h12*alpha*alpha*c2
	+4.*a2*g11*h11*alpha*alpha*c2
	+4.*a2*g11_2*g22*h12*alpha*alpha*cth*sth
	+4.*a2*g11*h11*alpha*alpha*cth*sth
	+8.*a2*g11_2*h11*g21*alpha*alpha*cth*sth
	+4.*a2*g12_3*h12*alpha*alpha*c2
	+4.*a2*g12*h12*alpha*alpha*c2
	+4.*a2*g11_3*h11*alpha*alpha*c2
	+4.*a2*g21*h11*g22*g12*alpha*alpha*s2
	+4.*a2*g21*h11*alpha*alpha*cth*sth
	+4.*a2*g11*h11*g22*g12*alpha*alpha*cth*sth
	+4.*a2*g11*h11*g12_2*alpha*alpha*c2
	+4.*a2*g21*h11*g12_2*alpha*alpha*cth*sth
	+4.*a2*g22*h12*alpha*alpha*cth*sth
	+4.*a2*g22_2*h12*g12*alpha*alpha*s2
	+4.*a2*g22*h12*alpha*alpha*s2
	+4.*a2*g21*g11*g22*h12*alpha*alpha*s2
	+4.*a2*g12*h12*alpha*alpha*cth*sth
	+4.*a2*g21*h11*alpha*alpha*s2
	+4.*a2*g21*g11*g12*h12*alpha*alpha*cth*sth
	+8.*a2*g12_2*h12*g22*alpha*alpha*cth*sth
	-2.*a4*a4*h11*g11
	-2.*a4*a4*h12*g12
	-2.*a2*g22_2*h12*g12*beta*beta*c2
	+2.*a2*g12*h12*beta*beta*sth*cth
	+2.*a2*g22*h12*beta*beta*sth*cth
	-2.*a2*g21_2*h11*g11*beta*beta*c2
	+4.*a2*g12_2*h12*g22*beta*beta*sth*cth
	+2.*a2*g21*g11*g12*h12*beta*beta*sth*cth
	+2.*a2*g21*h11*g12_2*beta*beta*sth*cth
	-2.*a2*g22*h12*beta*beta*c2
	-2.*a2*g21*h11*beta*beta*c2
	-2.*a2*g11*h11*g12_2*beta*beta*s2
	-2.*a2*g11_2*g12*h12*beta*beta*s2
	+2.*a2*g11_2*g22*h12*beta*beta*sth*cth
	+2.*a2*g11*h11*g22*g12*beta*beta*sth*cth
	-2.*a2*g12*h12*beta*beta*s2
	+2.*a2*g11*h11*beta*beta*sth*cth
	-2.*a2*g11*h11*beta*beta*s2
	-2.*a2*g11_3*h11*beta*beta*s2
	+4.*a2*g11_2*h11*g21*beta*beta*sth*cth
	-2.*a2*g21*g11*g22*h12*beta*beta*c2
	+2.*a2*g21*h11*beta*beta*sth*cth
	-2.*a2*g12_3*h12*beta*beta*s2
	-2.*a2*g21*h11*g22*g12*beta*beta*c2)
	-4.*alpha*cth*(2.*a2*g11_2*g22*h12*beta*s2*alpha
	-2.*a2*g11_2*g22*h12*alpha*c2*beta
	-2.*a2*g21*h11*g12_2*alpha*c2*beta
	+2.*a2*g21*g11*g12*h12*beta*s2*alpha
	-4.*a2*g12_2*h12*g22*alpha*c2*beta
	+4.*a2*g12_2*h12*g22*beta*s2*alpha
	+2.*a2*g21*h11*g12_2*beta*s2*alpha
	-2.*a2*g11*h11*g22*g12*alpha*c2*beta
	+2.*a2*g11*h11*g22*g12*beta*s2*alpha
	+4.*a2*g11*h11*g12_2*alpha*sth*beta*cth
	-2.*a2*g22*h12*alpha*c2*beta
	+2.*a2*g22*h12*beta*s2*alpha
	-4.*a2*g22_2*h12*g12*alpha*sth*beta*cth
	-2.*a2*g12*h12*alpha*c2*beta
	+2.*a2*g12*h12*beta*s2*alpha
	-4.*a2*g22*h12*alpha*sth*beta*cth
	-2.*a2*g11*h11*alpha*c2*beta
	-4.*a2*g11_2*h11*g21*alpha*c2*beta
	+4.*a2*g11_2*h11*g21*beta*s2*alpha
	-4.*a2*g21*h11*g22*g12*alpha*sth*beta*cth
	+4.*a2*g12*h12*alpha*sth*beta*cth
	+4.*a2*g12_3*h12*alpha*sth*beta*cth
	+4.*a2*g11_3*h11*alpha*sth*beta*cth
	-2.*a2*g21*h11*alpha*c2*beta
	+2.*a2*g21*h11*beta*s2*alpha
	-4.*a2*g21*h11*alpha*sth*beta*cth
	-2.*a2*g21*g11*g12*h12*alpha*c2*beta
	+4.*a2*g11_2*g12*h12*alpha*sth*beta*cth
	+4.*a2*g11*h11*alpha*sth*beta*cth
	-4.*a2*g21*g11*g22*h12*alpha*sth*beta*cth
	+2.*a2*g11*h11*beta*s2*alpha
	-4.*a2*g21_2*h11*g11*alpha*sth*beta*cth)
	+2.*beta*sth*(a2*g21*g11*g22*h12*beta*beta*c2
	+a2*g12_3*h12*beta*beta*s2
	+a2*g21*h11*g22*g12*beta*beta*c2
	+a2*g11*h11*g12_2*beta*beta*s2
	-a2*g11*h11*beta*beta*sth*cth
	-a4*a4*h11*g11
	+a2*g11*h11*beta*beta*s2
	+a2*g22_2*h12*g12*beta*beta*c2
	+a2*g22*h12*beta*beta*c2
	+a2*g12*h12*beta*beta*s2
	-a2*g11*h11*g22*g12*beta*beta*sth*cth
	-a2*g12*h12*beta*beta*sth*cth
	-a2*g21*h11*g12_2*beta*beta*sth*cth
	-a2*g21*h11*beta*beta*sth*cth
	-a2*g11_2*g22*h12*beta*beta*sth*cth
	-2.*a2*g12_2*h12*g22*beta*beta*sth*cth
	-2.*a2*g11_2*h11*g21*beta*beta*sth*cth
	-a2*g22*h12*beta*beta*sth*cth
	+a2*g11_3*h11*beta*beta*s2
	+a2*g21_2*h11*g11*beta*beta*c2
	+a2*g21*h11*beta*beta*c2
	-a2*g21*g11*g12*h12*beta*beta*sth*cth
	+a2*g11_2*g12*h12*beta*beta*s2
	-a4*a4*h12*g12);

	B3 = -2.*beta*sth*(-2.*a2*g11_2*g22*h12*beta*s2*alpha
	+2.*a2*g11_2*g22*h12*alpha*c2*beta
	+2.*a2*g21*h11*g12_2*alpha*c2*beta
	-2.*a2*g21*g11*g12*h12*beta*s2*alpha
	+4.*a2*g12_2*h12*g22*alpha*c2*beta
	-4.*a2*g12_2*h12*g22*beta*s2*alpha
	-2.*a2*g21*h11*g12_2*beta*s2*alpha
	+2.*a2*g11*h11*g22*g12*alpha*c2*beta
	-2.*a2*g11*h11*g22*g12*beta*s2*alpha
	-4.*a2*g11*h11*g12_2*alpha*sth*beta*cth
	+2.*a2*g22*h12*alpha*c2*beta
	-2.*a2*g22*h12*beta*s2*alpha
	+4.*a2*g22_2*h12*g12*alpha*sth*beta*cth
	+2.*a2*g12*h12*alpha*c2*beta
	-2.*a2*g12*h12*beta*s2*alpha
	+4.*a2*g22*h12*alpha*sth*beta*cth
	+2.*a2*g11*h11*alpha*c2*beta
	+4.*a2*g11_2*h11*g21*alpha*c2*beta
	-4.*a2*g11_2*h11*g21*beta*s2*alpha
	+4.*a2*g21*h11*g22*g12*alpha*sth*beta*cth
	-4.*a2*g12*h12*alpha*sth*beta*cth
	-4.*a2*g12_3*h12*alpha*sth*beta*cth
	-4.*a2*g11_3*h11*alpha*sth*beta*cth
	+2.*a2*g21*h11*alpha*c2*beta
	-2.*a2*g21*h11*beta*s2*alpha
	+4.*a2*g21*h11*alpha*sth*beta*cth
	+2.*a2*g21*g11*g12*h12*alpha*c2*beta
	-4.*a2*g11_2*g12*h12*alpha*sth*beta*cth
	-4.*a2*g11*h11*alpha*sth*beta*cth
	+4.*a2*g21*g11*g22*h12*alpha*sth*beta*cth
	-2.*a2*g11*h11*beta*s2*alpha
	+4.*a2*g21_2*h11*g11*alpha*sth*beta*cth)
	-4.*alpha*cth*(4.*a2*g21_2*h11*g11*alpha*alpha*s2
	+4.*a2*g11_2*g12*h12*alpha*alpha*c2
	+4.*a2*g11*h11*alpha*alpha*c2
	+4.*a2*g11_2*g22*h12*alpha*alpha*cth*sth
	+4.*a2*g11*h11*alpha*alpha*cth*sth
	+8.*a2*g11_2*h11*g21*alpha*alpha*cth*sth
	+4.*a2*g12_3*h12*alpha*alpha*c2
	+4.*a2*g12*h12*alpha*alpha*c2
	+4.*a2*g11_3*h11*alpha*alpha*c2
	+4.*a2*g21*h11*g22*g12*alpha*alpha*s2
	+4.*a2*g21*h11*alpha*alpha*cth*sth
	+4.*a2*g11*h11*g22*g12*alpha*alpha*cth*sth
	+4.*a2*g11*h11*g12_2*alpha*alpha*c2
	+4.*a2*g21*h11*g12_2*alpha*alpha*cth*sth
	+4.*a2*g22*h12*alpha*alpha*cth*sth
	+4.*a2*g22_2*h12*g12*alpha*alpha*s2
	+4.*a2*g22*h12*alpha*alpha*s2
	+4.*a2*g21*g11*g22*h12*alpha*alpha*s2
	+4.*a2*g12*h12*alpha*alpha*cth*sth
	+4.*a2*g21*h11*alpha*alpha*s2
	+4.*a2*g21*g11*g12*h12*alpha*alpha*cth*sth
	+8.*a2*g12_2*h12*g22*alpha*alpha*cth*sth
	-2.*a4*a4*h11*g11
	-2.*a4*a4*h12*g12
	-2.*a2*g22_2*h12*g12*beta*beta*c2
	+2.*a2*g12*h12*beta*beta*sth*cth
	+2.*a2*g22*h12*beta*beta*sth*cth
	-2.*a2*g21_2*h11*g11*beta*beta*c2
	+4.*a2*g12_2*h12*g22*beta*beta*sth*cth
	+2.*a2*g21*g11*g12*h12*beta*beta*sth*cth
	+2.*a2*g21*h11*g12_2*beta*beta*sth*cth
	-2.*a2*g22*h12*beta*beta*c2
	-2.*a2*g21*h11*beta*beta*c2
	-2.*a2*g11*h11*g12_2*beta*beta*s2
	-2.*a2*g11_2*g12*h12*beta*beta*s2
	+2.*a2*g11_2*g22*h12*beta*beta*sth*cth
	+2.*a2*g11*h11*g22*g12*beta*beta*sth*cth
	-2.*a2*g12*h12*beta*beta*s2
	+2.*a2*g11*h11*beta*beta*sth*cth
	-2.*a2*g11*h11*beta*beta*s2
	-2.*a2*g11_3*h11*beta*beta*s2
	+4.*a2*g11_2*h11*g21*beta*beta*sth*cth
	-2.*a2*g21*g11*g22*h12*beta*beta*c2
	+2.*a2*g21*h11*beta*beta*sth*cth
	-2.*a2*g12_3*h12*beta*beta*s2
	-2.*a2*g21*h11*g22*g12*beta*beta*c2)
	+2.*beta*sth*(2.*a2*g11_2*g22*h12*beta*s2*alpha
	-2.*a2*g11_2*g22*h12*alpha*c2*beta
	-2.*a2*g21*h11*g12_2*alpha*c2*beta
	+2.*a2*g21*g11*g12*h12*beta*s2*alpha
	-4.*a2*g12_2*h12*g22*alpha*c2*beta
	+4.*a2*g12_2*h12*g22*beta*s2*alpha
	+2.*a2*g21*h11*g12_2*beta*s2*alpha
	-2.*a2*g11*h11*g22*g12*alpha*c2*beta
	+2.*a2*g11*h11*g22*g12*beta*s2*alpha
	+4.*a2*g11*h11*g12_2*alpha*sth*beta*cth
	-2.*a2*g22*h12*alpha*c2*beta
	+2.*a2*g22*h12*beta*s2*alpha
	-4.*a2*g22_2*h12*g12*alpha*sth*beta*cth
	-2.*a2*g12*h12*alpha*c2*beta
	+2.*a2*g12*h12*beta*s2*alpha
	-4.*a2*g22*h12*alpha*sth*beta*cth
	-2.*a2*g11*h11*alpha*c2*beta
	-4.*a2*g11_2*h11*g21*alpha*c2*beta
	+4.*a2*g11_2*h11*g21*beta*s2*alpha
	-4.*a2*g21*h11*g22*g12*alpha*sth*beta*cth
	+4.*a2*g12*h12*alpha*sth*beta*cth
	+4.*a2*g12_3*h12*alpha*sth*beta*cth
	+4.*a2*g11_3*h11*alpha*sth*beta*cth
	-2.*a2*g21*h11*alpha*c2*beta
	+2.*a2*g21*h11*beta*s2*alpha
	-4.*a2*g21*h11*alpha*sth*beta*cth
	-2.*a2*g21*g11*g12*h12*alpha*c2*beta
	+4.*a2*g11_2*g12*h12*alpha*sth*beta*cth
	+4.*a2*g11*h11*alpha*sth*beta*cth
	-4.*a2*g21*g11*g22*h12*alpha*sth*beta*cth
	+2.*a2*g11*h11*beta*s2*alpha
	-4.*a2*g21_2*h11*g11*alpha*sth*beta*cth);

	C0 = -beta*beta*s2*(-a4*a4*h12*h12
	+2.*a2*g21*h11*g22*h12*beta*beta*c2
	-2.*a2*g11*h11*h11*g21*beta*beta*sth*cth
	-a4*a4*h11*h11
	-2.*a2*g12*h12*h12*g22*beta*beta*sth*cth
	+a2*g21_2*h11*h11*beta*beta*c2
	-2.*a2*g21*h11*g12*h12*beta*beta*sth*cth
	+a2*g12_2*h12*h12*beta*beta*s2
	+a2*g22_2*h12*h12*beta*beta*c2
	+2.*a2*g11*h11*g12*h12*beta*beta*s2
	-2.*a2*g11*h11*g22*h12*beta*beta*sth*cth
	+a2*g11_2*h11*h11*beta*beta*s2);

	C1 = -beta*beta*s2*(8.*a2*g11*h11*g12*h12*alpha*sth*beta*cth
	-4.*a2*g22_2*h12*h12*alpha*sth*beta*cth
	-4.*a2*g21_2*h11*h11*alpha*sth*beta*cth
	-4.*a2*g21*h11*g12*h12*alpha*c2*beta
	-4.*a2*g12*h12*h12*g22*alpha*c2*beta
	+4.*a2*g21*h11*g12*h12*beta*s2*alpha
	-4.*a2*g11*h11*h11*g21*alpha*c2*beta
	-4.*a2*g11*h11*g22*h12*alpha*c2*beta
	+4.*a2*g11_2*h11*h11*alpha*sth*beta*cth
	+4.*a2*g11*h11*g22*h12*beta*s2*alpha
	+4.*a2*g12*h12*h12*g22*beta*s2*alpha
	-8.*a2*g21*h11*g22*h12*alpha*sth*beta*cth
	+4.*a2*g11*h11*h11*g21*beta*s2*alpha
	+4.*a2*g12_2*h12*h12*alpha*sth*beta*cth)
	-4.*beta*sth*alpha*cth*(-a4*a4*h12*h12
	+2.*a2*g21*h11*g22*h12*beta*beta*c2
	-2.*a2*g11*h11*h11*g21*beta*beta*sth*cth
	-a4*a4*h11*h11
	-2.*a2*g12*h12*h12*g22*beta*beta*sth*cth
	+a2*g21_2*h11*h11*beta*beta*c2
	-2.*a2*g21*h11*g12*h12*beta*beta*sth*cth
	+a2*g12_2*h12*h12*beta*beta*s2
	+a2*g22_2*h12*h12*beta*beta*c2
	+2.*a2*g11*h11*g12*h12*beta*beta*s2
	-2.*a2*g11*h11*g22*h12*beta*beta*sth*cth
	+a2*g11_2*h11*h11*beta*beta*s2);

	C2 = -beta*beta*s2*(-4.*a2*g11*h11*g12*h12*beta*beta*s2
	+4.*a2*g12*h12*h12*g22*beta*beta*sth*cth
	+4.*a2*g11*h11*h11*g21*beta*beta*sth*cth
	-2.*a4*a4*h12*h12
	+4.*a2*g21_2*h11*h11*alpha*alpha*s2
	+8.*a2*g21*h11*g12*h12*alpha*alpha*cth*sth
	+4.*a2*g12_2*h12*h12*alpha*alpha*c2
	-2.*a2*g22_2*h12*h12*beta*beta*c2
	-4.*a2*g21*h11*g22*h12*beta*beta*c2
	+8.*a2*g11*h11*g22*h12*alpha*alpha*cth*sth
	+4.*a2*g11*h11*g22*h12*beta*beta*sth*cth
	+8.*a2*g11*h11*g12*h12*alpha*alpha*c2
	+8.*a2*g21*h11*g22*h12*alpha*alpha*s2
	-2.*a2*g11_2*h11*h11*beta*beta*s2
	+4.*a2*g11_2*h11*h11*alpha*alpha*c2
	-2.*a4*a4*h11*h11
	+4.*a2*g22_2*h12*h12*alpha*alpha*s2
	+8.*a2*g12*h12*h12*g22*alpha*alpha*cth*sth
	+8.*a2*g11*h11*h11*g21*alpha*alpha*cth*sth
	+4.*a2*g21*h11*g12*h12*beta*beta*sth*cth
	-2.*a2*g12_2*h12*h12*beta*beta*s2
	-2.*a2*g21_2*h11*h11*beta*beta*c2)
	-4.*beta*sth*alpha*cth*(8.*a2*g11*h11*g12*h12*alpha*sth*beta*cth
	-4.*a2*g22_2*h12*h12*alpha*sth*beta*cth
	-4.*a2*g21_2*h11*h11*alpha*sth*beta*cth
	-4.*a2*g21*h11*g12*h12*alpha*c2*beta
	-4.*a2*g12*h12*h12*g22*alpha*c2*beta
	+4.*a2*g21*h11*g12*h12*beta*s2*alpha
	-4.*a2*g11*h11*h11*g21*alpha*c2*beta
	-4.*a2*g11*h11*g22*h12*alpha*c2*beta
	+4.*a2*g11_2*h11*h11*alpha*sth*beta*cth
	+4.*a2*g11*h11*g22*h12*beta*s2*alpha
	+4.*a2*g12*h12*h12*g22*beta*s2*alpha
	-8.*a2*g21*h11*g22*h12*alpha*sth*beta*cth
	+4.*a2*g11*h11*h11*g21*beta*s2*alpha
	+4.*a2*g12_2*h12*h12*alpha*sth*beta*cth)
	-(-2.*beta*beta*s2
	+4.*alpha*alpha*c2)*(-a4*a4*h12*h12
	+2.*a2*g21*h11*g22*h12*beta*beta*c2
	-2.*a2*g11*h11*h11*g21*beta*beta*sth*cth
	-a4*a4*h11*h11
	-2.*a2*g12*h12*h12*g22*beta*beta*sth*cth
	+a2*g21_2*h11*h11*beta*beta*c2
	-2.*a2*g21*h11*g12*h12*beta*beta*sth*cth
	+a2*g12_2*h12*h12*beta*beta*s2
	+a2*g22_2*h12*h12*beta*beta*c2
	+2.*a2*g11*h11*g12*h12*beta*beta*s2
	-2.*a2*g11*h11*g22*h12*beta*beta*sth*cth
	+a2*g11_2*h11*h11*beta*beta*s2);

	C3 = -beta*beta*s2*(4.*a2*g21*h11*g12*h12*alpha*c2*beta
	-8.*a2*g11*h11*g12*h12*alpha*sth*beta*cth
	+8.*a2*g21*h11*g22*h12*alpha*sth*beta*cth
	-4.*a2*g12*h12*h12*g22*beta*s2*alpha
	-4.*a2*g12_2*h12*h12*alpha*sth*beta*cth
	+4.*a2*g12*h12*h12*g22*alpha*c2*beta
	+4.*a2*g21_2*h11*h11*alpha*sth*beta*cth
	-4.*a2*g11*h11*h11*g21*beta*s2*alpha
	+4.*a2*g11*h11*g22*h12*alpha*c2*beta
	-4.*a2*g21*h11*g12*h12*beta*s2*alpha
	-4.*a2*g11_2*h11*h11*alpha*sth*beta*cth
	+4.*a2*g11*h11*h11*g21*alpha*c2*beta
	-4.*a2*g11*h11*g22*h12*beta*s2*alpha
	+4.*a2*g22_2*h12*h12*alpha*sth*beta*cth)
	-4.*beta*sth*alpha*cth*(-4.*a2*g11*h11*g12*h12*beta*beta*s2
	+4.*a2*g12*h12*h12*g22*beta*beta*sth*cth
	+4.*a2*g11*h11*h11*g21*beta*beta*sth*cth
	-2.*a4*a4*h12*h12
	+4.*a2*g21_2*h11*h11*alpha*alpha*s2
	+8.*a2*g21*h11*g12*h12*alpha*alpha*cth*sth
	+4.*a2*g12_2*h12*h12*alpha*alpha*c2
	-2.*a2*g22_2*h12*h12*beta*beta*c2
	-4.*a2*g21*h11*g22*h12*beta*beta*c2
	+8.*a2*g11*h11*g22*h12*alpha*alpha*cth*sth
	+4.*a2*g11*h11*g22*h12*beta*beta*sth*cth
	+8.*a2*g11*h11*g12*h12*alpha*alpha*c2
	+8.*a2*g21*h11*g22*h12*alpha*alpha*s2
	-2.*a2*g11_2*h11*h11*beta*beta*s2
	+4.*a2*g11_2*h11*h11*alpha*alpha*c2
	-2.*a4*a4*h11*h11
	+4.*a2*g22_2*h12*h12*alpha*alpha*s2
	+8.*a2*g12*h12*h12*g22*alpha*alpha*cth*sth
	+8.*a2*g11*h11*h11*g21*alpha*alpha*cth*sth
	+4.*a2*g21*h11*g12*h12*beta*beta*sth*cth
	-2.*a2*g12_2*h12*h12*beta*beta*s2
	-2.*a2*g21_2*h11*h11*beta*beta*c2)
	-(-2.*beta*beta*s2
	+4.*alpha*alpha*c2)*(8.*a2*g11*h11*g12*h12*alpha*sth*beta*cth
	-4.*a2*g22_2*h12*h12*alpha*sth*beta*cth
	-4.*a2*g21_2*h11*h11*alpha*sth*beta*cth
	-4.*a2*g21*h11*g12*h12*alpha*c2*beta
	-4.*a2*g12*h12*h12*g22*alpha*c2*beta
	+4.*a2*g21*h11*g12*h12*beta*s2*alpha
	-4.*a2*g11*h11*h11*g21*alpha*c2*beta
	-4.*a2*g11*h11*g22*h12*alpha*c2*beta
	+4.*a2*g11_2*h11*h11*alpha*sth*beta*cth
	+4.*a2*g11*h11*g22*h12*beta*s2*alpha
	+4.*a2*g12*h12*h12*g22*beta*s2*alpha
	-8.*a2*g21*h11*g22*h12*alpha*sth*beta*cth
	+4.*a2*g11*h11*h11*g21*beta*s2*alpha
	+4.*a2*g12_2*h12*h12*alpha*sth*beta*cth)
	+4.*beta*sth*alpha*cth*(-a4*a4*h12*h12
	+2.*a2*g21*h11*g22*h12*beta*beta*c2
	-2.*a2*g11*h11*h11*g21*beta*beta*sth*cth
	-a4*a4*h11*h11
	-2.*a2*g12*h12*h12*g22*beta*beta*sth*cth
	+a2*g21_2*h11*h11*beta*beta*c2
	-2.*a2*g21*h11*g12*h12*beta*beta*sth*cth
	+a2*g12_2*h12*h12*beta*beta*s2
	+a2*g22_2*h12*h12*beta*beta*c2
	+2.*a2*g11*h11*g12*h12*beta*beta*s2
	-2.*a2*g11*h11*g22*h12*beta*beta*sth*cth
	+a2*g11_2*h11*h11*beta*beta*s2);

	C4 = -2.*beta*beta*s2*(-a4*a4*h12*h12
	+2.*a2*g21*h11*g22*h12*beta*beta*c2
	-2.*a2*g11*h11*h11*g21*beta*beta*sth*cth
	-a4*a4*h11*h11
	-2.*a2*g12*h12*h12*g22*beta*beta*sth*cth
	+a2*g21_2*h11*h11*beta*beta*c2
	-2.*a2*g21*h11*g12*h12*beta*beta*sth*cth
	+a2*g12_2*h12*h12*beta*beta*s2
	+a2*g22_2*h12*h12*beta*beta*c2
	+2.*a2*g11*h11*g12*h12*beta*beta*s2
	-2.*a2*g11*h11*g22*h12*beta*beta*sth*cth
	+a2*g11_2*h11*h11*beta*beta*s2)
	-4.*beta*sth*alpha*cth*(4.*a2*g21*h11*g12*h12*alpha*c2*beta
	-8.*a2*g11*h11*g12*h12*alpha*sth*beta*cth
	+8.*a2*g21*h11*g22*h12*alpha*sth*beta*cth
	-4.*a2*g12*h12*h12*g22*beta*s2*alpha
	-4.*a2*g12_2*h12*h12*alpha*sth*beta*cth
	+4.*a2*g12*h12*h12*g22*alpha*c2*beta
	+4.*a2*g21_2*h11*h11*alpha*sth*beta*cth
	-4.*a2*g11*h11*h11*g21*beta*s2*alpha
	+4.*a2*g11*h11*g22*h12*alpha*c2*beta
	-4.*a2*g21*h11*g12*h12*beta*s2*alpha
	-4.*a2*g11_2*h11*h11*alpha*sth*beta*cth
	+4.*a2*g11*h11*h11*g21*alpha*c2*beta
	-4.*a2*g11*h11*g22*h12*beta*s2*alpha
	+4.*a2*g22_2*h12*h12*alpha*sth*beta*cth)
	-(-2.*beta*beta*s2
	+4.*alpha*alpha*c2)*(-4.*a2*g11*h11*g12*h12*beta*beta*s2
	+4.*a2*g12*h12*h12*g22*beta*beta*sth*cth
	+4.*a2*g11*h11*h11*g21*beta*beta*sth*cth
	-2.*a4*a4*h12*h12
	+4.*a2*g21_2*h11*h11*alpha*alpha*s2
	+8.*a2*g21*h11*g12*h12*alpha*alpha*cth*sth
	+4.*a2*g12_2*h12*h12*alpha*alpha*c2
	-2.*a2*g22_2*h12*h12*beta*beta*c2
	-4.*a2*g21*h11*g22*h12*beta*beta*c2
	+8.*a2*g11*h11*g22*h12*alpha*alpha*cth*sth
	+4.*a2*g11*h11*g22*h12*beta*beta*sth*cth
	+8.*a2*g11*h11*g12*h12*alpha*alpha*c2
	+8.*a2*g21*h11*g22*h12*alpha*alpha*s2
	-2.*a2*g11_2*h11*h11*beta*beta*s2
	+4.*a2*g11_2*h11*h11*alpha*alpha*c2
	-2.*a4*a4*h11*h11
	+4.*a2*g22_2*h12*h12*alpha*alpha*s2
	+8.*a2*g12*h12*h12*g22*alpha*alpha*cth*sth
	+8.*a2*g11*h11*h11*g21*alpha*alpha*cth*sth
	+4.*a2*g21*h11*g12*h12*beta*beta*sth*cth
	-2.*a2*g12_2*h12*h12*beta*beta*s2
	-2.*a2*g21_2*h11*h11*beta*beta*c2)
	+4.*beta*sth*alpha*cth*(8.*a2*g11*h11*g12*h12*alpha*sth*beta*cth
	-4.*a2*g22_2*h12*h12*alpha*sth*beta*cth
	-4.*a2*g21_2*h11*h11*alpha*sth*beta*cth
	-4.*a2*g21*h11*g12*h12*alpha*c2*beta
	-4.*a2*g12*h12*h12*g22*alpha*c2*beta
	+4.*a2*g21*h11*g12*h12*beta*s2*alpha
	-4.*a2*g11*h11*h11*g21*alpha*c2*beta
	-4.*a2*g11*h11*g22*h12*alpha*c2*beta
	+4.*a2*g11_2*h11*h11*alpha*sth*beta*cth
	+4.*a2*g11*h11*g22*h12*beta*s2*alpha
	+4.*a2*g12*h12*h12*g22*beta*s2*alpha
	-8.*a2*g21*h11*g22*h12*alpha*sth*beta*cth
	+4.*a2*g11*h11*h11*g21*beta*s2*alpha
	+4.*a2*g12_2*h12*h12*alpha*sth*beta*cth);

	E0 = 2.*a3*g21_2*g12*g22*beta*beta*cth*sth
	+2.*a3*g12*g22*beta*beta*cth*sth
	+2.*a3*g12*g22_3*beta*beta*cth*sth
	-a3*beta*beta*s2
	-2.*a3*g12*g22*beta*beta*s2
	-2.*a3*g11*g21*beta*beta*s2
	-a3*beta*beta*c2
	+2.*a3*g11*g21*beta*beta*cth*sth
	+2.*a3*g22_2*beta*beta*cth*sth
	-a3*g21_4*beta*beta*c2
	+2.*a3*g21_2*beta*beta*cth*sth
	-a3*g12_2*g22_2*beta*beta*s2
	-a3*g11_2*g21_2*beta*beta*s2
	-2.*a3*g21_2*g22_2*beta*beta*c2
	+a6*a6*g21_2
	+a6*a6*g22_2
	+2.*a3*g11*g21_3*beta*beta*cth*sth
	+2.*a3*g11*g21*g22_2*beta*beta*cth*sth
	-2.*a3*g11*g21*g12*g22*beta*beta*s2
	+a6*a6
	-2.*a3*g21_2*beta*beta*c2
	-2.*a3*g22_2*beta*beta*c2
	-a3*g22_4*beta*beta*c2
	+2.*a3*beta*beta*cth*sth;

	E1 = -4.*a3*g11_2*g21_2*alpha*sth*beta*cth
	+8.*a3*g21_2*g22_2*alpha*sth*beta*cth
	-4.*a3*g12_2*g22_2*alpha*sth*beta*cth
	+4.*a3*g22_2*beta*c2*alpha
	-4.*a3*g22_2*alpha*s2*beta
	-8.*a3*g11*g21*alpha*sth*beta*cth
	-4.*a3*g12*g22*alpha*s2*beta
	+4.*a3*g12*g22*beta*c2*alpha
	+4.*a3*g21_2*g12*g22*beta*c2*alpha
	-4.*a3*g12*g22_3*alpha*s2*beta
	+4.*a3*g12*g22_3*beta*c2*alpha
	+8.*a3*g21_2*alpha*sth*beta*cth
	+4.*a3*g21_4*alpha*sth*beta*cth
	-4.*a3*g11*g21*alpha*s2*beta
	+4.*a3*g11*g21*beta*c2*alpha
	-8.*a3*g12*g22*alpha*sth*beta*cth
	-4.*a3*g21_2*g12*g22*alpha*s2*beta
	+4.*a3*g11*g21*g22_2*beta*c2*alpha
	-8.*a3*g11*g21*g12*g22*alpha*sth*beta*cth
	-4.*a3*g11*g21_3*alpha*s2*beta
	+4.*a3*g11*g21_3*beta*c2*alpha
	+4.*a3*g21_2*beta*c2*alpha
	+8.*a3*g22_2*alpha*sth*beta*cth
	-4.*a3*alpha*s2*beta
	-4.*a3*g21_2*alpha*s2*beta
	+4.*a3*g22_4*alpha*sth*beta*cth
	-4.*a3*g11*g21*g22_2*alpha*s2*beta
	+4.*a3*beta*c2*alpha;

	E2 = -4.*a3*g21_2*g12*g22*beta*beta*cth*sth
	-4.*a3*g12*g22*beta*beta*cth*sth
	-4.*a3*g12*g22_3*beta*beta*cth*sth
	+2.*a3*beta*beta*s2
	+4.*a3*g12*g22*beta*beta*s2
	+4.*a3*g11*g21*beta*beta*s2
	+2.*a3*beta*beta*c2
	-4.*a3*g11*g21*beta*beta*cth*sth
	-4.*a3*g22_2*beta*beta*cth*sth
	+2.*a3*g21_4*beta*beta*c2
	-4.*a3*g21_2*beta*beta*cth*sth
	+2.*a3*g12_2*g22_2*beta*beta*s2
	+2.*a3*g11_2*g21_2*beta*beta*s2
	+4.*a3*g21_2*g22_2*beta*beta*c2
	+2.*a6*a6*g21_2
	+2.*a6*a6*g22_2
	-4.*a3*g11*g21_3*beta*beta*cth*sth
	-4.*a3*g11*g21*g22_2*beta*beta*cth*sth
	+4.*a3*g11*g21*g12*g22*beta*beta*s2
	+2.*a6*a6
	-4.*a3*g11_2*g21_2*alpha*alpha*c2
	-8.*a3*g12*g22*alpha*alpha*c2
	-8.*a3*g11*g21*alpha*alpha*c2
	-8.*a3*g12*g22*alpha*alpha*sth*cth
	-8.*a3*g12*g22_3*alpha*alpha*sth*cth
	-8.*a3*g21_2*g22_2*alpha*alpha*s2
	-8.*a3*g22_2*alpha*alpha*sth*cth
	-8.*a3*g21_2*g12*g22*alpha*alpha*sth*cth
	-8.*a3*g11*g21_3*alpha*alpha*sth*cth
	-4.*a3*alpha*alpha*s2
	-8.*a3*g11*g21*g12*g22*alpha*alpha*c2
	-4.*a3*g12_2*g22_2*alpha*alpha*c2
	-8.*a3*g21_2*alpha*alpha*sth*cth
	-8.*a3*g11*g21*alpha*alpha*sth*cth
	-4.*a3*alpha*alpha*c2
	-8.*a3*g11*g21*g22_2*alpha*alpha*sth*cth
	+4.*a3*g21_2*beta*beta*c2
	+4.*a3*g22_2*beta*beta*c2
	+2.*a3*g22_4*beta*beta*c2
	-4.*a3*beta*beta*cth*sth
	-4.*a3*g21_4*alpha*alpha*s2
	-8.*a3*g21_2*alpha*alpha*s2
	-8.*a3*alpha*alpha*sth*cth
	-4.*a3*g22_4*alpha*alpha*s2
	-8.*a3*g22_2*alpha*alpha*s2;

	F0 = -2.*beta*cth*(-a6*a6*h22*g22
	-a6*a6*h21*g21
	-a3*g11*h21*g22_2*beta*beta*cth*sth
	+a3*g11*h21*g12*g22*beta*beta*s2
	-2.*a3*g11*h21*g21_2*beta*beta*cth*sth
	+a3*g22_3*h22*beta*beta*c2
	+a3*g22*h22*beta*beta*c2
	-a3*g21*h21*g12*g22*beta*beta*cth*sth
	+a3*g12*h22*beta*beta*s2
	+a3*g12_2*h22*g22*beta*beta*s2
	+a3*g21*h21*g22_2*beta*beta*c2
	-a3*g21_2*g12*h22*beta*beta*cth*sth
	-a3*g11*h21*beta*beta*cth*sth
	+a3*g11*h21*beta*beta*s2
	-a3*g11*g21*g22*h22*beta*beta*cth*sth
	-a3*g22*h22*beta*beta*cth*sth
	-2.*a3*g12*h22*g22_2*beta*beta*cth*sth
	+a3*g21*h21*beta*beta*c2
	+a3*g21_2*g22*h22*beta*beta*c2
	-a3*g21*h21*beta*beta*cth*sth
	-a3*g12*h22*beta*beta*cth*sth
	+a3*g11_2*h21*g21*beta*beta*s2
	+a3*g21_3*h21*beta*beta*c2
	+a3*g11*g21*g12*h22*beta*beta*s2);

	F1 = -2.*beta*cth*(-4.*a3*g21*h21*alpha*sth*beta*cth
	-4.*a3*g21_2*g22*h22*alpha*sth*beta*cth
	-4.*a3*g21_3*h21*alpha*sth*beta*cth
	+2.*a3*g21*h21*alpha*s2*beta
	-2.*a3*g21*h21*beta*c2*alpha
	+2.*a3*g12*h22*alpha*s2*beta
	+4.*a3*g12_2*h22*g22*alpha*sth*beta*cth
	+4.*a3*g12*h22*g22_2*alpha*s2*beta
	-4.*a3*g12*h22*g22_2*beta*c2*alpha
	-4.*a3*g21*h21*g22_2*alpha*sth*beta*cth
	+2.*a3*g21_2*g12*h22*alpha*s2*beta
	-2.*a3*g21_2*g12*h22*beta*c2*alpha
	-2.*a3*g11*h21*beta*c2*alpha
	+4.*a3*g11*h21*alpha*sth*beta*cth
	-4.*a3*g22*h22*alpha*sth*beta*cth
	+2.*a3*g21*h21*g12*g22*alpha*s2*beta
	-2.*a3*g21*h21*g12*g22*beta*c2*alpha
	+4.*a3*g11*g21*g12*h22*alpha*sth*beta*cth
	+2.*a3*g11*g21*g22*h22*alpha*s2*beta
	-2.*a3*g11*g21*g22*h22*beta*c2*alpha
	+4.*a3*g12*h22*alpha*sth*beta*cth
	+4.*a3*g11*h21*g12*g22*alpha*sth*beta*cth
	-4.*a3*g11*h21*g21_2*beta*c2*alpha
	-4.*a3*g22_3*h22*alpha*sth*beta*cth
	-2.*a3*g11*h21*g22_2*beta*c2*alpha
	+2.*a3*g11*h21*alpha*s2*beta
	-2.*a3*g12*h22*beta*c2*alpha
	+4.*a3*g11*h21*g21_2*alpha*s2*beta
	+2.*a3*g11*h21*g22_2*alpha*s2*beta
	+4.*a3*g11_2*h21*g21*alpha*sth*beta*cth
	+2.*a3*g22*h22*alpha*s2*beta
	-2.*a3*g22*h22*beta*c2*alpha)
	+4.*alpha*sth*(-a6*a6*h22*g22
	-a6*a6*h21*g21
	-a3*g11*h21*g22_2*beta*beta*cth*sth
	+a3*g11*h21*g12*g22*beta*beta*s2
	-2.*a3*g11*h21*g21_2*beta*beta*cth*sth
	+a3*g22_3*h22*beta*beta*c2
	+a3*g22*h22*beta*beta*c2
	-a3*g21*h21*g12*g22*beta*beta*cth*sth
	+a3*g12*h22*beta*beta*s2
	+a3*g12_2*h22*g22*beta*beta*s2
	+a3*g21*h21*g22_2*beta*beta*c2
	-a3*g21_2*g12*h22*beta*beta*cth*sth
	-a3*g11*h21*beta*beta*cth*sth
	+a3*g11*h21*beta*beta*s2
	-a3*g11*g21*g22*h22*beta*beta*cth*sth
	-a3*g22*h22*beta*beta*cth*sth
	-2.*a3*g12*h22*g22_2*beta*beta*cth*sth
	+a3*g21*h21*beta*beta*c2
	+a3*g21_2*g22*h22*beta*beta*c2
	-a3*g21*h21*beta*beta*cth*sth
	-a3*g12*h22*beta*beta*cth*sth
	+a3*g11_2*h21*g21*beta*beta*s2
	+a3*g21_3*h21*beta*beta*c2
	+a3*g11*g21*g12*h22*beta*beta*s2);

	F2 = -2.*beta*cth*(-(2*a6*a6*h22*g22)
	-(2*a6*a6*h21*g21)
	+2.*a3*g11*h21*(g22_2)*beta*beta*cth*sth
	-2.*a3*g11*h21*g12*g22*beta*beta*s2
	+4.*a3*g11*h21*(g21_2)*beta*beta*cth*sth
	-2.*a3*g22_3*h22*beta*beta*c2
	-2.*a3*g22*h22*beta*beta*c2
	+2.*a3*g21*h21*g12*g22*beta*beta*cth*sth
	-2.*a3*g12*h22*beta*beta*s2
	-2.*a3*g12_2*h22*g22*beta*beta*s2
	-2.*a3*g21*h21*(g22_2)*beta*beta*c2
	+2.*a3*(g21_2)*g12*h22*beta*beta*cth*sth
	+2.*a3*g11*h21*beta*beta*cth*sth
	-2.*a3*g11*h21*beta*beta*s2
	+2.*a3*g11*g21*g22*h22*beta*beta*cth*sth
	+2.*a3*g22*h22*beta*beta*cth*sth
	+4.*a3*g12*h22*(g22_2)*beta*beta*cth*sth
	-2.*a3*g21*h21*beta*beta*c2
	-2.*a3*(g21_2)*g22*h22*beta*beta*c2
	+2.*a3*g21*h21*beta*beta*cth*sth
	+2.*a3*g12*h22*beta*beta*cth*sth
	-2.*a3*g11_2*h21*g21*beta*beta*s2
	-2.*a3*g21_3*h21*beta*beta*c2
	-2.*a3*g11*g21*g12*h22*beta*beta*s2
	+4.*a3*g21*h21*alpha*alpha*s2
	+4.*a3*(g21_2)*g22*h22*alpha*alpha*s2
	+4.*a3*g21_3*h21*alpha*alpha*s2
	+4.*a3*g21*h21*alpha*alpha*sth*cth
	+4.*a3*g12*h22*alpha*alpha*sth*cth
	+4.*a3*g12_2*h22*g22*alpha*alpha*c2
	+4.*a3*(g21_2)*g12*h22*alpha*alpha*sth*cth
	+4.*a3*g11*h21*alpha*alpha*c2
	+4.*a3*g11*h21*alpha*alpha*sth*cth
	+4.*a3*g21*h21*g12*g22*alpha*alpha*sth*cth
	+4.*a3*g11*g21*g12*h22*alpha*alpha*c2
	+4.*a3*g11*h21*g12*g22*alpha*alpha*c2
	+4.*a3*g11*h21*(g22_2)*alpha*alpha*sth*cth
	+8.*a3*g11*h21*(g21_2)*alpha*alpha*sth*cth
	+4.*a3*g22*h22*alpha*alpha*s2
	+4.*a3*g21*h21*(g22_2)*alpha*alpha*s2
	+4.*a3*g22_3*h22*alpha*alpha*s2
	+4.*a3*g11_2*h21*g21*alpha*alpha*c2
	+4.*a3*g12*h22*alpha*alpha*c2
	+4.*a3*g11*g21*g22*h22*alpha*alpha*sth*cth
	+4.*a3*g22*h22*alpha*alpha*sth*cth
	+8.*a3*g12*h22*(g22_2)*alpha*alpha*sth*cth)
	+4.*alpha*sth*(-4.*a3*g21*h21*alpha*sth*beta*cth
	-4.*a3*(g21_2)*g22*h22*alpha*sth*beta*cth
	-4.*a3*g21_3*h21*alpha*sth*beta*cth
	+2.*a3*g21*h21*alpha*s2*beta
	-2.*a3*g21*h21*beta*c2*alpha
	+2.*a3*g12*h22*alpha*s2*beta
	+4.*a3*g12_2*h22*g22*alpha*sth*beta*cth
	+4.*a3*g12*h22*(g22_2)*alpha*s2*beta
	-4.*a3*g12*h22*(g22_2)*beta*c2*alpha
	-4.*a3*g21*h21*(g22_2)*alpha*sth*beta*cth
	+2.*a3*(g21_2)*g12*h22*alpha*s2*beta
	-2.*a3*(g21_2)*g12*h22*beta*c2*alpha
	-2.*a3*g11*h21*beta*c2*alpha
	+4.*a3*g11*h21*alpha*sth*beta*cth
	-4.*a3*g22*h22*alpha*sth*beta*cth
	+2.*a3*g21*h21*g12*g22*alpha*s2*beta
	-2.*a3*g21*h21*g12*g22*beta*c2*alpha
	+4.*a3*g11*g21*g12*h22*alpha*sth*beta*cth
	+2.*a3*g11*g21*g22*h22*alpha*s2*beta
	-2.*a3*g11*g21*g22*h22*beta*c2*alpha
	+4.*a3*g12*h22*alpha*sth*beta*cth
	+4.*a3*g11*h21*g12*g22*alpha*sth*beta*cth
	-4.*a3*g11*h21*(g21_2)*beta*c2*alpha
	-4.*a3*g22_3*h22*alpha*sth*beta*cth
	-2.*a3*g11*h21*(g22_2)*beta*c2*alpha
	+2.*a3*g11*h21*alpha*s2*beta
	-2.*a3*g12*h22*beta*c2*alpha
	+4.*a3*g11*h21*(g21_2)*alpha*s2*beta
	+2.*a3*g11*h21*(g22_2)*alpha*s2*beta
	+4.*a3*g11_2*h21*g21*alpha*sth*beta*cth
	+2.*a3*g22*h22*alpha*s2*beta
	-2.*a3*g22*h22*beta*c2*alpha)
	+2.*beta*cth*(-(a6*a6*h22*g22)
	-(a6*a6*h21*g21)
	-a3*g11*h21*(g22_2)*beta*beta*cth*sth
	+a3*g11*h21*g12*g22*beta*beta*s2
	-2.*a3*g11*h21*(g21_2)*beta*beta*cth*sth
	+a3*g22_3*h22*beta*beta*c2
	+a3*g22*h22*beta*beta*c2
	-a3*g21*h21*g12*g22*beta*beta*cth*sth
	+a3*g12*h22*beta*beta*s2
	+a3*g12_2*h22*g22*beta*beta*s2
	+a3*g21*h21*(g22_2)*beta*beta*c2
	-a3*(g21_2)*g12*h22*beta*beta*cth*sth
	-a3*g11*h21*beta*beta*cth*sth
	+a3*g11*h21*beta*beta*s2
	-a3*g11*g21*g22*h22*beta*beta*cth*sth
	-a3*g22*h22*beta*beta*cth*sth
	-2.*a3*g12*h22*(g22_2)*beta*beta*cth*sth
	+a3*g21*h21*beta*beta*c2
	+a3*(g21_2)*g22*h22*beta*beta*c2
	-a3*g21*h21*beta*beta*cth*sth
	-a3*g12*h22*beta*beta*cth*sth
	+a3*g11_2*h21*g21*beta*beta*s2
	+a3*g21_3*h21*beta*beta*c2
	+a3*g11*g21*g12*h22*beta*beta*s2);

	F3 = -2.*beta*cth*(4.*a3*g21*h21*alpha*sth*beta*cth
	+4.*a3*g21_2*g22*h22*alpha*sth*beta*cth
	+4.*a3*g21_3*h21*alpha*sth*beta*cth
	-2.*a3*g21*h21*alpha*s2*beta
	+2.*a3*g21*h21*beta*c2*alpha
	-2.*a3*g12*h22*alpha*s2*beta
	-4.*a3*g12_2*h22*g22*alpha*sth*beta*cth
	-4.*a3*g12*h22*g22_2*alpha*s2*beta
	+4.*a3*g12*h22*g22_2*beta*c2*alpha
	+4.*a3*g21*h21*g22_2*alpha*sth*beta*cth
	-2.*a3*g21_2*g12*h22*alpha*s2*beta
	+2.*a3*g21_2*g12*h22*beta*c2*alpha
	+2.*a3*g11*h21*beta*c2*alpha
	-4.*a3*g11*h21*alpha*sth*beta*cth
	+4.*a3*g22*h22*alpha*sth*beta*cth
	-2.*a3*g21*h21*g12*g22*alpha*s2*beta
	+2.*a3*g21*h21*g12*g22*beta*c2*alpha
	-4.*a3*g11*g21*g12*h22*alpha*sth*beta*cth
	-2.*a3*g11*g21*g22*h22*alpha*s2*beta
	+2.*a3*g11*g21*g22*h22*beta*c2*alpha
	-4.*a3*g12*h22*alpha*sth*beta*cth
	-4.*a3*g11*h21*g12*g22*alpha*sth*beta*cth
	+4.*a3*g11*h21*g21_2*beta*c2*alpha
	+4.*a3*g22_3*h22*alpha*sth*beta*cth
	+2.*a3*g11*h21*g22_2*beta*c2*alpha
	-2.*a3*g11*h21*alpha*s2*beta
	+2.*a3*g12*h22*beta*c2*alpha
	-4.*a3*g11*h21*g21_2*alpha*s2*beta
	-2.*a3*g11*h21*g22_2*alpha*s2*beta
	-4.*a3*g11_2*h21*g21*alpha*sth*beta*cth
	-2.*a3*g22*h22*alpha*s2*beta
	+2.*a3*g22*h22*beta*c2*alpha)
	+4.*alpha*sth*(-2.*a6*a6*h22*g22
	-2.*a6*a6*h21*g21
	+2.*a3*g11*h21*g22_2*beta*beta*cth*sth
	-2.*a3*g11*h21*g12*g22*beta*beta*s2
	+4.*a3*g11*h21*g21_2*beta*beta*cth*sth
	-2.*a3*g22_3*h22*beta*beta*c2
	-2.*a3*g22*h22*beta*beta*c2
	+2.*a3*g21*h21*g12*g22*beta*beta*cth*sth
	-2.*a3*g12*h22*beta*beta*s2
	-2.*a3*g12_2*h22*g22*beta*beta*s2
	-2.*a3*g21*h21*g22_2*beta*beta*c2
	+2.*a3*g21_2*g12*h22*beta*beta*cth*sth
	+2.*a3*g11*h21*beta*beta*cth*sth
	-2.*a3*g11*h21*beta*beta*s2
	+2.*a3*g11*g21*g22*h22*beta*beta*cth*sth
	+2.*a3*g22*h22*beta*beta*cth*sth
	+4.*a3*g12*h22*g22_2*beta*beta*cth*sth
	-2.*a3*g21*h21*beta*beta*c2
	-2.*a3*g21_2*g22*h22*beta*beta*c2
	+2.*a3*g21*h21*beta*beta*cth*sth
	+2.*a3*g12*h22*beta*beta*cth*sth
	-2.*a3*g11_2*h21*g21*beta*beta*s2
	-2.*a3*g21_3*h21*beta*beta*c2
	-2.*a3*g11*g21*g12*h22*beta*beta*s2
	+4.*a3*g21*h21*alpha*alpha*s2
	+4.*a3*g21_2*g22*h22*alpha*alpha*s2
	+4.*a3*g21_3*h21*alpha*alpha*s2
	+4.*a3*g21*h21*alpha*alpha*sth*cth
	+4.*a3*g12*h22*alpha*alpha*sth*cth
	+4.*a3*g12_2*h22*g22*alpha*alpha*c2
	+4.*a3*g21_2*g12*h22*alpha*alpha*sth*cth
	+4.*a3*g11*h21*alpha*alpha*c2
	+4.*a3*g11*h21*alpha*alpha*sth*cth
	+4.*a3*g21*h21*g12*g22*alpha*alpha*sth*cth
	+4.*a3*g11*g21*g12*h22*alpha*alpha*c2
	+4.*a3*g11*h21*g12*g22*alpha*alpha*c2
	+4.*a3*g11*h21*g22_2*alpha*alpha*sth*cth
	+8.*a3*g11*h21*g21_2*alpha*alpha*sth*cth
	+4.*a3*g22*h22*alpha*alpha*s2
	+4.*a3*g21*h21*g22_2*alpha*alpha*s2
	+4.*a3*g22_3*h22*alpha*alpha*s2
	+4.*a3*g11_2*h21*g21*alpha*alpha*c2
	+4.*a3*g12*h22*alpha*alpha*c2
	+4.*a3*g11*g21*g22*h22*alpha*alpha*sth*cth
	+4.*a3*g22*h22*alpha*alpha*sth*cth
	+8.*a3*g12*h22*g22_2*alpha*alpha*sth*cth)
	+2.*beta*cth*(-4.*a3*g21*h21*alpha*sth*beta*cth
	-4.*a3*g21_2*g22*h22*alpha*sth*beta*cth
	-4.*a3*g21_3*h21*alpha*sth*beta*cth
	+2.*a3*g21*h21*alpha*s2*beta
	-2.*a3*g21*h21*beta*c2*alpha
	+2.*a3*g12*h22*alpha*s2*beta
	+4.*a3*g12_2*h22*g22*alpha*sth*beta*cth
	+4.*a3*g12*h22*g22_2*alpha*s2*beta
	-4.*a3*g12*h22*g22_2*beta*c2*alpha
	-4.*a3*g21*h21*g22_2*alpha*sth*beta*cth
	+2.*a3*g21_2*g12*h22*alpha*s2*beta
	-2.*a3*g21_2*g12*h22*beta*c2*alpha
	-2.*a3*g11*h21*beta*c2*alpha
	+4.*a3*g11*h21*alpha*sth*beta*cth
	-4.*a3*g22*h22*alpha*sth*beta*cth
	+2.*a3*g21*h21*g12*g22*alpha*s2*beta
	-2.*a3*g21*h21*g12*g22*beta*c2*alpha
	+4.*a3*g11*g21*g12*h22*alpha*sth*beta*cth
	+2.*a3*g11*g21*g22*h22*alpha*s2*beta
	-2.*a3*g11*g21*g22*h22*beta*c2*alpha
	+4.*a3*g12*h22*alpha*sth*beta*cth
	+4.*a3*g11*h21*g12*g22*alpha*sth*beta*cth
	-4.*a3*g11*h21*g21_2*beta*c2*alpha
	-4.*a3*g22_3*h22*alpha*sth*beta*cth
	-2.*a3*g11*h21*g22_2*beta*c2*alpha
	+2.*a3*g11*h21*alpha*s2*beta
	-2.*a3*g12*h22*beta*c2*alpha
	+4.*a3*g11*h21*g21_2*alpha*s2*beta
	+2.*a3*g11*h21*g22_2*alpha*s2*beta
	+4.*a3*g11_2*h21*g21*alpha*sth*beta*cth
	+2.*a3*g22*h22*alpha*s2*beta
	-2.*a3*g22*h22*beta*c2*alpha);

	G0 = -beta*beta*c2*(-a6*a6*h21*h21
	-a6*a6*h22*h22
	+a3*g12_2*h22*h22*beta*beta*s2
	+2.*a3*g11*h21*g12*h22*beta*beta*s2
	+a3*g11_2*h21*h21*beta*beta*s2
	-2.*a3*g11*h21*h21*g21*beta*beta*cth*sth
	+a3*g21_2*h21*h21*beta*beta*c2
	-2.*a3*g12*h22*h22*g22*beta*beta*cth*sth
	-2.*a3*g21*h21*g12*h22*beta*beta*cth*sth
	+2.*a3*g21*h21*g22*h22*beta*beta*c2
	+a3*g22_2*h22*h22*beta*beta*c2
	-2.*a3*g11*h21*g22*h22*beta*beta*cth*sth);

	G1 = -beta*beta*c2*(-4.*a3*g21*h21*g12*h22*beta*c2*alpha
	-4.*a3*g11*h21*h21*g21*beta*c2*alpha
	+4.*a3*g11_2*h21*h21*alpha*sth*beta*cth
	-4.*a3*g21_2*h21*h21*alpha*sth*beta*cth
	+4.*a3*g12_2*h22*h22*alpha*sth*beta*cth
	+4.*a3*g11*h21*g22*h22*alpha*s2*beta
	-4.*a3*g12*h22*h22*g22*beta*c2*alpha
	+4.*a3*g12*h22*h22*g22*alpha*s2*beta
	-4.*a3*g11*h21*g22*h22*beta*c2*alpha
	-8.*a3*g21*h21*g22*h22*alpha*sth*beta*cth
	-4.*a3*g22_2*h22*h22*alpha*sth*beta*cth
	+8.*a3*g11*h21*g12*h22*alpha*sth*beta*cth
	+4.*a3*g21*h21*g12*h22*alpha*s2*beta
	+4.*a3*g11*h21*h21*g21*alpha*s2*beta)
	+4.*beta*cth*alpha*sth*(-a6*a6*h21*h21
	-a6*a6*h22*h22
	+a3*g12_2*h22*h22*beta*beta*s2
	+2.*a3*g11*h21*g12*h22*beta*beta*s2
	+a3*g11_2*h21*h21*beta*beta*s2
	-2.*a3*g11*h21*h21*g21*beta*beta*cth*sth
	+a3*g21_2*h21*h21*beta*beta*c2
	-2.*a3*g12*h22*h22*g22*beta*beta*cth*sth
	-2.*a3*g21*h21*g12*h22*beta*beta*cth*sth
	+2.*a3*g21*h21*g22*h22*beta*beta*c2
	+a3*g22_2*h22*h22*beta*beta*c2
	-2.*a3*g11*h21*g22*h22*beta*beta*cth*sth);

	G2 = -beta*beta*c2*(-4.*a3*g11*h21*g12*h22*beta*beta*s2
	+8.*a3*g11*h21*g12*h22*alpha*alpha*c2
	-2.*a3*g22_2*h22*h22*beta*beta*c2
	+4.*a3*g12*h22*h22*g22*beta*beta*cth*sth
	-2.*a6*a6*h21*h21
	+8.*a3*g21*h21*g12*h22*alpha*alpha*sth*cth
	-2.*a3*g11_2*h21*h21*beta*beta*s2
	-2.*a6*a6*h22*h22
	+8.*a3*g11*h21*h21*g21*alpha*alpha*sth*cth
	+4.*a3*g11*h21*h21*g21*beta*beta*cth*sth
	+4.*a3*g22_2*h22*h22*alpha*alpha*s2
	+4.*a3*g12_2*h22*h22*alpha*alpha*c2
	+4.*a3*g21_2*h21*h21*alpha*alpha*s2
	-4.*a3*g21*h21*g22*h22*beta*beta*c2
	-2.*a3*g12_2*h22*h22*beta*beta*s2
	+8.*a3*g12*h22*h22*g22*alpha*alpha*sth*cth
	-2.*a3*g21_2*h21*h21*beta*beta*c2
	+4.*a3*g11*h21*g22*h22*beta*beta*cth*sth
	+8.*a3*g21*h21*g22*h22*alpha*alpha*s2
	+4.*a3*g21*h21*g12*h22*beta*beta*cth*sth
	+4.*a3*g11_2*h21*h21*alpha*alpha*c2
	+8.*a3*g11*h21*g22*h22*alpha*alpha*sth*cth)
	+4.*beta*cth*alpha*sth*(-4.*a3*g21*h21*g12*h22*beta*c2*alpha
	-4.*a3*g11*h21*h21*g21*beta*c2*alpha
	+4.*a3*g11_2*h21*h21*alpha*sth*beta*cth
	-4.*a3*g21_2*h21*h21*alpha*sth*beta*cth
	+4.*a3*g12_2*h22*h22*alpha*sth*beta*cth
	+4.*a3*g11*h21*g22*h22*alpha*s2*beta
	-4.*a3*g12*h22*h22*g22*beta*c2*alpha
	+4.*a3*g12*h22*h22*g22*alpha*s2*beta
	-4.*a3*g11*h21*g22*h22*beta*c2*alpha
	-8.*a3*g21*h21*g22*h22*alpha*sth*beta*cth
	-4.*a3*g22_2*h22*h22*alpha*sth*beta*cth
	+8.*a3*g11*h21*g12*h22*alpha*sth*beta*cth
	+4.*a3*g21*h21*g12*h22*alpha*s2*beta
	+4.*a3*g11*h21*h21*g21*alpha*s2*beta)
	-(-2.*beta*beta*c2
	+4.*alpha*alpha*s2)*(-a6*a6*h21*h21
	-a6*a6*h22*h22
	+a3*g12_2*h22*h22*beta*beta*s2
	+2.*a3*g11*h21*g12*h22*beta*beta*s2
	+a3*g11_2*h21*h21*beta*beta*s2
	-2.*a3*g11*h21*h21*g21*beta*beta*cth*sth
	+a3*g21_2*h21*h21*beta*beta*c2
	-2.*a3*g12*h22*h22*g22*beta*beta*cth*sth
	-2.*a3*g21*h21*g12*h22*beta*beta*cth*sth
	+2.*a3*g21*h21*g22*h22*beta*beta*c2
	+a3*g22_2*h22*h22*beta*beta*c2
	-2.*a3*g11*h21*g22*h22*beta*beta*cth*sth);

	G3 = -beta*beta*c2*(4.*a3*g22_2*h22*h22*alpha*sth*beta*cth
	+4.*a3*g21_2*h21*h21*alpha*sth*beta*cth
	+4.*a3*g11*h21*h21*g21*beta*c2*alpha
	-4.*a3*g12*h22*h22*g22*alpha*s2*beta
	+4.*a3*g11*h21*g22*h22*beta*c2*alpha
	-4.*a3*g11_2*h21*h21*alpha*sth*beta*cth
	+8.*a3*g21*h21*g22*h22*alpha*sth*beta*cth
	-4.*a3*g12_2*h22*h22*alpha*sth*beta*cth
	+4.*a3*g12*h22*h22*g22*beta*c2*alpha
	-4.*a3*g11*h21*g22*h22*alpha*s2*beta
	-8.*a3*g11*h21*g12*h22*alpha*sth*beta*cth
	-4.*a3*g11*h21*h21*g21*alpha*s2*beta
	-4.*a3*g21*h21*g12*h22*alpha*s2*beta
	+4.*a3*g21*h21*g12*h22*beta*c2*alpha)
	+4.*beta*cth*alpha*sth*(-4.*a3*g11*h21*g12*h22*beta*beta*s2
	+8.*a3*g11*h21*g12*h22*alpha*alpha*c2
	-2.*a3*g22_2*h22*h22*beta*beta*c2
	+4.*a3*g12*h22*h22*g22*beta*beta*cth*sth
	-2.*a6*a6*h21*h21
	+8.*a3*g21*h21*g12*h22*alpha*alpha*sth*cth
	-2.*a3*g11_2*h21*h21*beta*beta*s2
	-2.*a6*a6*h22*h22
	+8.*a3*g11*h21*h21*g21*alpha*alpha*sth*cth
	+4.*a3*g11*h21*h21*g21*beta*beta*cth*sth
	+4.*a3*g22_2*h22*h22*alpha*alpha*s2
	+4.*a3*g12_2*h22*h22*alpha*alpha*c2
	+4.*a3*g21_2*h21*h21*alpha*alpha*s2
	-4.*a3*g21*h21*g22*h22*beta*beta*c2
	-2.*a3*g12_2*h22*h22*beta*beta*s2
	+8.*a3*g12*h22*h22*g22*alpha*alpha*sth*cth
	-2.*a3*g21_2*h21*h21*beta*beta*c2
	+4.*a3*g11*h21*g22*h22*beta*beta*cth*sth
	+8.*a3*g21*h21*g22*h22*alpha*alpha*s2
	+4.*a3*g21*h21*g12*h22*beta*beta*cth*sth
	+4.*a3*g11_2*h21*h21*alpha*alpha*c2
	+8.*a3*g11*h21*g22*h22*alpha*alpha*sth*cth)
	-(-2.*beta*beta*c2
	+4.*alpha*alpha*s2)*(-4.*a3*g21*h21*g12*h22*beta*c2*alpha
	-4.*a3*g11*h21*h21*g21*beta*c2*alpha
	+4.*a3*g11_2*h21*h21*alpha*sth*beta*cth
	-4.*a3*g21_2*h21*h21*alpha*sth*beta*cth
	+4.*a3*g12_2*h22*h22*alpha*sth*beta*cth
	+4.*a3*g11*h21*g22*h22*alpha*s2*beta
	-4.*a3*g12*h22*h22*g22*beta*c2*alpha
	+4.*a3*g12*h22*h22*g22*alpha*s2*beta
	-4.*a3*g11*h21*g22*h22*beta*c2*alpha
	-8.*a3*g21*h21*g22*h22*alpha*sth*beta*cth
	-4.*a3*g22_2*h22*h22*alpha*sth*beta*cth
	+8.*a3*g11*h21*g12*h22*alpha*sth*beta*cth
	+4.*a3*g21*h21*g12*h22*alpha*s2*beta
	+4.*a3*g11*h21*h21*g21*alpha*s2*beta)
	-4.*beta*cth*alpha*sth*(-a6*a6*h21*h21
	-a6*a6*h22*h22
	+a3*g12_2*h22*h22*beta*beta*s2
	+2.*a3*g11*h21*g12*h22*beta*beta*s2
	+a3*g11_2*h21*h21*beta*beta*s2
	-2.*a3*g11*h21*h21*g21*beta*beta*cth*sth
	+a3*g21_2*h21*h21*beta*beta*c2
	-2.*a3*g12*h22*h22*g22*beta*beta*cth*sth
	-2.*a3*g21*h21*g12*h22*beta*beta*cth*sth
	+2.*a3*g21*h21*g22*h22*beta*beta*c2
	+a3*g22_2*h22*h22*beta*beta*c2
	-2.*a3*g11*h21*g22*h22*beta*beta*cth*sth);

	G4 = -2.*beta*beta*c2*(-a6*a6*h21*h21
	-a6*a6*h22*h22
	+a3*g12_2*h22*h22*beta*beta*s2
	+2.*a3*g11*h21*g12*h22*beta*beta*s2
	+a3*g11_2*h21*h21*beta*beta*s2
	-2.*a3*g11*h21*h21*g21*beta*beta*cth*sth
	+a3*g21_2*h21*h21*beta*beta*c2
	-2.*a3*g12*h22*h22*g22*beta*beta*cth*sth
	-2.*a3*g21*h21*g12*h22*beta*beta*cth*sth
	+2.*a3*g21*h21*g22*h22*beta*beta*c2
	+a3*g22_2*h22*h22*beta*beta*c2
	-2.*a3*g11*h21*g22*h22*beta*beta*cth*sth)
	+4.*beta*cth*alpha*sth*(4.*a3*g22_2*h22*h22*alpha*sth*beta*cth
	+4.*a3*g21_2*h21*h21*alpha*sth*beta*cth
	+4.*a3*g11*h21*h21*g21*beta*c2*alpha
	-4.*a3*g12*h22*h22*g22*alpha*s2*beta
	+4.*a3*g11*h21*g22*h22*beta*c2*alpha
	-4.*a3*g11_2*h21*h21*alpha*sth*beta*cth
	+8.*a3*g21*h21*g22*h22*alpha*sth*beta*cth
	-4.*a3*g12_2*h22*h22*alpha*sth*beta*cth
	+4.*a3*g12*h22*h22*g22*beta*c2*alpha
	-4.*a3*g11*h21*g22*h22*alpha*s2*beta
	-8.*a3*g11*h21*g12*h22*alpha*sth*beta*cth
	-4.*a3*g11*h21*h21*g21*alpha*s2*beta
	-4.*a3*g21*h21*g12*h22*alpha*s2*beta
	+4.*a3*g21*h21*g12*h22*beta*c2*alpha)
	-(-2.*beta*beta*c2
	+4.*alpha*alpha*s2)*(-4.*a3*g11*h21*g12*h22*beta*beta*s2
	+8.*a3*g11*h21*g12*h22*alpha*alpha*c2
	-2.*a3*g22_2*h22*h22*beta*beta*c2
	+4.*a3*g12*h22*h22*g22*beta*beta*cth*sth
	-2.*a6*a6*h21*h21
	+8.*a3*g21*h21*g12*h22*alpha*alpha*sth*cth
	-2.*a3*g11_2*h21*h21*beta*beta*s2
	-2.*a6*a6*h22*h22
	+8.*a3*g11*h21*h21*g21*alpha*alpha*sth*cth
	+4.*a3*g11*h21*h21*g21*beta*beta*cth*sth
	+4.*a3*g22_2*h22*h22*alpha*alpha*s2
	+4.*a3*g12_2*h22*h22*alpha*alpha*c2
	+4.*a3*g21_2*h21*h21*alpha*alpha*s2
	-4.*a3*g21*h21*g22*h22*beta*beta*c2
	-2.*a3*g12_2*h22*h22*beta*beta*s2
	+8.*a3*g12*h22*h22*g22*alpha*alpha*sth*cth
	-2.*a3*g21_2*h21*h21*beta*beta*c2
	+4.*a3*g11*h21*g22*h22*beta*beta*cth*sth
	+8.*a3*g21*h21*g22*h22*alpha*alpha*s2
	+4.*a3*g21*h21*g12*h22*beta*beta*cth*sth
	+4.*a3*g11_2*h21*h21*alpha*alpha*c2
	+8.*a3*g11*h21*g22*h22*alpha*alpha*sth*cth)
	-4.*beta*cth*alpha*sth*(-4.*a3*g21*h21*g12*h22*beta*c2*alpha
	-4.*a3*g11*h21*h21*g21*beta*c2*alpha
	+4.*a3*g11_2*h21*h21*alpha*sth*beta*cth
	-4.*a3*g21_2*h21*h21*alpha*sth*beta*cth
	+4.*a3*g12_2*h22*h22*alpha*sth*beta*cth
	+4.*a3*g11*h21*g22*h22*alpha*s2*beta
	-4.*a3*g12*h22*h22*g22*beta*c2*alpha
	+4.*a3*g12*h22*h22*g22*alpha*s2*beta
	-4.*a3*g11*h21*g22*h22*beta*c2*alpha
	-8.*a3*g21*h21*g22*h22*alpha*sth*beta*cth
	-4.*a3*g22_2*h22*h22*alpha*sth*beta*cth
	+8.*a3*g11*h21*g12*h22*alpha*sth*beta*cth
	+4.*a3*g21*h21*g12*h22*alpha*s2*beta
	+4.*a3*g11*h21*h21*g21*alpha*s2*beta);

	H0 = -beta*beta*sth*cth*(a5*g11_2*h11*h21*beta*beta*s2
	-a4*a6*h11*h21
	-a4*a6*h12*h22
	-a5*g22*h12*g11*h21*beta*beta*sth*cth
	+a5*g22*h12*g21*h21*beta*beta*c2
	+a5*g22_2*h12*h22*beta*beta*c2
	-2.*a5*g11*h11*g21*h21*beta*beta*sth*cth
	+a5*g11*h11*g12*h22*beta*beta*s2
	-a5*g11*h11*g22*h22*beta*beta*sth*cth
	+a5*g21_2*h11*h21*beta*beta*c2
	-a5*g21*h11*g12*h22*beta*beta*sth*cth
	+a5*g12*h12*g11*h21*beta*beta*s2
	+a5*g21*h11*g22*h22*beta*beta*c2
	-2.*a5*g12*h12*g22*h22*beta*beta*sth*cth
	-a5*g12*h12*g21*h21*beta*beta*sth*cth
	+a5*g12_2*h12*h22*beta*beta*s2);

	H1 = -beta*beta*sth*cth*(4.*a5*g11_2*h11*h21*alpha*cth*beta*sth
	+4.*a5*g11*h11*g12*h22*alpha*cth*beta*sth
	+2.*a5*g21*h11*g12*h22*beta*s2*alpha
	+4.*a5*g12*h12*g11*h21*alpha*cth*beta*sth
	+4.*a5*g12_2*h12*h22*alpha*cth*beta*sth
	+4.*a5*g12*h12*g22*h22*beta*s2*alpha
	-4.*a5*g22*h12*g21*h21*alpha*cth*beta*sth
	+2.*a5*g11*h11*g22*h22*beta*s2*alpha
	-4.*a5*g21_2*h11*h21*alpha*cth*beta*sth
	+2.*a5*g22*h12*g11*h21*beta*s2*alpha
	-4.*a5*g21*h11*g22*h22*alpha*cth*beta*sth
	-2.*a5*g11*h11*g22*h22*alpha*c2*beta
	-4.*a5*g11*h11*g21*h21*alpha*c2*beta
	+2.*a5*g12*h12*g21*h21*beta*s2*alpha
	-2.*a5*g12*h12*g21*h21*alpha*c2*beta
	-2.*a5*g21*h11*g12*h22*alpha*c2*beta
	+4.*a5*g11*h11*g21*h21*beta*s2*alpha
	-2.*a5*g22*h12*g11*h21*alpha*c2*beta
	-4.*a5*g22_2*h12*h22*alpha*cth*beta*sth
	-4.*a5*g12*h12*g22*h22*alpha*c2*beta)
	-(-2.*beta*s2*alpha
	+2.*alpha*c2*beta)*(a5*g11_2*h11*h21*beta*beta*s2
	-a4*a6*h11*h21
	-a4*a6*h12*h22
	-a5*g22*h12*g11*h21*beta*beta*sth*cth
	+a5*g22*h12*g21*h21*beta*beta*c2
	+a5*g22_2*h12*h22*beta*beta*c2
	-2.*a5*g11*h11*g21*h21*beta*beta*sth*cth
	+a5*g11*h11*g12*h22*beta*beta*s2
	-a5*g11*h11*g22*h22*beta*beta*sth*cth
	+a5*g21_2*h11*h21*beta*beta*c2
	-a5*g21*h11*g12*h22*beta*beta*sth*cth
	+a5*g12*h12*g11*h21*beta*beta*s2
	+a5*g21*h11*g22*h22*beta*beta*c2
	-2.*a5*g12*h12*g22*h22*beta*beta*sth*cth
	-a5*g12*h12*g21*h21*beta*beta*sth*cth
	+a5*g12_2*h12*h22*beta*beta*s2);

	H2 = -beta*beta*sth*cth*(-2.*a5*g11_2*h11*h21*beta*beta*s2
	-2.*a4*a6*h11*h21
	-2.*a4*a6*h12*h22
	+2.*a5*g22*h12*g11*h21*beta*beta*sth*cth
	-2.*a5*g22*h12*g21*h21*beta*beta*c2
	-2.*a5*g22_2*h12*h22*beta*beta*c2
	+4.*a5*g11*h11*g21*h21*beta*beta*sth*cth
	-2.*a5*g11*h11*g12*h22*beta*beta*s2
	+2.*a5*g11*h11*g22*h22*beta*beta*sth*cth
	-2.*a5*g21_2*h11*h21*beta*beta*c2
	+2.*a5*g21*h11*g12*h22*beta*beta*sth*cth
	-2.*a5*g12*h12*g11*h21*beta*beta*s2
	-2.*a5*g21*h11*g22*h22*beta*beta*c2
	+4.*a5*g12*h12*g22*h22*beta*beta*sth*cth
	+2.*a5*g12*h12*g21*h21*beta*beta*sth*cth
	-2.*a5*g12_2*h12*h22*beta*beta*s2
	+4.*a5*g11*h11*g22*h22*alpha*alpha*cth*sth
	+4.*a5*g11*h11*g12*h22*alpha*alpha*c2
	+4.*a5*g21_2*h11*h21*alpha*alpha*s2
	+8.*a5*g11*h11*g21*h21*alpha*alpha*cth*sth
	+8.*a5*g12*h12*g22*h22*alpha*alpha*cth*sth
	+4.*a5*g22*h12*g21*h21*alpha*alpha*s2
	+4.*a5*g22_2*h12*h22*alpha*alpha*s2
	+4.*a5*g22*h12*g11*h21*alpha*alpha*cth*sth
	+4.*a5*g12_2*h12*h22*alpha*alpha*c2
	+4.*a5*g11_2*h11*h21*alpha*alpha*c2
	+4.*a5*g12*h12*g21*h21*alpha*alpha*cth*sth
	+4.*a5*g21*h11*g22*h22*alpha*alpha*s2
	+4.*a5*g12*h12*g11*h21*alpha*alpha*c2
	+4.*a5*g21*h11*g12*h22*alpha*alpha*cth*sth)
	-(-2.*beta*s2*alpha
	+2.*alpha*c2*beta)*(4.*a5*g11_2*h11*h21*alpha*cth*beta*sth
	+4.*a5*g11*h11*g12*h22*alpha*cth*beta*sth
	+2.*a5*g21*h11*g12*h22*beta*s2*alpha
	+4.*a5*g12*h12*g11*h21*alpha*cth*beta*sth
	+4.*a5*g12_2*h12*h22*alpha*cth*beta*sth
	+4.*a5*g12*h12*g22*h22*beta*s2*alpha
	-4.*a5*g22*h12*g21*h21*alpha*cth*beta*sth
	+2.*a5*g11*h11*g22*h22*beta*s2*alpha
	-4.*a5*g21_2*h11*h21*alpha*cth*beta*sth
	+2.*a5*g22*h12*g11*h21*beta*s2*alpha
	-4.*a5*g21*h11*g22*h22*alpha*cth*beta*sth
	-2.*a5*g11*h11*g22*h22*alpha*c2*beta
	-4.*a5*g11*h11*g21*h21*alpha*c2*beta
	+2.*a5*g12*h12*g21*h21*beta*s2*alpha
	-2.*a5*g12*h12*g21*h21*alpha*c2*beta
	-2.*a5*g21*h11*g12*h22*alpha*c2*beta
	+4.*a5*g11*h11*g21*h21*beta*s2*alpha
	-2.*a5*g22*h12*g11*h21*alpha*c2*beta
	-4.*a5*g22_2*h12*h22*alpha*cth*beta*sth
	-4.*a5*g12*h12*g22*h22*alpha*c2*beta)
	-(-2.*beta*beta*sth*cth
	-4.*alpha*alpha*cth*sth)*(a5*g11_2*h11*h21*beta*beta*s2
	-a4*a6*h11*h21
	-a4*a6*h12*h22
	-a5*g22*h12*g11*h21*beta*beta*sth*cth
	+a5*g22*h12*g21*h21*beta*beta*c2
	+a5*g22_2*h12*h22*beta*beta*c2
	-2.*a5*g11*h11*g21*h21*beta*beta*sth*cth
	+a5*g11*h11*g12*h22*beta*beta*s2
	-a5*g11*h11*g22*h22*beta*beta*sth*cth
	+a5*g21_2*h11*h21*beta*beta*c2
	-a5*g21*h11*g12*h22*beta*beta*sth*cth
	+a5*g12*h12*g11*h21*beta*beta*s2
	+a5*g21*h11*g22*h22*beta*beta*c2
	-2.*a5*g12*h12*g22*h22*beta*beta*sth*cth
	-a5*g12*h12*g21*h21*beta*beta*sth*cth
	+a5*g12_2*h12*h22*beta*beta*s2);

	H3 = -beta*beta*sth*cth*(4.*a5*g21_2*h11*h21*alpha*cth*beta*sth
	-4.*a5*g11*h11*g12*h22*alpha*cth*beta*sth
	-4.*a5*g12*h12*g22*h22*beta*s2*alpha
	-4.*a5*g12*h12*g11*h21*alpha*cth*beta*sth
	-2.*a5*g11*h11*g22*h22*beta*s2*alpha
	-2.*a5*g12*h12*g21*h21*beta*s2*alpha
	+4.*a5*g12*h12*g22*h22*alpha*c2*beta
	+2.*a5*g11*h11*g22*h22*alpha*c2*beta
	+2.*a5*g21*h11*g12*h22*alpha*c2*beta
	+2.*a5*g12*h12*g21*h21*alpha*c2*beta
	+4.*a5*g21*h11*g22*h22*alpha*cth*beta*sth
	-2.*a5*g22*h12*g11*h21*beta*s2*alpha
	+4.*a5*g22*h12*g21*h21*alpha*cth*beta*sth
	+4.*a5*g22_2*h12*h22*alpha*cth*beta*sth
	-4.*a5*g11_2*h11*h21*alpha*cth*beta*sth
	+4.*a5*g11*h11*g21*h21*alpha*c2*beta
	-4.*a5*g11*h11*g21*h21*beta*s2*alpha
	-2.*a5*g21*h11*g12*h22*beta*s2*alpha
	-4.*a5*g12_2*h12*h22*alpha*cth*beta*sth
	+2.*a5*g22*h12*g11*h21*alpha*c2*beta)
	-(-2.*beta*s2*alpha
	+2.*alpha*c2*beta)*(-2.*a5*g11_2*h11*h21*beta*beta*s2
	-2.*a4*a6*h11*h21
	-2.*a4*a6*h12*h22
	+2.*a5*g22*h12*g11*h21*beta*beta*sth*cth
	-2.*a5*g22*h12*g21*h21*beta*beta*c2
	-2.*a5*g22_2*h12*h22*beta*beta*c2
	+4.*a5*g11*h11*g21*h21*beta*beta*sth*cth
	-2.*a5*g11*h11*g12*h22*beta*beta*s2
	+2.*a5*g11*h11*g22*h22*beta*beta*sth*cth
	-2.*a5*g21_2*h11*h21*beta*beta*c2
	+2.*a5*g21*h11*g12*h22*beta*beta*sth*cth
	-2.*a5*g12*h12*g11*h21*beta*beta*s2
	-2.*a5*g21*h11*g22*h22*beta*beta*c2
	+4.*a5*g12*h12*g22*h22*beta*beta*sth*cth
	+2.*a5*g12*h12*g21*h21*beta*beta*sth*cth
	-2.*a5*g12_2*h12*h22*beta*beta*s2
	+4.*a5*g11*h11*g22*h22*alpha*alpha*cth*sth
	+4.*a5*g11*h11*g12*h22*alpha*alpha*c2
	+4.*a5*g21_2*h11*h21*alpha*alpha*s2
	+8.*a5*g11*h11*g21*h21*alpha*alpha*cth*sth
	+8.*a5*g12*h12*g22*h22*alpha*alpha*cth*sth
	+4.*a5*g22*h12*g21*h21*alpha*alpha*s2
	+4.*a5*g22_2*h12*h22*alpha*alpha*s2
	+4.*a5*g22*h12*g11*h21*alpha*alpha*cth*sth
	+4.*a5*g12_2*h12*h22*alpha*alpha*c2
	+4.*a5*g11_2*h11*h21*alpha*alpha*c2
	+4.*a5*g12*h12*g21*h21*alpha*alpha*cth*sth
	+4.*a5*g21*h11*g22*h22*alpha*alpha*s2
	+4.*a5*g12*h12*g11*h21*alpha*alpha*c2
	+4.*a5*g21*h11*g12*h22*alpha*alpha*cth*sth)
	-(-2.*beta*beta*sth*cth
	-4.*alpha*alpha*cth*sth)*(4.*a5*g11_2*h11*h21*alpha*cth*beta*sth
	+4.*a5*g11*h11*g12*h22*alpha*cth*beta*sth
	+2.*a5*g21*h11*g12*h22*beta*s2*alpha
	+4.*a5*g12*h12*g11*h21*alpha*cth*beta*sth
	+4.*a5*g12_2*h12*h22*alpha*cth*beta*sth
	+4.*a5*g12*h12*g22*h22*beta*s2*alpha
	-4.*a5*g22*h12*g21*h21*alpha*cth*beta*sth
	+2.*a5*g11*h11*g22*h22*beta*s2*alpha
	-4.*a5*g21_2*h11*h21*alpha*cth*beta*sth
	+2.*a5*g22*h12*g11*h21*beta*s2*alpha
	-4.*a5*g21*h11*g22*h22*alpha*cth*beta*sth
	-2.*a5*g11*h11*g22*h22*alpha*c2*beta
	-4.*a5*g11*h11*g21*h21*alpha*c2*beta
	+2.*a5*g12*h12*g21*h21*beta*s2*alpha
	-2.*a5*g12*h12*g21*h21*alpha*c2*beta
	-2.*a5*g21*h11*g12*h22*alpha*c2*beta
	+4.*a5*g11*h11*g21*h21*beta*s2*alpha
	-2.*a5*g22*h12*g11*h21*alpha*c2*beta
	-4.*a5*g22_2*h12*h22*alpha*cth*beta*sth
	-4.*a5*g12*h12*g22*h22*alpha*c2*beta)
	-(-2.*alpha*c2*beta
	+2.*beta*s2*alpha)*(a5*g11_2*h11*h21*beta*beta*s2
	-a4*a6*h11*h21
	-a4*a6*h12*h22
	-a5*g22*h12*g11*h21*beta*beta*sth*cth
	+a5*g22*h12*g21*h21*beta*beta*c2
	+a5*g22_2*h12*h22*beta*beta*c2
	-2.*a5*g11*h11*g21*h21*beta*beta*sth*cth
	+a5*g11*h11*g12*h22*beta*beta*s2
	-a5*g11*h11*g22*h22*beta*beta*sth*cth
	+a5*g21_2*h11*h21*beta*beta*c2
	-a5*g21*h11*g12*h22*beta*beta*sth*cth
	+a5*g12*h12*g11*h21*beta*beta*s2
	+a5*g21*h11*g22*h22*beta*beta*c2
	-2.*a5*g12*h12*g22*h22*beta*beta*sth*cth
	-a5*g12*h12*g21*h21*beta*beta*sth*cth
	+a5*g12_2*h12*h22*beta*beta*s2);

	H4 = -2.*beta*beta*sth*cth*(a5*g11_2*h11*h21*beta*beta*s2
	-a4*a6*h11*h21
	-a4*a6*h12*h22
	-a5*g22*h12*g11*h21*beta*beta*sth*cth
	+a5*g22*h12*g21*h21*beta*beta*c2
	+a5*g22_2*h12*h22*beta*beta*c2
	-2.*a5*g11*h11*g21*h21*beta*beta*sth*cth
	+a5*g11*h11*g12*h22*beta*beta*s2
	-a5*g11*h11*g22*h22*beta*beta*sth*cth
	+a5*g21_2*h11*h21*beta*beta*c2
	-a5*g21*h11*g12*h22*beta*beta*sth*cth
	+a5*g12*h12*g11*h21*beta*beta*s2
	+a5*g21*h11*g22*h22*beta*beta*c2
	-2.*a5*g12*h12*g22*h22*beta*beta*sth*cth
	-a5*g12*h12*g21*h21*beta*beta*sth*cth
	+a5*g12_2*h12*h22*beta*beta*s2)
	-(-2.*beta*s2*alpha
	+2.*alpha*c2*beta)*(4.*a5*g21_2*h11*h21*alpha*cth*beta*sth
	-4.*a5*g11*h11*g12*h22*alpha*cth*beta*sth
	-4.*a5*g12*h12*g22*h22*beta*s2*alpha
	-4.*a5*g12*h12*g11*h21*alpha*cth*beta*sth
	-2.*a5*g11*h11*g22*h22*beta*s2*alpha
	-2.*a5*g12*h12*g21*h21*beta*s2*alpha
	+4.*a5*g12*h12*g22*h22*alpha*c2*beta
	+2.*a5*g11*h11*g22*h22*alpha*c2*beta
	+2.*a5*g21*h11*g12*h22*alpha*c2*beta
	+2.*a5*g12*h12*g21*h21*alpha*c2*beta
	+4.*a5*g21*h11*g22*h22*alpha*cth*beta*sth
	-2.*a5*g22*h12*g11*h21*beta*s2*alpha
	+4.*a5*g22*h12*g21*h21*alpha*cth*beta*sth
	+4.*a5*g22_2*h12*h22*alpha*cth*beta*sth
	-4.*a5*g11_2*h11*h21*alpha*cth*beta*sth
	+4.*a5*g11*h11*g21*h21*alpha*c2*beta
	-4.*a5*g11*h11*g21*h21*beta*s2*alpha
	-2.*a5*g21*h11*g12*h22*beta*s2*alpha
	-4.*a5*g12_2*h12*h22*alpha*cth*beta*sth
	+2.*a5*g22*h12*g11*h21*alpha*c2*beta)
	-(-2.*beta*beta*sth*cth
	-4.*alpha*alpha*cth*sth)*(-2.*a5*g11_2*h11*h21*beta*beta*s2
	-2.*a4*a6*h11*h21
	-2.*a4*a6*h12*h22
	+2.*a5*g22*h12*g11*h21*beta*beta*sth*cth
	-2.*a5*g22*h12*g21*h21*beta*beta*c2
	-2.*a5*g22_2*h12*h22*beta*beta*c2
	+4.*a5*g11*h11*g21*h21*beta*beta*sth*cth
	-2.*a5*g11*h11*g12*h22*beta*beta*s2
	+2.*a5*g11*h11*g22*h22*beta*beta*sth*cth
	-2.*a5*g21_2*h11*h21*beta*beta*c2
	+2.*a5*g21*h11*g12*h22*beta*beta*sth*cth
	-2.*a5*g12*h12*g11*h21*beta*beta*s2
	-2.*a5*g21*h11*g22*h22*beta*beta*c2
	+4.*a5*g12*h12*g22*h22*beta*beta*sth*cth
	+2.*a5*g12*h12*g21*h21*beta*beta*sth*cth
	-2.*a5*g12_2*h12*h22*beta*beta*s2
	+4.*a5*g11*h11*g22*h22*alpha*alpha*cth*sth
	+4.*a5*g11*h11*g12*h22*alpha*alpha*c2
	+4.*a5*g21_2*h11*h21*alpha*alpha*s2
	+8.*a5*g11*h11*g21*h21*alpha*alpha*cth*sth
	+8.*a5*g12*h12*g22*h22*alpha*alpha*cth*sth
	+4.*a5*g22*h12*g21*h21*alpha*alpha*s2
	+4.*a5*g22_2*h12*h22*alpha*alpha*s2
	+4.*a5*g22*h12*g11*h21*alpha*alpha*cth*sth
	+4.*a5*g12_2*h12*h22*alpha*alpha*c2
	+4.*a5*g11_2*h11*h21*alpha*alpha*c2
	+4.*a5*g12*h12*g21*h21*alpha*alpha*cth*sth
	+4.*a5*g21*h11*g22*h22*alpha*alpha*s2
	+4.*a5*g12*h12*g11*h21*alpha*alpha*c2
	+4.*a5*g21*h11*g12*h22*alpha*alpha*cth*sth)
	-(-2.*alpha*c2*beta
	+2.*beta*s2*alpha)*(4.*a5*g11_2*h11*h21*alpha*cth*beta*sth
	+4.*a5*g11*h11*g12*h22*alpha*cth*beta*sth
	+2.*a5*g21*h11*g12*h22*beta*s2*alpha
	+4.*a5*g12*h12*g11*h21*alpha*cth*beta*sth
	+4.*a5*g12_2*h12*h22*alpha*cth*beta*sth
	+4.*a5*g12*h12*g22*h22*beta*s2*alpha
	-4.*a5*g22*h12*g21*h21*alpha*cth*beta*sth
	+2.*a5*g11*h11*g22*h22*beta*s2*alpha
	-4.*a5*g21_2*h11*h21*alpha*cth*beta*sth
	+2.*a5*g22*h12*g11*h21*beta*s2*alpha
	-4.*a5*g21*h11*g22*h22*alpha*cth*beta*sth
	-2.*a5*g11*h11*g22*h22*alpha*c2*beta
	-4.*a5*g11*h11*g21*h21*alpha*c2*beta
	+2.*a5*g12*h12*g21*h21*beta*s2*alpha
	-2.*a5*g12*h12*g21*h21*alpha*c2*beta
	-2.*a5*g21*h11*g12*h22*alpha*c2*beta
	+4.*a5*g11*h11*g21*h21*beta*s2*alpha
	-2.*a5*g22*h12*g11*h21*alpha*c2*beta
	-4.*a5*g22_2*h12*h22*alpha*cth*beta*sth
	-4.*a5*g12*h12*g22*h22*alpha*c2*beta);

	J0 = -beta*cth*(-a4*a6*g12*h22
	-a4*a6*g11*h21
	+a5*g12*h22*beta*beta*s2
	+a5*g22*g12*g21*h21*beta*beta*c2
	+a5*g12_2*g11*h21*beta*beta*s2
	-a5*g12_2*g21*h21*beta*beta*sth*cth
	-a5*g22*g12*g11*h21*beta*beta*sth*cth
	-a5*g22*h22*beta*beta*sth*cth
	+a5*g11*h21*beta*beta*s2
	-2.*a5*g11_2*g21*h21*beta*beta*sth*cth
	-a5*g12*h22*beta*beta*sth*cth
	-a5*g11*h21*beta*beta*sth*cth
	-a5*g21*g11*g12*h22*beta*beta*sth*cth
	+a5*g12_3*h22*beta*beta*s2
	+a5*g22_2*g12*h22*beta*beta*c2
	-a5*g21*h21*beta*beta*sth*cth
	+a5*g21*h21*beta*beta*c2
	+a5*g11_3*h21*beta*beta*s2
	+a5*g11_2*g12*h22*beta*beta*s2
	-a5*g11_2*g22*h22*beta*beta*sth*cth
	+a5*g21_2*g11*h21*beta*beta*c2
	+a5*g21*g11*g22*h22*beta*beta*c2
	+a5*g22*h22*beta*beta*c2
	-2.*a5*g12_2*g22*h22*beta*beta*sth*cth);

	J1 = -beta*cth*(4.*a5*g11_2*g12*h22*alpha*cth*beta*sth
	-2.*a5*g11_2*g22*h22*alpha*c2*beta
	-2.*a5*g21*g11*g12*h22*alpha*c2*beta
	+2.*a5*g21*g11*g12*h22*beta*s2*alpha
	-4.*a5*g11_2*g21*h21*alpha*c2*beta
	+4.*a5*g11_2*g21*h21*beta*s2*alpha
	-4.*a5*g22*g12*g21*h21*alpha*cth*beta*sth
	+4.*a5*g12_3*h22*alpha*cth*beta*sth
	+2.*a5*g12*h22*beta*s2*alpha
	-4.*a5*g21*h21*alpha*cth*beta*sth
	+4.*a5*g12*h22*alpha*cth*beta*sth
	+4.*a5*g11_3*h21*alpha*cth*beta*sth
	-2.*a5*g12*h22*alpha*c2*beta
	+2.*a5*g12_2*g21*h21*beta*s2*alpha
	-2.*a5*g22*g12*g11*h21*alpha*c2*beta
	+2.*a5*g22*g12*g11*h21*beta*s2*alpha
	+4.*a5*g11*h21*alpha*cth*beta*sth
	-2.*a5*g11*h21*alpha*c2*beta
	-2.*a5*g22*h22*alpha*c2*beta
	+4.*a5*g12_2*g11*h21*alpha*cth*beta*sth
	-2.*a5*g12_2*g21*h21*alpha*c2*beta
	+2.*a5*g22*h22*beta*s2*alpha
	-4.*a5*g22_2*g12*h22*alpha*cth*beta*sth
	-2.*a5*g21*h21*alpha*c2*beta
	+2.*a5*g21*h21*beta*s2*alpha
	-4.*a5*g21*g11*g22*h22*alpha*cth*beta*sth
	-4.*a5*g22*h22*alpha*cth*beta*sth
	-4.*a5*g12_2*g22*h22*alpha*c2*beta
	+4.*a5*g12_2*g22*h22*beta*s2*alpha
	+2.*a5*g11*h21*beta*s2*alpha
	+2.*a5*g11_2*g22*h22*beta*s2*alpha
	-4.*a5*g21_2*g11*h21*alpha*cth*beta*sth)
	+2.*alpha*sth*(-a4*a6*g12*h22
	-a4*a6*g11*h21
	+a5*g12*h22*beta*beta*s2
	+a5*g22*g12*g21*h21*beta*beta*c2
	+a5*g12_2*g11*h21*beta*beta*s2
	-a5*g12_2*g21*h21*beta*beta*sth*cth
	-a5*g22*g12*g11*h21*beta*beta*sth*cth
	-a5*g22*h22*beta*beta*sth*cth
	+a5*g11*h21*beta*beta*s2
	-2.*a5*g11_2*g21*h21*beta*beta*sth*cth
	-a5*g12*h22*beta*beta*sth*cth
	-a5*g11*h21*beta*beta*sth*cth
	-a5*g21*g11*g12*h22*beta*beta*sth*cth
	+a5*g12_3*h22*beta*beta*s2
	+a5*g22_2*g12*h22*beta*beta*c2
	-a5*g21*h21*beta*beta*sth*cth
	+a5*g21*h21*beta*beta*c2
	+a5*g11_3*h21*beta*beta*s2
	+a5*g11_2*g12*h22*beta*beta*s2
	-a5*g11_2*g22*h22*beta*beta*sth*cth
	+a5*g21_2*g11*h21*beta*beta*c2
	+a5*g21*g11*g22*h22*beta*beta*c2
	+a5*g22*h22*beta*beta*c2
	-2.*a5*g12_2*g22*h22*beta*beta*sth*cth);

	J2 = -beta*cth*(-(2*a4*a6*g12*h22)
	-(2*a4*a6*g11*h21)
	-2.*a5*g12*h22*beta*beta*s2
	-2.*a5*g22*g12*g21*h21*beta*beta*c2
	-2.*a5*(g12_2)*g11*h21*beta*beta*s2
	+2.*a5*(g12_2)*g21*h21*beta*beta*sth*cth
	+2.*a5*g22*g12*g11*h21*beta*beta*sth*cth
	+2.*a5*g22*h22*beta*beta*sth*cth
	-2.*a5*g11*h21*beta*beta*s2
	+4.*a5*(g11_2)*g21*h21*beta*beta*sth*cth
	+2.*a5*g12*h22*beta*beta*sth*cth
	+2.*a5*g11*h21*beta*beta*sth*cth
	+2.*a5*g21*g11*g12*h22*beta*beta*sth*cth
	-2.*a5*g12_3*h22*beta*beta*s2
	-2.*a5*g22_2*g12*h22*beta*beta*c2
	+2.*a5*g21*h21*beta*beta*sth*cth
	-2.*a5*g21*h21*beta*beta*c2
	-2.*a5*g11_3*h21*beta*beta*s2
	-2.*a5*(g11_2)*g12*h22*beta*beta*s2
	+2.*a5*(g11_2)*g22*h22*beta*beta*sth*cth
	-2.*a5*g21_2*g11*h21*beta*beta*c2
	-2.*a5*g21*g11*g22*h22*beta*beta*c2
	-2.*a5*g22*h22*beta*beta*c2
	+4.*a5*(g12_2)*g22*h22*beta*beta*sth*cth
	+4.*a5*(g11_2)*g12*h22*alpha*alpha*c2
	+4.*a5*(g11_2)*g22*h22*alpha*alpha*cth*sth
	+4.*a5*g21*g11*g12*h22*alpha*alpha*cth*sth
	+8.*a5*(g11_2)*g21*h21*alpha*alpha*cth*sth
	+4.*a5*g12_3*h22*alpha*alpha*c2
	+4.*a5*g22_2*g12*h22*alpha*alpha*s2
	+4.*a5*g11*h21*alpha*alpha*c2
	+4.*a5*g21*h21*alpha*alpha*s2
	+4.*a5*g12*h22*alpha*alpha*c2
	+4.*a5*g11_3*h21*alpha*alpha*c2
	+4.*a5*g12*h22*alpha*alpha*cth*sth
	+4.*a5*g22*g12*g21*h21*alpha*alpha*s2
	+4.*a5*g22*g12*g11*h21*alpha*alpha*cth*sth
	+4.*a5*g11*h21*alpha*alpha*cth*sth
	+4.*a5*g21*g11*g22*h22*alpha*alpha*s2
	+4.*a5*g22*h22*alpha*alpha*cth*sth
	+4.*a5*(g12_2)*g21*h21*alpha*alpha*cth*sth
	+4.*a5*g21*h21*alpha*alpha*cth*sth
	+4.*a5*(g12_2)*g11*h21*alpha*alpha*c2
	+4.*a5*g22*h22*alpha*alpha*s2
	+8.*a5*(g12_2)*g22*h22*alpha*alpha*cth*sth
	+4.*a5*g21_2*g11*h21*alpha*alpha*s2)
	+2.*alpha*sth*(4.*a5*(g11_2)*g12*h22*alpha*cth*beta*sth
	-2.*a5*(g11_2)*g22*h22*alpha*c2*beta
	-2.*a5*g21*g11*g12*h22*alpha*c2*beta
	+2.*a5*g21*g11*g12*h22*beta*s2*alpha
	-4.*a5*(g11_2)*g21*h21*alpha*c2*beta
	+4.*a5*(g11_2)*g21*h21*beta*s2*alpha
	-4.*a5*g22*g12*g21*h21*alpha*cth*beta*sth
	+4.*a5*g12_3*h22*alpha*cth*beta*sth
	+2.*a5*g12*h22*beta*s2*alpha
	-4.*a5*g21*h21*alpha*cth*beta*sth
	+4.*a5*g12*h22*alpha*cth*beta*sth
	+4.*a5*g11_3*h21*alpha*cth*beta*sth
	-2.*a5*g12*h22*alpha*c2*beta
	+2.*a5*(g12_2)*g21*h21*beta*s2*alpha
	-2.*a5*g22*g12*g11*h21*alpha*c2*beta
	+2.*a5*g22*g12*g11*h21*beta*s2*alpha
	+4.*a5*g11*h21*alpha*cth*beta*sth
	-2.*a5*g11*h21*alpha*c2*beta
	-2.*a5*g22*h22*alpha*c2*beta
	+4.*a5*(g12_2)*g11*h21*alpha*cth*beta*sth
	-2.*a5*(g12_2)*g21*h21*alpha*c2*beta
	+2.*a5*g22*h22*beta*s2*alpha
	-4.*a5*g22_2*g12*h22*alpha*cth*beta*sth
	-2.*a5*g21*h21*alpha*c2*beta
	+2.*a5*g21*h21*beta*s2*alpha
	-4.*a5*g21*g11*g22*h22*alpha*cth*beta*sth
	-4.*a5*g22*h22*alpha*cth*beta*sth
	-4.*a5*(g12_2)*g22*h22*alpha*c2*beta
	+4.*a5*(g12_2)*g22*h22*beta*s2*alpha
	+2.*a5*g11*h21*beta*s2*alpha
	+2.*a5*(g11_2)*g22*h22*beta*s2*alpha
	-4.*a5*g21_2*g11*h21*alpha*cth*beta*sth)
	+beta*cth*(-(a4*a6*g12*h22)
	-(a4*a6*g11*h21)
	+a5*g12*h22*beta*beta*s2
	+a5*g22*g12*g21*h21*beta*beta*c2
	+a5*(g12_2)*g11*h21*beta*beta*s2
	-a5*(g12_2)*g21*h21*beta*beta*sth*cth
	-a5*g22*g12*g11*h21*beta*beta*sth*cth
	-a5*g22*h22*beta*beta*sth*cth
	+a5*g11*h21*beta*beta*s2
	-2.*a5*(g11_2)*g21*h21*beta*beta*sth*cth
	-a5*g12*h22*beta*beta*sth*cth
	-a5*g11*h21*beta*beta*sth*cth
	-a5*g21*g11*g12*h22*beta*beta*sth*cth
	+a5*g12_3*h22*beta*beta*s2
	+a5*g22_2*g12*h22*beta*beta*c2
	-a5*g21*h21*beta*beta*sth*cth
	+a5*g21*h21*beta*beta*c2
	+a5*g11_3*h21*beta*beta*s2
	+a5*(g11_2)*g12*h22*beta*beta*s2
	-a5*(g11_2)*g22*h22*beta*beta*sth*cth
	+a5*g21_2*g11*h21*beta*beta*c2
	+a5*g21*g11*g22*h22*beta*beta*c2
	+a5*g22*h22*beta*beta*c2
	-2.*a5*(g12_2)*g22*h22*beta*beta*sth*cth);

	J3 = -beta*cth*(-4.*a5*g11_2*g12*h22*alpha*cth*beta*sth
	+2.*a5*g11_2*g22*h22*alpha*c2*beta
	+2.*a5*g21*g11*g12*h22*alpha*c2*beta
	-2.*a5*g21*g11*g12*h22*beta*s2*alpha
	+4.*a5*g11_2*g21*h21*alpha*c2*beta
	-4.*a5*g11_2*g21*h21*beta*s2*alpha
	+4.*a5*g22*g12*g21*h21*alpha*cth*beta*sth
	-4.*a5*g12_3*h22*alpha*cth*beta*sth
	-2.*a5*g12*h22*beta*s2*alpha
	+4.*a5*g21*h21*alpha*cth*beta*sth
	-4.*a5*g12*h22*alpha*cth*beta*sth
	-4.*a5*g11_3*h21*alpha*cth*beta*sth
	+2.*a5*g12*h22*alpha*c2*beta
	-2.*a5*g12_2*g21*h21*beta*s2*alpha
	+2.*a5*g22*g12*g11*h21*alpha*c2*beta
	-2.*a5*g22*g12*g11*h21*beta*s2*alpha
	-4.*a5*g11*h21*alpha*cth*beta*sth
	+2.*a5*g11*h21*alpha*c2*beta
	+2.*a5*g22*h22*alpha*c2*beta
	-4.*a5*g12_2*g11*h21*alpha*cth*beta*sth
	+2.*a5*g12_2*g21*h21*alpha*c2*beta
	-2.*a5*g22*h22*beta*s2*alpha
	+4.*a5*g22_2*g12*h22*alpha*cth*beta*sth
	+2.*a5*g21*h21*alpha*c2*beta
	-2.*a5*g21*h21*beta*s2*alpha
	+4.*a5*g21*g11*g22*h22*alpha*cth*beta*sth
	+4.*a5*g22*h22*alpha*cth*beta*sth
	+4.*a5*g12_2*g22*h22*alpha*c2*beta
	-4.*a5*g12_2*g22*h22*beta*s2*alpha
	-2.*a5*g11*h21*beta*s2*alpha
	-2.*a5*g11_2*g22*h22*beta*s2*alpha
	+4.*a5*g21_2*g11*h21*alpha*cth*beta*sth)
	+2.*alpha*sth*(-2.*a4*a6*g12*h22
	-2.*a4*a6*g11*h21
	-2.*a5*g12*h22*beta*beta*s2
	-2.*a5*g22*g12*g21*h21*beta*beta*c2
	-2.*a5*g12_2*g11*h21*beta*beta*s2
	+2.*a5*g12_2*g21*h21*beta*beta*sth*cth
	+2.*a5*g22*g12*g11*h21*beta*beta*sth*cth
	+2.*a5*g22*h22*beta*beta*sth*cth
	-2.*a5*g11*h21*beta*beta*s2
	+4.*a5*g11_2*g21*h21*beta*beta*sth*cth
	+2.*a5*g12*h22*beta*beta*sth*cth
	+2.*a5*g11*h21*beta*beta*sth*cth
	+2.*a5*g21*g11*g12*h22*beta*beta*sth*cth
	-2.*a5*g12_3*h22*beta*beta*s2
	-2.*a5*g22_2*g12*h22*beta*beta*c2
	+2.*a5*g21*h21*beta*beta*sth*cth
	-2.*a5*g21*h21*beta*beta*c2
	-2.*a5*g11_3*h21*beta*beta*s2
	-2.*a5*g11_2*g12*h22*beta*beta*s2
	+2.*a5*g11_2*g22*h22*beta*beta*sth*cth
	-2.*a5*g21_2*g11*h21*beta*beta*c2
	-2.*a5*g21*g11*g22*h22*beta*beta*c2
	-2.*a5*g22*h22*beta*beta*c2
	+4.*a5*g12_2*g22*h22*beta*beta*sth*cth
	+4.*a5*g11_2*g12*h22*alpha*alpha*c2
	+4.*a5*g11_2*g22*h22*alpha*alpha*cth*sth
	+4.*a5*g21*g11*g12*h22*alpha*alpha*cth*sth
	+8.*a5*g11_2*g21*h21*alpha*alpha*cth*sth
	+4.*a5*g12_3*h22*alpha*alpha*c2
	+4.*a5*g22_2*g12*h22*alpha*alpha*s2
	+4.*a5*g11*h21*alpha*alpha*c2
	+4.*a5*g21*h21*alpha*alpha*s2
	+4.*a5*g12*h22*alpha*alpha*c2
	+4.*a5*g11_3*h21*alpha*alpha*c2
	+4.*a5*g12*h22*alpha*alpha*cth*sth
	+4.*a5*g22*g12*g21*h21*alpha*alpha*s2
	+4.*a5*g22*g12*g11*h21*alpha*alpha*cth*sth
	+4.*a5*g11*h21*alpha*alpha*cth*sth
	+4.*a5*g21*g11*g22*h22*alpha*alpha*s2
	+4.*a5*g22*h22*alpha*alpha*cth*sth
	+4.*a5*g12_2*g21*h21*alpha*alpha*cth*sth
	+4.*a5*g21*h21*alpha*alpha*cth*sth
	+4.*a5*g12_2*g11*h21*alpha*alpha*c2
	+4.*a5*g22*h22*alpha*alpha*s2
	+8.*a5*g12_2*g22*h22*alpha*alpha*cth*sth
	+4.*a5*g21_2*g11*h21*alpha*alpha*s2)
	+beta*cth*(4.*a5*g11_2*g12*h22*alpha*cth*beta*sth
	-2.*a5*g11_2*g22*h22*alpha*c2*beta
	-2.*a5*g21*g11*g12*h22*alpha*c2*beta
	+2.*a5*g21*g11*g12*h22*beta*s2*alpha
	-4.*a5*g11_2*g21*h21*alpha*c2*beta
	+4.*a5*g11_2*g21*h21*beta*s2*alpha
	-4.*a5*g22*g12*g21*h21*alpha*cth*beta*sth
	+4.*a5*g12_3*h22*alpha*cth*beta*sth
	+2.*a5*g12*h22*beta*s2*alpha
	-4.*a5*g21*h21*alpha*cth*beta*sth
	+4.*a5*g12*h22*alpha*cth*beta*sth
	+4.*a5*g11_3*h21*alpha*cth*beta*sth
	-2.*a5*g12*h22*alpha*c2*beta
	+2.*a5*g12_2*g21*h21*beta*s2*alpha
	-2.*a5*g22*g12*g11*h21*alpha*c2*beta
	+2.*a5*g22*g12*g11*h21*beta*s2*alpha
	+4.*a5*g11*h21*alpha*cth*beta*sth
	-2.*a5*g11*h21*alpha*c2*beta
	-2.*a5*g22*h22*alpha*c2*beta
	+4.*a5*g12_2*g11*h21*alpha*cth*beta*sth
	-2.*a5*g12_2*g21*h21*alpha*c2*beta
	+2.*a5*g22*h22*beta*s2*alpha
	-4.*a5*g22_2*g12*h22*alpha*cth*beta*sth
	-2.*a5*g21*h21*alpha*c2*beta
	+2.*a5*g21*h21*beta*s2*alpha
	-4.*a5*g21*g11*g22*h22*alpha*cth*beta*sth
	-4.*a5*g22*h22*alpha*cth*beta*sth
	-4.*a5*g12_2*g22*h22*alpha*c2*beta
	+4.*a5*g12_2*g22*h22*beta*s2*alpha
	+2.*a5*g11*h21*beta*s2*alpha
	+2.*a5*g11_2*g22*h22*beta*s2*alpha
	-4.*a5*g21_2*g11*h21*alpha*cth*beta*sth);

	K0 = -beta*sth*(a5*g21_3*h11*beta*beta*c2
	-a5*g12*h12*beta*beta*sth*cth
	+a5*g11*h11*beta*beta*s2
	+a5*g21*h11*beta*beta*c2
	-a5*g22*h12*g11*g21*beta*beta*sth*cth
	-a5*g11*h11*beta*beta*sth*cth
	+a5*g12*h12*beta*beta*s2
	+a5*g22*h12*beta*beta*c2
	+a5*g22_3*h12*beta*beta*c2
	+a5*g11_2*h11*g21*beta*beta*s2
	-a5*g21*h11*g12*g22*beta*beta*sth*cth
	+a5*g22*h12*g21_2*beta*beta*c2
	+a5*g21*h11*g22_2*beta*beta*c2
	-a5*g21*h11*beta*beta*sth*cth
	+a5*g12*h12*g11*g21*beta*beta*s2
	-a5*g12*h12*g21_2*beta*beta*sth*cth
	-2.*a5*g11*h11*g21_2*beta*beta*sth*cth
	+a5*g11*h11*g12*g22*beta*beta*s2
	-a5*g11*h11*g22_2*beta*beta*sth*cth
	+a5*g12_2*h12*g22*beta*beta*s2
	-a5*g22*h12*beta*beta*sth*cth
	-2.*a5*g12*h12*g22_2*beta*beta*sth*cth
	-a4*a6*h11*g21
	-a4*a6*h12*g22);

	K1 = -beta*sth*(2.*a5*g22*h12*g11*g21*beta*s2*alpha
	+2.*a5*g12*h12*beta*s2*alpha
	-2.*a5*g22*h12*g11*g21*alpha*c2*beta
	-4.*a5*g22*h12*g21_2*alpha*cth*beta*sth
	+4.*a5*g12_2*h12*g22*alpha*cth*beta*sth
	-4.*a5*g21*h11*g22_2*alpha*cth*beta*sth
	+4.*a5*g12*h12*g11*g21*alpha*cth*beta*sth
	+4.*a5*g11_2*h11*g21*alpha*cth*beta*sth
	+4.*a5*g12*h12*alpha*cth*beta*sth
	-4.*a5*g21*h11*alpha*cth*beta*sth
	+2.*a5*g21*h11*g12*g22*beta*s2*alpha
	-2.*a5*g12*h12*alpha*c2*beta
	-4.*a5*g22*h12*alpha*cth*beta*sth
	+2.*a5*g11*h11*beta*s2*alpha
	-4.*a5*g21_3*h11*alpha*cth*beta*sth
	-2.*a5*g11*h11*alpha*c2*beta
	+2.*a5*g12*h12*g21_2*beta*s2*alpha
	+4.*a5*g11*h11*alpha*cth*beta*sth
	+4.*a5*g11*h11*g21_2*beta*s2*alpha
	-4.*a5*g11*h11*g21_2*alpha*c2*beta
	+2.*a5*g11*h11*g22_2*beta*s2*alpha
	-2.*a5*g11*h11*g22_2*alpha*c2*beta
	-2.*a5*g21*h11*g12*g22*alpha*c2*beta
	+4.*a5*g11*h11*g12*g22*alpha*cth*beta*sth
	+2.*a5*g22*h12*beta*s2*alpha
	-2.*a5*g22*h12*alpha*c2*beta
	-2.*a5*g12*h12*g21_2*alpha*c2*beta
	-4.*a5*g12*h12*g22_2*alpha*c2*beta
	-4.*a5*g22_3*h12*alpha*cth*beta*sth
	+2.*a5*g21*h11*beta*s2*alpha
	-2.*a5*g21*h11*alpha*c2*beta
	+4.*a5*g12*h12*g22_2*beta*s2*alpha)
	-2.*alpha*cth*(a5*g21_3*h11*beta*beta*c2
	-a5*g12*h12*beta*beta*sth*cth
	+a5*g11*h11*beta*beta*s2
	+a5*g21*h11*beta*beta*c2
	-a5*g22*h12*g11*g21*beta*beta*sth*cth
	-a5*g11*h11*beta*beta*sth*cth
	+a5*g12*h12*beta*beta*s2
	+a5*g22*h12*beta*beta*c2
	+a5*g22_3*h12*beta*beta*c2
	+a5*g11_2*h11*g21*beta*beta*s2
	-a5*g21*h11*g12*g22*beta*beta*sth*cth
	+a5*g22*h12*g21_2*beta*beta*c2
	+a5*g21*h11*g22_2*beta*beta*c2
	-a5*g21*h11*beta*beta*sth*cth
	+a5*g12*h12*g11*g21*beta*beta*s2
	-a5*g12*h12*g21_2*beta*beta*sth*cth
	-2.*a5*g11*h11*g21_2*beta*beta*sth*cth
	+a5*g11*h11*g12*g22*beta*beta*s2
	-a5*g11*h11*g22_2*beta*beta*sth*cth
	+a5*g12_2*h12*g22*beta*beta*s2
	-a5*g22*h12*beta*beta*sth*cth
	-2.*a5*g12*h12*g22_2*beta*beta*sth*cth
	-a4*a6*h11*g21
	-a4*a6*h12*g22);

	K2 = -beta*sth*(-2.*a5*g21_3*h11*beta*beta*c2
	+2.*a5*g12*h12*beta*beta*sth*cth
	-2.*a5*g11*h11*beta*beta*s2
	-2.*a5*g21*h11*beta*beta*c2
	+2.*a5*g22*h12*g11*g21*beta*beta*sth*cth
	+2.*a5*g11*h11*beta*beta*sth*cth
	-2.*a5*g12*h12*beta*beta*s2
	-2.*a5*g22*h12*beta*beta*c2
	-2.*a5*g22_3*h12*beta*beta*c2
	-2.*a5*g11_2*h11*g21*beta*beta*s2
	+2.*a5*g21*h11*g12*g22*beta*beta*sth*cth
	-2.*a5*g22*h12*g21_2*beta*beta*c2
	-2.*a5*g21*h11*g22_2*beta*beta*c2
	+2.*a5*g21*h11*beta*beta*sth*cth
	-2.*a5*g12*h12*g11*g21*beta*beta*s2
	+2.*a5*g12*h12*g21_2*beta*beta*sth*cth
	+4.*a5*g11*h11*g21_2*beta*beta*sth*cth
	-2.*a5*g11*h11*g12*g22*beta*beta*s2
	+2.*a5*g11*h11*g22_2*beta*beta*sth*cth
	-2.*a5*g12_2*h12*g22*beta*beta*s2
	+2.*a5*g22*h12*beta*beta*sth*cth
	+4.*a5*g12*h12*g22_2*beta*beta*sth*cth
	-2.*a4*a6*h11*g21
	-2.*a4*a6*h12*g22
	+4.*a5*g22_3*h12*alpha*alpha*s2
	+4.*a5*g21_3*h11*alpha*alpha*s2
	+4.*a5*g22*h12*alpha*alpha*s2
	+4.*a5*g12*h12*alpha*alpha*c2
	+4.*a5*g22*h12*g11*g21*alpha*alpha*cth*sth
	+4.*a5*g11*h11*g12*g22*alpha*alpha*c2
	+4.*a5*g12*h12*alpha*alpha*cth*sth
	+4.*a5*g21*h11*alpha*alpha*s2
	+4.*a5*g12_2*h12*g22*alpha*alpha*c2
	+4.*a5*g11*h11*alpha*alpha*cth*sth
	+4.*a5*g12*h12*g11*g21*alpha*alpha*c2
	+8.*a5*g11*h11*g21_2*alpha*alpha*cth*sth
	+4.*a5*g11*h11*alpha*alpha*c2
	+4.*a5*g11*h11*g22_2*alpha*alpha*cth*sth
	+4.*a5*g21*h11*g12*g22*alpha*alpha*cth*sth
	+4.*a5*g22*h12*alpha*alpha*cth*sth
	+4.*a5*g22*h12*g21_2*alpha*alpha*s2
	+4.*a5*g12*h12*g21_2*alpha*alpha*cth*sth
	+4.*a5*g11_2*h11*g21*alpha*alpha*c2
	+8.*a5*g12*h12*g22_2*alpha*alpha*cth*sth
	+4.*a5*g21*h11*alpha*alpha*cth*sth
	+4.*a5*g21*h11*g22_2*alpha*alpha*s2)
	-2.*alpha*cth*(2.*a5*g22*h12*g11*g21*beta*s2*alpha
	+2.*a5*g12*h12*beta*s2*alpha
	-2.*a5*g22*h12*g11*g21*alpha*c2*beta
	-4.*a5*g22*h12*g21_2*alpha*cth*beta*sth
	+4.*a5*g12_2*h12*g22*alpha*cth*beta*sth
	-4.*a5*g21*h11*g22_2*alpha*cth*beta*sth
	+4.*a5*g12*h12*g11*g21*alpha*cth*beta*sth
	+4.*a5*g11_2*h11*g21*alpha*cth*beta*sth
	+4.*a5*g12*h12*alpha*cth*beta*sth
	-4.*a5*g21*h11*alpha*cth*beta*sth
	+2.*a5*g21*h11*g12*g22*beta*s2*alpha
	-2.*a5*g12*h12*alpha*c2*beta
	-4.*a5*g22*h12*alpha*cth*beta*sth
	+2.*a5*g11*h11*beta*s2*alpha
	-4.*a5*g21_3*h11*alpha*cth*beta*sth
	-2.*a5*g11*h11*alpha*c2*beta
	+2.*a5*g12*h12*g21_2*beta*s2*alpha
	+4.*a5*g11*h11*alpha*cth*beta*sth
	+4.*a5*g11*h11*g21_2*beta*s2*alpha
	-4.*a5*g11*h11*g21_2*alpha*c2*beta
	+2.*a5*g11*h11*g22_2*beta*s2*alpha
	-2.*a5*g11*h11*g22_2*alpha*c2*beta
	-2.*a5*g21*h11*g12*g22*alpha*c2*beta
	+4.*a5*g11*h11*g12*g22*alpha*cth*beta*sth
	+2.*a5*g22*h12*beta*s2*alpha
	-2.*a5*g22*h12*alpha*c2*beta
	-2.*a5*g12*h12*g21_2*alpha*c2*beta
	-4.*a5*g12*h12*g22_2*alpha*c2*beta
	-4.*a5*g22_3*h12*alpha*cth*beta*sth
	+2.*a5*g21*h11*beta*s2*alpha
	-2.*a5*g21*h11*alpha*c2*beta
	+4.*a5*g12*h12*g22_2*beta*s2*alpha)
	+beta*sth*(a5*g21_3*h11*beta*beta*c2
	-a5*g12*h12*beta*beta*sth*cth
	+a5*g11*h11*beta*beta*s2
	+a5*g21*h11*beta*beta*c2
	-a5*g22*h12*g11*g21*beta*beta*sth*cth
	-a5*g11*h11*beta*beta*sth*cth
	+a5*g12*h12*beta*beta*s2
	+a5*g22*h12*beta*beta*c2
	+a5*g22_3*h12*beta*beta*c2
	+a5*g11_2*h11*g21*beta*beta*s2
	-a5*g21*h11*g12*g22*beta*beta*sth*cth
	+a5*g22*h12*g21_2*beta*beta*c2
	+a5*g21*h11*g22_2*beta*beta*c2
	-a5*g21*h11*beta*beta*sth*cth
	+a5*g12*h12*g11*g21*beta*beta*s2
	-a5*g12*h12*g21_2*beta*beta*sth*cth
	-2.*a5*g11*h11*g21_2*beta*beta*sth*cth
	+a5*g11*h11*g12*g22*beta*beta*s2
	-a5*g11*h11*g22_2*beta*beta*sth*cth
	+a5*g12_2*h12*g22*beta*beta*s2
	-a5*g22*h12*beta*beta*sth*cth
	-2.*a5*g12*h12*g22_2*beta*beta*sth*cth
	-a4*a6*h11*g21
	-a4*a6*h12*g22);

	K3 = -beta*sth*(-2.*a5*g22*h12*g11*g21*beta*s2*alpha
	-2.*a5*g12*h12*beta*s2*alpha
	+2.*a5*g22*h12*g11*g21*alpha*c2*beta
	+4.*a5*g22*h12*g21_2*alpha*cth*beta*sth
	-4.*a5*g12_2*h12*g22*alpha*cth*beta*sth
	+4.*a5*g21*h11*g22_2*alpha*cth*beta*sth
	-4.*a5*g12*h12*g11*g21*alpha*cth*beta*sth
	-4.*a5*g11_2*h11*g21*alpha*cth*beta*sth
	-4.*a5*g12*h12*alpha*cth*beta*sth
	+4.*a5*g21*h11*alpha*cth*beta*sth
	-2.*a5*g21*h11*g12*g22*beta*s2*alpha
	+2.*a5*g12*h12*alpha*c2*beta
	+4.*a5*g22*h12*alpha*cth*beta*sth
	-2.*a5*g11*h11*beta*s2*alpha
	+4.*a5*g21_3*h11*alpha*cth*beta*sth
	+2.*a5*g11*h11*alpha*c2*beta
	-2.*a5*g12*h12*g21_2*beta*s2*alpha
	-4.*a5*g11*h11*alpha*cth*beta*sth
	-4.*a5*g11*h11*g21_2*beta*s2*alpha
	+4.*a5*g11*h11*g21_2*alpha*c2*beta
	-2.*a5*g11*h11*g22_2*beta*s2*alpha
	+2.*a5*g11*h11*g22_2*alpha*c2*beta
	+2.*a5*g21*h11*g12*g22*alpha*c2*beta
	-4.*a5*g11*h11*g12*g22*alpha*cth*beta*sth
	-2.*a5*g22*h12*beta*s2*alpha
	+2.*a5*g22*h12*alpha*c2*beta
	+2.*a5*g12*h12*g21_2*alpha*c2*beta
	+4.*a5*g12*h12*g22_2*alpha*c2*beta
	+4.*a5*g22_3*h12*alpha*cth*beta*sth
	-2.*a5*g21*h11*beta*s2*alpha
	+2.*a5*g21*h11*alpha*c2*beta
	-4.*a5*g12*h12*g22_2*beta*s2*alpha)
	-2.*alpha*cth*(-2.*a5*g21_3*h11*beta*beta*c2
	+2.*a5*g12*h12*beta*beta*sth*cth
	-2.*a5*g11*h11*beta*beta*s2
	-2.*a5*g21*h11*beta*beta*c2
	+2.*a5*g22*h12*g11*g21*beta*beta*sth*cth
	+2.*a5*g11*h11*beta*beta*sth*cth
	-2.*a5*g12*h12*beta*beta*s2
	-2.*a5*g22*h12*beta*beta*c2
	-2.*a5*g22_3*h12*beta*beta*c2
	-2.*a5*g11_2*h11*g21*beta*beta*s2
	+2.*a5*g21*h11*g12*g22*beta*beta*sth*cth
	-2.*a5*g22*h12*g21_2*beta*beta*c2
	-2.*a5*g21*h11*g22_2*beta*beta*c2
	+2.*a5*g21*h11*beta*beta*sth*cth
	-2.*a5*g12*h12*g11*g21*beta*beta*s2
	+2.*a5*g12*h12*g21_2*beta*beta*sth*cth
	+4.*a5*g11*h11*g21_2*beta*beta*sth*cth
	-2.*a5*g11*h11*g12*g22*beta*beta*s2
	+2.*a5*g11*h11*g22_2*beta*beta*sth*cth
	-2.*a5*g12_2*h12*g22*beta*beta*s2
	+2.*a5*g22*h12*beta*beta*sth*cth
	+4.*a5*g12*h12*g22_2*beta*beta*sth*cth
	-2.*a4*a6*h11*g21
	-2.*a4*a6*h12*g22
	+4.*a5*g22_3*h12*alpha*alpha*s2
	+4.*a5*g21_3*h11*alpha*alpha*s2
	+4.*a5*g22*h12*alpha*alpha*s2
	+4.*a5*g12*h12*alpha*alpha*c2
	+4.*a5*g22*h12*g11*g21*alpha*alpha*cth*sth
	+4.*a5*g11*h11*g12*g22*alpha*alpha*c2
	+4.*a5*g12*h12*alpha*alpha*cth*sth
	+4.*a5*g21*h11*alpha*alpha*s2
	+4.*a5*g12_2*h12*g22*alpha*alpha*c2
	+4.*a5*g11*h11*alpha*alpha*cth*sth
	+4.*a5*g12*h12*g11*g21*alpha*alpha*c2
	+8.*a5*g11*h11*g21_2*alpha*alpha*cth*sth
	+4.*a5*g11*h11*alpha*alpha*c2
	+4.*a5*g11*h11*g22_2*alpha*alpha*cth*sth
	+4.*a5*g21*h11*g12*g22*alpha*alpha*cth*sth
	+4.*a5*g22*h12*alpha*alpha*cth*sth
	+4.*a5*g22*h12*g21_2*alpha*alpha*s2
	+4.*a5*g12*h12*g21_2*alpha*alpha*cth*sth
	+4.*a5*g11_2*h11*g21*alpha*alpha*c2
	+8.*a5*g12*h12*g22_2*alpha*alpha*cth*sth
	+4.*a5*g21*h11*alpha*alpha*cth*sth
	+4.*a5*g21*h11*g22_2*alpha*alpha*s2)
	+beta*sth*(2.*a5*g22*h12*g11*g21*beta*s2*alpha
	+2.*a5*g12*h12*beta*s2*alpha
	-2.*a5*g22*h12*g11*g21*alpha*c2*beta
	-4.*a5*g22*h12*g21_2*alpha*cth*beta*sth
	+4.*a5*g12_2*h12*g22*alpha*cth*beta*sth
	-4.*a5*g21*h11*g22_2*alpha*cth*beta*sth
	+4.*a5*g12*h12*g11*g21*alpha*cth*beta*sth
	+4.*a5*g11_2*h11*g21*alpha*cth*beta*sth
	+4.*a5*g12*h12*alpha*cth*beta*sth
	-4.*a5*g21*h11*alpha*cth*beta*sth
	+2.*a5*g21*h11*g12*g22*beta*s2*alpha
	-2.*a5*g12*h12*alpha*c2*beta
	-4.*a5*g22*h12*alpha*cth*beta*sth
	+2.*a5*g11*h11*beta*s2*alpha
	-4.*a5*g21_3*h11*alpha*cth*beta*sth
	-2.*a5*g11*h11*alpha*c2*beta
	+2.*a5*g12*h12*g21_2*beta*s2*alpha
	+4.*a5*g11*h11*alpha*cth*beta*sth
	+4.*a5*g11*h11*g21_2*beta*s2*alpha
	-4.*a5*g11*h11*g21_2*alpha*c2*beta
	+2.*a5*g11*h11*g22_2*beta*s2*alpha
	-2.*a5*g11*h11*g22_2*alpha*c2*beta
	-2.*a5*g21*h11*g12*g22*alpha*c2*beta
	+4.*a5*g11*h11*g12*g22*alpha*cth*beta*sth
	+2.*a5*g22*h12*beta*s2*alpha
	-2.*a5*g22*h12*alpha*c2*beta
	-2.*a5*g12*h12*g21_2*alpha*c2*beta
	-4.*a5*g12*h12*g22_2*alpha*c2*beta
	-4.*a5*g22_3*h12*alpha*cth*beta*sth
	+2.*a5*g21*h11*beta*s2*alpha
	-2.*a5*g21*h11*alpha*c2*beta
	+4.*a5*g12*h12*g22_2*beta*s2*alpha);

  L0 = a4*a6*g11*g21
	+a4*a6*g12*g22
	+a4*a6
	-a5*g22_2*beta*beta*c2
	-a5*g21_2*beta*beta*c2
	+a5*g11_2*g22_2*beta*beta*sth*cth
	-a5*beta*beta*s2
	-a5*beta*beta*c2
	-a5*g12_3*g22*beta*beta*s2
	-a5*g12_2*g11*g21*beta*beta*s2
	-a5*g21_3*g11*beta*beta*c2
	-a5*g22*g12*beta*beta*c2
	-a5*g21*g11*beta*beta*s2
	-a5*g22_3*g12*beta*beta*c2
	+2.*a5*g21*g11*beta*beta*sth*cth
	-a5*g22*g12*beta*beta*s2
	+2.*a5*g21*g11*g12*g22*beta*beta*sth*cth
	-a5*g21*g11*beta*beta*c2
	+a5*g12_2*beta*beta*sth*cth
	-a5*g12_2*beta*beta*s2
	+2.*a5*g11_2*g21_2*beta*beta*sth*cth
	+2.*a5*g22*g12*beta*beta*sth*cth
	+2.*a5*beta*beta*sth*cth
	-a5*g22*g12*g21_2*beta*beta*c2
	+2.*a5*g12_2*g22_2*beta*beta*sth*cth
	+a5*g12_2*g21_2*beta*beta*sth*cth
	+a5*g11_2*beta*beta*sth*cth
	-a5*g11_2*g12*g22*beta*beta*s2
	-a5*g11_3*g21*beta*beta*s2
	+a5*g21_2*beta*beta*sth*cth
	+a5*g22_2*beta*beta*sth*cth
	-a5*g21*g11*g22_2*beta*beta*c2
	-a5*g11_2*beta*beta*s2;

  L1 = 4.*a5*g22_3*g12*alpha*cth*beta*sth
	+4.*a5*g21_3*g11*alpha*cth*beta*sth
	+4.*a5*g21_2*alpha*cth*beta*sth
	+4.*a5*g21*g11*g12*g22*alpha*c2*beta
	-4.*a5*g21*g11*g12*g22*beta*s2*alpha
	+4.*a5*g21*g11*alpha*c2*beta
	-4.*a5*g21*g11*beta*s2*alpha
	-4.*a5*g12_2*g11*g21*alpha*cth*beta*sth
	-4.*a5*g12_3*g22*alpha*cth*beta*sth
	+2.*a5*g11_2*alpha*c2*beta
	-2.*a5*g11_2*beta*s2*alpha
	+4.*a5*g12_2*g22_2*alpha*c2*beta
	+4.*a5*alpha*c2*beta
	-4.*a5*beta*s2*alpha
	-4.*a5*g12_2*alpha*cth*beta*sth
	+4.*a5*g22*g12*g21_2*alpha*cth*beta*sth
	+4.*a5*g22*g12*alpha*c2*beta
	-4.*a5*g22*g12*beta*s2*alpha
	-2.*a5*g21_2*beta*s2*alpha
	-4.*a5*g11_2*g21_2*beta*s2*alpha
	+4.*a5*g11_2*g21_2*alpha*c2*beta
	+4.*a5*g22_2*alpha*cth*beta*sth
	+2.*a5*g12_2*alpha*c2*beta
	-2.*a5*g12_2*beta*s2*alpha
	+2.*a5*g11_2*g22_2*alpha*c2*beta
	-2.*a5*g11_2*g22_2*beta*s2*alpha
	+4.*a5*g21*g11*g22_2*alpha*cth*beta*sth
	-4.*a5*g11_3*g21*alpha*cth*beta*sth
	+2.*a5*g21_2*alpha*c2*beta
	+2.*a5*g22_2*alpha*c2*beta
	-2.*a5*g22_2*beta*s2*alpha
	-4.*a5*g11_2*g12*g22*alpha*cth*beta*sth
	-4.*a5*g11_2*alpha*cth*beta*sth
	-4.*a5*g12_2*g22_2*beta*s2*alpha
	+2.*a5*g12_2*g21_2*alpha*c2*beta
	-2.*a5*g12_2*g21_2*beta*s2*alpha;

  L2 = (2*a4*a6*g11*g21)
	+(2*a4*a6*g12*g22)
	+(2*a4*a6)
	+2.*a5*(g22_2)*beta*beta*c2
	+2.*a5*(g21_2)*beta*beta*c2
	-2.*a5*(g11_2)*(g22_2)*beta*beta*sth*cth
	+2.*a5*beta*beta*s2
	+2.*a5*beta*beta*c2
	+2.*a5*g12_3*g22*beta*beta*s2
	+2.*a5*(g12_2)*g11*g21*beta*beta*s2
	+2.*a5*g21_3*g11*beta*beta*c2
	+2.*a5*g22*g12*beta*beta*c2
	+2.*a5*g21*g11*beta*beta*s2
	+2.*a5*g22_3*g12*beta*beta*c2
	-4.*a5*g21*g11*beta*beta*sth*cth
	+2.*a5*g22*g12*beta*beta*s2
	-4.*a5*g21*g11*g12*g22*beta*beta*sth*cth
	+2.*a5*g21*g11*beta*beta*c2
	-2.*a5*(g12_2)*beta*beta*sth*cth
	+2.*a5*(g12_2)*beta*beta*s2
	-4.*a5*(g11_2)*(g21_2)*beta*beta*sth*cth
	-4.*a5*g22*g12*beta*beta*sth*cth
	-4.*a5*beta*beta*sth*cth
	+2.*a5*g22*g12*(g21_2)*beta*beta*c2
	-4.*a5*(g12_2)*(g22_2)*beta*beta*sth*cth
	-2.*a5*(g12_2)*(g21_2)*beta*beta*sth*cth
	-2.*a5*(g11_2)*beta*beta*sth*cth
	+2.*a5*(g11_2)*g12*g22*beta*beta*s2
	+2.*a5*g11_3*g21*beta*beta*s2
	-2.*a5*(g21_2)*beta*beta*sth*cth
	-2.*a5*(g22_2)*beta*beta*sth*cth
	+2.*a5*g21*g11*(g22_2)*beta*beta*c2
	+2.*a5*(g11_2)*beta*beta*s2
	-4.*a5*alpha*alpha*c2
	-4.*a5*alpha*alpha*s2
	-4.*a5*g22_3*g12*alpha*alpha*s2
	-4.*a5*(g21_2)*alpha*alpha*s2
	-8.*a5*(g11_2)*(g21_2)*alpha*alpha*cth*sth
	-4.*a5*g22*g12*alpha*alpha*c2
	-4.*a5*(g11_2)*alpha*alpha*c2
	-4.*a5*(g22_2)*alpha*alpha*s2
	-8.*a5*alpha*alpha*cth*sth
	-4.*a5*g21_3*g11*alpha*alpha*s2
	-4.*a5*g22*g12*alpha*alpha*s2
	-4.*a5*(g12_2)*g11*g21*alpha*alpha*c2
	-4.*a5*g21*g11*alpha*alpha*c2
	-4.*a5*g12_3*g22*alpha*alpha*c2
	-4.*a5*(g11_2)*alpha*alpha*cth*sth
	-8.*a5*(g12_2)*(g22_2)*alpha*alpha*cth*sth
	-4.*a5*(g12_2)*alpha*alpha*c2
	-8.*a5*g21*g11*g12*g22*alpha*alpha*cth*sth
	-8.*a5*g22*g12*alpha*alpha*cth*sth
	-4.*a5*(g21_2)*alpha*alpha*cth*sth
	-4.*a5*g21*g11*alpha*alpha*s2
	-4.*a5*(g11_2)*(g22_2)*alpha*alpha*cth*sth
	-4.*a5*(g12_2)*(g21_2)*alpha*alpha*cth*sth
	-4.*a5*g11_3*g21*alpha*alpha*c2
	-4.*a5*(g22_2)*alpha*alpha*cth*sth
	-4.*a5*(g11_2)*g12*g22*alpha*alpha*c2
	-8.*a5*g21*g11*alpha*alpha*cth*sth
	-4.*a5*g22*g12*(g21_2)*alpha*alpha*s2
	-4.*a5*(g12_2)*alpha*alpha*cth*sth
	-4.*a5*g21*g11*(g22_2)*alpha*alpha*s2;
}
  
// sigmas are basically depth derivatives which are intermediate variables
// needed to build R
template<typename T>
void pose_poly<T>::
get_sigmas(const unsigned char ts_len, const T (&ts)[TS_MAX_LEN],
	T (*sigmas)[2][TS_MAX_LEN][4], unsigned char sigmas_len[TS_MAX_LEN])
{
	T   (&sigmas1)[TS_MAX_LEN][4] = (*sigmas)[0];
	T   (&sigmas2)[TS_MAX_LEN][4] = (*sigmas)[1];
	T p[10];
	for (unsigned char i = 0; i < ts_len; i++) {
		sigmas_len[i] = 0; 

		fn_t(ts[i], p); // TODO: perhaps normalize H J K L

		T &A = p[0], &B = p[1], &C = p[2], 
      &E = p[3], &F = p[4], &G = p[5], &H = p[6],
      &J = p[7], &K = p[8], &L = p[9];
    A += A; E += E;

    T delta = B*B - 2*A*C;
    if (delta < -1e-4) continue;
    delta = (delta < 0)? 0 : sqrt(delta);
    B = -B;
		T sigma1_m = (B - delta)/A;
		T sigma1_p = (B + delta)/A;

		delta = F*F - 2*E*G;
    if (delta < -1e-4) continue;
    delta = (delta < 0)? 0 : sqrt(delta);
    F = -F;
		T sigma2_m = (F - delta)/E;
		T sigma2_p = (F + delta)/E;

		static constexpr T my_eps = 1.0; // TODO optimize, and use my_eps a fraction
                                     // of norm H J K L norm
		if (std::fabs(H + J*sigma1_m + K*sigma2_m + L*sigma1_m*sigma2_m) < my_eps) {
			sigmas1[i][sigmas_len[i]]   = sigma1_m;
			sigmas2[i][sigmas_len[i]++] = sigma2_m;
		}
		if (std::fabs(H + J*sigma1_p + K*sigma2_m + L*sigma1_p*sigma2_m) < my_eps) {
			sigmas1[i][sigmas_len[i]]   = sigma1_p;
			sigmas2[i][sigmas_len[i]++] = sigma2_m;
		}
		if (std::fabs(H + J*sigma1_p + K*sigma2_p + L*sigma1_p*sigma2_p) < my_eps) {
			sigmas1[i][sigmas_len[i]]   = sigma1_p;
			sigmas2[i][sigmas_len[i]++] = sigma2_p;
		}
		if (std::fabs(H + J*sigma1_m + K*sigma2_p + L*sigma1_m*sigma2_p) < my_eps) {
			sigmas1[i][sigmas_len[i]]   = sigma1_m;
			sigmas2[i][sigmas_len[i]++] = sigma2_p;
		}
	}
}

template<typename T>
inline void
invm3x3(T (&M)[3][3])
{
	// 3x3 MATRIX INVERSION ALGORITHM
	//             -1               T
	//  -1  [a b c]      1   [A B C]      1   [A D G]
	// M  = [d e f] = ------ [D E F] = ------ [B E H]
	//      [g h i]   det(M) [G H I]   det(M) [C F I]
	//
	// A =  (ei - fh), D = -(bi - ch), G =  (bf - ce),
	// B = -(di - fg), E =  (ai - cg), H = -(af - cd),
	// C =  (dh - eg), F = -(ah - bg), I =  (ae - bd).

	const T 
	a = M[0][0], b = M[0][1], c = M[0][2],
	d = M[1][0], e = M[1][1], f = M[1][2],
	g = M[2][0], h = M[2][1], i = M[2][2];

	const T 
	A =  (e*i - f*h), B = -(d*i - f*g), C =  (d*h - e*g),
	D = -(b*i - c*h), E =  (a*i - c*g), F = -(a*h - b*g),
	G =  (b*f - c*e), H = -(a*f - c*d), I =  (a*e - b*d);

	const T invdet_M = 1. / (a*A + b*B + c*C);
	M[0][0] = invdet_M * A; M[0][1] = invdet_M * D; M[0][2] = invdet_M * G;
	M[1][0] = invdet_M * B; M[1][1] = invdet_M * E; M[1][2] = invdet_M * H;
	M[2][0] = invdet_M * C; M[2][1] = invdet_M * F; M[2][2] = invdet_M * I;
}

template<typename T>
inline void
multm3x3(const T (&m1)[3][3], const T (&m2)[3][3], T output_m[][3])
{
	for (unsigned char i = 0; i < 3; i++)
		for (unsigned char j = 0; j < 3; j++) {
			output_m[i][j] = 0;
			for (unsigned char k = 0; k < 3; k++)
				output_m[i][j] += m1[i][k] * m2[k][j];
		}
}

template<typename T>
void pose_poly<T>::
get_r_t_from_rhos(
	const unsigned char ts_len,
	const T sigmas1[TS_MAX_LEN][4], const unsigned char sigmas_len[TS_MAX_LEN],
	const T sigmas2[TS_MAX_LEN][4],
	const T rhos1[TS_MAX_LEN], const T rhos2[TS_MAX_LEN],
	const T gama1[3], const T tgt1[3],
	const T gama2[3], const T tgt2[3],
	const T Gama1[3], const T Tgt1[3],
	const T Gama2[3], const T Tgt2[3],
	T (*output)[RT_MAX_LEN][4][3], unsigned char *output_len
)
{
	T lambdas1[TS_MAX_LEN][TS_MAX_LEN]; T lambdas2[TS_MAX_LEN][TS_MAX_LEN];
	const T DGama[3] = {Gama1[0]-Gama2[0], Gama1[1]-Gama2[1], Gama1[2]-Gama2[2]};
  
	for (unsigned char i = 0; i < ts_len; i++) {
    const T dgamas_rhos[3] = {
     rhos1[i]*gama1[0] - rhos2[i]*gama2[0],
     rhos1[i]*gama1[1] - rhos2[i]*gama2[1],
     rhos1[i]*gama1[2] - rhos2[i]*gama2[2]};
		for (unsigned char j = 0; j < sigmas_len[i]; j++) {
			lambdas1[i][j] = 
        (DGama[0]*Tgt1[0]+DGama[1]*Tgt1[1] + DGama[2]*Tgt1[2]) / 
        (dgamas_rhos[0]*(rhos1[i]*tgt1[0]  + sigmas1[i][j]*gama1[0]) + 
        dgamas_rhos[1]*(rhos1[i]*tgt1[1]   + sigmas1[i][j]*gama1[1]) +
        dgamas_rhos[2]*(rhos1[i]*tgt1[2]   + sigmas1[i][j]*gama1[2]));
			lambdas2[i][j] = 
        (DGama[0]*Tgt2[0]+DGama[1]*Tgt2[1] + DGama[2]*Tgt2[2]) /
        (dgamas_rhos[0]*(rhos2[i]*tgt2[0]  + sigmas2[i][j]*gama2[0]) + 
        dgamas_rhos[1]*(rhos2[i]*tgt2[1]   + sigmas2[i][j]*gama2[1]) +
        dgamas_rhos[2]*(rhos2[i]*tgt2[2]   + sigmas2[i][j]*gama2[2]));
		}
	}

	//% Rotation:
	T A[3][3] = {
		DGama[0], Tgt1[0], Tgt2[0],
		DGama[1], Tgt1[1], Tgt2[1],
		DGama[2], Tgt1[2], Tgt2[2]};
	invm3x3(A);

	// Matrix containing Rotations and Translations
	T (&RT)[RT_MAX_LEN][4][3] = *output;
	unsigned char &RT_len     = *output_len; RT_len = 0;
	for (unsigned char i = 0; i < ts_len; i++) {
		for (unsigned char j = 0; j < sigmas_len[i]; j++, RT_len++) {
      assert(RT_len < RT_MAX_LEN);
			T (&Rots)[4][3] = RT[RT_len]; T (&Transls)[3] = RT[RT_len][3];

			#define B_row(r) \
				rhos1[i]*gama1[(r)] - rhos2[i]*gama2[(r)], \
				lambdas1[i][j]*(rhos1[i]*tgt1[(r)] + sigmas1[i][j]*gama1[(r)]), \
				lambdas2[i][j]*(rhos2[i]*tgt2[(r)] + sigmas2[i][j]*gama2[(r)])

			const T B[3][3] = { B_row(0), B_row(1), B_row(2) };
			multm3x3(B, A, Rots);
      Transls[0] = rhos1[i]*gama1[0] - Rots[0][0]*Gama1[0] - Rots[0][1]*Gama1[1] - Rots[0][2]*Gama1[2];
      Transls[1] = rhos1[i]*gama1[1] - Rots[1][0]*Gama1[0] - Rots[1][1]*Gama1[1] - Rots[1][2]*Gama1[2];
      Transls[2] = rhos1[i]*gama1[2] - Rots[2][0]*Gama1[0] - Rots[2][1]*Gama1[1] - Rots[2][2]*Gama1[2];
		}
	}
}

void P2PtSolver_Fabbri::Solve(
    const Mat &bearing_vectors,
    const Mat &tangent_vectors,
    const Mat &X, // 3D points
    const Mat &T, // 3D tangents
    std::vector<Mat34> *models)
{
//  OPENMVG_LOG_INFO << bearing_vectors.rows() << std::endl;
  assert(3 == bearing_vectors.rows());
  assert(3 == X.rows());
  assert(bearing_vectors.cols() == X.cols());
  assert(bearing_vectors.cols() == 2);

  unsigned char nsols;
  double degen;
	double rotation_translation_solutions[RT_MAX_LEN][4][3];

  // if (!
    p2pt<double>::pose_from_point_tangents(
    bearing_vectors.col(0).data(), tangent_vectors.col(0).data(),
    bearing_vectors.col(1).data(), tangent_vectors.col(1).data(),
    X.col(0).data(), T.col(0).data(),
    X.col(1).data(), T.col(1).data(),
    &rotation_translation_solutions, &nsols, &degen);
  //  )
    // OPENMVG_LOG_ERROR << "degeneracy"; 

  OPENMVG_LOG_INFO << "Number of models returned p2pt: " << (int)nsols;
	for (unsigned char i = 0; i < nsols; ++i) {
    Mat34 P;
    for (unsigned char j = 0 ; j < 3; ++j)
      for (unsigned char k = 0 ; k < 3; ++k)
        P(j,k) = rotation_translation_solutions[i][j][k];

    for (unsigned char k = 0 ; k < 3; ++k)
      P(k,3) = rotation_translation_solutions[i][3][k];
    models->push_back(P);
  }
};

} // namespace euclidean_resection
} // namespace openMVG
