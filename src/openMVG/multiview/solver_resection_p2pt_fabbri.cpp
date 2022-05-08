// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2020 Ricardo Fabbri, Ariel Kovaljski
//
// Author: Ariel Kovaljski and Ricardo Fabbri
// Rio de Janeiro State University

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/multiview/solver_resection_p2pt_fabbri.hpp"
#include "openMVG/multiview/projection.hpp"

#include <array>
#include <complex>

namespace openMVG
{
namespace euclidean_resection
{
  
// At most 8 solutions with positive depth, TODO: assert if longer
static constexpr unsigned TS_MAX_LEN = 8;
static constexpr unsigned RT_MAX_LEN = (TS_MAX_LEN * TS_MAX_LEN);

template <typename T=double>
class p2pt { // fully static, not to be instantiated - just used for templating
	public:
  static void pose_from_point_tangents(
    const T gama1[3], const T tgt1[3],
    const T gama2[3], const T tgt2[3],
    const T Gama1[3], const T Tgt1[3],
    const T Gama2[3], const T Tgt2[3],
    T (*output_RT)[RT_MAX_LEN][4][3],
    unsigned *output_RT_len,
    T *output_degen
  );
};

template<typename T>
struct pose_poly {
	T A0, A1, A2, B0, B1, B2, B3, C0, C1, C2, C3, C4,
		E0, E1, E2, F0, F1, F2, F3, G0, G1, G2, G3, G4,
		H0, H1, H2, H3, H4, J0, J1, J2, J3, K0, K1, K2, K3,
		L0, L1, L2, alpha, beta, theta, sth, cth;

  static constexpr unsigned T_LEN = 1900, ROOT_IDS_LEN = T_LEN - 1;
  static constexpr double T_LEN_2 = 2./T_LEN;
  inline T t_vec(unsigned i) { return T_LEN_2*i -1.; }

	void pose_from_point_tangents_2(
		const T gama1[3], const T tgt1[3], const T gama2[3], const T tgt2[3],
		const T Gama1[3], const T Tgt1[3], const T Gama2[3], const T Tgt2[3]);
  
	inline void find_bounded_root_intervals(T (*root_ids_out)[ROOT_IDS_LEN]) {
    T curr_val = fn_t(t_vec(0)), next_val;
    for (unsigned i = 0; i < ROOT_IDS_LEN; i++) {
      next_val = fn_t(t_vec(i+1));
      (*root_ids_out)[i] = (curr_val * next_val) < 0;
      curr_val = next_val;
    }
  }
  
	inline T fn_t(const T t, T o[10]) {  //o: output
    T &A = o[0]; T &B = o[1]; T &C = o[2]; T &E = o[3]; T &F = o[4]; 
    T &G = o[5]; T &H = o[6]; T &J = o[7]; T &K = o[8]; T &L = o[9];

    //%function of t part :
    const T t2 = t*t, t3 = t2*t, t4 = t3*t, t5 = t4*t, t6 = t5*t, t7 = t6*t, t8 = t7*t;
    const T t2p12 = (t2 + 1.) * (t2 + 1.), t2p13 = t2p12 * (t2 + 1.), t2p14 = t2p13 * (t2 + 1.);

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

    const T A2=A*A, B2=B*B, C2=C*C, E2=E*E, F2=F*F, G2=G*G, H2=H*H, 
            H3=H2*H, H4=H3*H, J2=J*J, J3=J2*J, K2=K*K, K3=K2*K, L2=L*L, L3=L2*L;
    
    return E2*B2*H2*J2 +G2*C2*L*L*L*L +G2*A2*K3*K +E2*A2*H4 +E2*C2*J3*J
    +-2.*E*A*H2*G*C*L2 +2.*E2*A*H2*C*J2 +-2.*E2*C*J3*B*H +2.*E*C2*J2*G*L2
    +2.*E*A2*H2*G*K2 +-2.*E2*A*H3*B*J +-2.*E*A*H2*G*B*K*L +-2.*E*C*J2*G*B*K*L
    +-2.*E*C*J2*G*A*K2 +-2.*E*B*H*J*G*C*L2 +-2.*E*B*H*J*G*A*K2 +G2*B2*K2*L2
    +-2.*G2*B*K*L3*C +-2.*G2*B*K3*L*A +2.*G2*C*L2*A*K2 +-2.*F*E*A2*H3*K
    +-2.*F*E*A*H*K*C*J2 +3.*F*E*A*H2*K*B*J +3.*F*A*H*K2*G*B*L
    +-2.*F*A*H*K*G*C*L2 +-2.*F*A2*H*K3*G +F*E*B*H3*L*A +3.*F*E*B*H*L*C*J2
    +-1.*F*E*B2*H2*L*J +-1.*F*B2*H*L2*G*K +F*B*H*L3*G*C +F*E*B*K*J3*C
    +-1.*F*E*B2*K*J2*H +-1.*F*B2*K2*J*G*L +3.*F*B*K*J*G*C*L2 +F*B*K3*J*G*A
    +-2.*F*E*C*J*L*A*H2 +-2.*F*E*C2*J3*L +-2.*F*C2*J*L3*G +-2.*F*C*J*L*G*A*K2
    +F2*A2*K2*H2 +F2*A*K2*C*J2 +-1.*F2*A*K2*B*H*J +-1.*F2*B*K*L*A*H2
    +-1.*F2*B*K*L*C*J2 +F2*B2*K*L*H*J +F2*C*L2*A*H2 +F2*C2*L2*J2
    +-1.*F2*C*L2*B*H*J +G*E*B2*H2*L2 +G*E*B2*K2*J2 +8.*G*E*A*H*K*C*J*L;
  }
  
	inline T fn_t(const T t) { T b[10]; return fn_t(t, b); }
	inline T operator()(T t) { return fn_t(t); }
  
	inline void rhos_from_root_ids(
      const T (&root_ids)[ROOT_IDS_LEN], T (*out)[3][ROOT_IDS_LEN], unsigned *out_ts_len) { 
    T (&ts)[ROOT_IDS_LEN] = (*out)[0];
    unsigned &ts_end = *out_ts_len; ts_end = 0;
    for (unsigned i = 0; i < ROOT_IDS_LEN; i++) {
      if (!root_ids[i]) continue;
      T t0 = t_vec(i), t1 = t_vec(i+1), &t2 = ts[ts_end++];
      T f0 = fn_t(t_vec(i)), f1 = fn_t(t_vec(i+1));
      for (unsigned k = 0; k < 3; ++k) {
        t2 = t1 - f1*(t1-t0)/(f1-f0); t0 = t1; t1 = t2;
        f0 = f1; if (k + 1 < 3) f1 = fn_t(t2);
      }
    }
    //% Each root is now ts(i), plus minus t_stddev. Now get rho1(t):
    T (&rhos1)[ROOT_IDS_LEN] = (*out)[1]; T (&rhos2)[ROOT_IDS_LEN] = (*out)[2];
    const T alpha_times_2 = 2.*alpha;
    for (unsigned i = 0; i < ts_end; i++) {
      const T ts_new = ts[i],
      x2 = ts_new * ts_new,
      ts_den = 1. + x2,
      alpha_ts_new2 = alpha_times_2 * ts_new,
      beta_1_minus_x2 = beta * (1. - x2);
      rhos1[i] = ( alpha_ts_new2 * cth + beta_1_minus_x2 * sth) / ts_den;
      rhos2[i] = (-alpha_ts_new2 * sth + beta_1_minus_x2 * cth) / ts_den;
    }
  }
  
	void get_sigmas(const unsigned ts_len, const T (&ts)[ROOT_IDS_LEN], 
      T (*out)[2][TS_MAX_LEN][TS_MAX_LEN], unsigned (*out_len)[2][TS_MAX_LEN]);
  
	void get_r_t_from_rhos(
		const unsigned ts_len,
		const T sigmas1[TS_MAX_LEN][TS_MAX_LEN], const unsigned sigmas1_len[TS_MAX_LEN],
		const T sigmas2[TS_MAX_LEN][TS_MAX_LEN],
		const T rhos1[ROOT_IDS_LEN], const T rhos2[ROOT_IDS_LEN],
		const T gama1[3], const T tgt1[3], const T gama2[3], const T tgt2[3],
		const T Gama1[3], const T Tgt1[3], const T Gama2[3], const T Tgt2[3],
		T (*out)[RT_MAX_LEN][4][3], unsigned *out_len
	);
};
  
// This is the main routine ----------------------------------------------------
template <typename T>
void p2pt<T>::
pose_from_point_tangents(
	const T gama1[3], const T tgt1[3], const T gama2[3], const T tgt2[3],
	const T Gama1[3], const T Tgt1[3], const T Gama2[3], const T Tgt2[3],
	T (*output_RT)[RT_MAX_LEN][4][3], unsigned *output_RT_len, T *output_degen
)
{
	T DGama[3] = { Gama1[0] - Gama2[0], Gama1[1] - Gama2[1], Gama1[2] - Gama2[2] };
  { // % test for geometric degeneracy -------------------------------
  const T norm = sqrt(DGama[0]*DGama[0] + DGama[1]*DGama[1] + DGama[2]*DGama[2]);
	// Matrix for degeneracy calculation
	const T d[3][3] = {
		DGama[0]/norm, Tgt1[0], Tgt2[0],
		DGama[1]/norm, Tgt1[1], Tgt2[1],
		DGama[2]/norm, Tgt1[2], Tgt2[2]
	};
	T &degen = *output_degen;
	degen = (d[0][0]*d[1][1]*d[2][2]+d[0][1]*d[1][2]*d[2][0]+d[0][2]*d[1][0]*d[2][1]) // det(d)
		     -(d[2][0]*d[1][1]*d[0][2]+d[2][1]*d[1][2]*d[0][0]+d[2][2]*d[1][0]*d[0][1]);

	if (std::abs(degen) < 1.0e-3) {
		*output_RT_len = 0;
		return;
	}
  }

	// % compute roots -------------------------------
	pose_poly<T> p;
	p.pose_from_point_tangents_2( gama1, tgt1, gama2, tgt2, Gama1, Tgt1, Gama2, Tgt2);

	T root_ids[pose_poly<T>::ROOT_IDS_LEN] __attribute__((aligned (16)));
	p.find_bounded_root_intervals(&root_ids);

	// % compute rhos, r, t --------------------------
	T rhos[3][pose_poly<T>::ROOT_IDS_LEN];
	unsigned ts_len;
	p.rhos_from_root_ids(root_ids, &rhos, &ts_len);

	const T (&ts)[pose_poly<T>::ROOT_IDS_LEN]    = rhos[0]; 
  const T (&rhos1)[pose_poly<T>::ROOT_IDS_LEN] = rhos[1]; 
  const T (&rhos2)[pose_poly<T>::ROOT_IDS_LEN] = rhos[2];
	T sigmas[2][TS_MAX_LEN][TS_MAX_LEN]; unsigned sigmas_len[2][TS_MAX_LEN];

	p.get_sigmas(ts_len, ts, &sigmas, &sigmas_len);

	T (&sigmas1)[TS_MAX_LEN][TS_MAX_LEN] = sigmas[0];
	T (&sigmas2)[TS_MAX_LEN][TS_MAX_LEN] = sigmas[1];
	const unsigned (&sigmas1_len)[TS_MAX_LEN] = sigmas_len[0];

	T (&RT)[RT_MAX_LEN][4][3] = *output_RT;
	unsigned &RT_len               = *output_RT_len;

	p.get_r_t_from_rhos( ts_len, sigmas1, sigmas1_len, sigmas2,
		rhos1, rhos2, gama1, tgt1, gama2, tgt2, Gama1, Tgt1, Gama2, Tgt2, 
    &RT, &RT_len);
}

template<typename T>
void pose_poly<T>::
pose_from_point_tangents_2(
	const T gama1[3], const T tgt1[3],
	const T gama2[3], const T tgt2[3],
	const T Gama1[3], const T Tgt1[3],
	const T Gama2[3], const T Tgt2[3]
)
{
	static constexpr T PI_OVER_2 = 3.141592653589793*0.5;

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
	if (theta < 0) theta += PI_OVER_2;
	sth = sin(theta); const T s2 = sth*sth; cth = cos(theta);
  const T c2 = cth*cth, c2th = 2.*c2-1., s2th = 2.*sth*cth;
  
	const T den1 = 2.*s2th*(g11*g21+g12*g22+1.)+c2th*(g11_2+g12_2-g21_2-g22_2),
	        den2 = g11_2 + g12_2 + g21_2 + g22_2 + 2.;
	beta  = sqrt(-2.*a1 / (den1 - den2)); // sqrt(t25)
	alpha = sqrt(2.*a1 / (den1 + den2));  // sqrt(t24)

	//% Coefficient code adapted from Maple ::: can be further cleaned up but works
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
get_sigmas(const unsigned ts_len, const T (&ts)[ROOT_IDS_LEN],
	T (*out)[2][TS_MAX_LEN][TS_MAX_LEN], unsigned (*out_len)[2][TS_MAX_LEN])
{
	/* `out`
	       [0] -> sigmas1[TS_MAX_LEN][TS_MAX_LEN]
	       [1] -> sigmas2[TS_MAX_LEN][TS_MAX_LEN]

	       `sigmasX` (can contain single values or array of values)
	             [0][0] -> float/double
	             [1][0] -> float/double
	             [2][0] -> float/double
	             [3][ ] -> [0] = flt/dbl, [1] = flt/dbl, [2] = flt/dbl, ...
	               .
	               .
	               .
	   `out_len`
	       [0] -> sigmas1_len[TS_MAX_LEN]
	       [1] -> sigmas2_len[TS_MAX_LEN]

	       `sigmasX_len` (single values)
	             [0] = int, [1] = int, [2] = int, ...  */

	T   (&sigmas1)[TS_MAX_LEN][TS_MAX_LEN] = (*out)[0];
	T   (&sigmas2)[TS_MAX_LEN][TS_MAX_LEN] = (*out)[1];
	unsigned (&sigmas1_len)[TS_MAX_LEN]         = (*out_len)[0];
	T pose_out[10];
	for (unsigned i = 0; i < ts_len; i++) {
		sigmas1_len[i] = 0; 

		fn_t(ts[i], pose_out);

		//T &fvalue = pose_out[0]; // double-checked: not used in matlab
		const T &A = pose_out[0], &B = pose_out[1], &C = pose_out[2], 
            &E = pose_out[3], &F = pose_out[4], &G = pose_out[5], &H = pose_out[6],
            &J = pose_out[7], &K = pose_out[8], &L = pose_out[9];

		std::complex<T> delta = sqrt(B*B - 4*A*C);
		std::complex<T> sigma1_m = (-B - delta)/(2*A);
		std::complex<T> sigma1_p = (-B + delta)/(2*A);

		delta = sqrt(F*F - 4*E*G);
		std::complex<T> sigma2_m = (-F - delta)/(2*E);
		std::complex<T> sigma2_p = (-F + delta)/(2*E);

		//% handle case of negative delta
		if (std::abs(std::imag(sigma1_m)) < 1e-4) {
			sigma1_m = std::real(sigma1_m);
			sigma1_p = std::real(sigma1_p);
		} // else
			// std::cerr << "Ignoring t = " << ts[i] << std::endl;

		if (std::abs(std::imag(sigma2_m)) < 1e-4) {
			sigma2_m = std::real(sigma2_m);
			sigma2_p = std::real(sigma2_p);
		} // else
			// std::cerr << "Ignoring t = " << ts[i] << std::endl;

		//% Now check to see which pair pass. Only a single pair should pass, in theory.
		//% If not, issue a warning.
		static constexpr T my_eps = 1.0;

		if (std::abs(H + J*sigma1_m + K*sigma2_m + L*sigma1_m*sigma2_m) < my_eps) {
			sigmas1[i][sigmas1_len[i]] = sigma1_m.real();
			sigmas2[i][sigmas1_len[i]++] = sigma2_m.real();
		}
		if (std::abs(H + J*sigma1_p + K*sigma2_m + L*sigma1_p*sigma2_m) < my_eps) {
			// if (sigmas1_len[i] != 0) // !isempty(sigmas1[i])
		  //		std::cerr << "more than one sigma1, sigma2 pair satisfies the 3rd constraint" << std::endl;
			sigmas1[i][sigmas1_len[i]] = sigma1_p.real();
			sigmas2[i][sigmas1_len[i]++] = sigma2_m.real();
		}
		if (std::abs(H + J*sigma1_p + K*sigma2_p + L*sigma1_p*sigma2_p) < my_eps) {
			// if (sigmas1_len[i] != 0)
      //	std::cerr << "more than one sigma1, sigma2 pair satisfies the 3rd constraint" << std::endl;
			sigmas1[i][sigmas1_len[i]] = sigma1_p.real();
			sigmas2[i][sigmas1_len[i]++] = sigma2_p.real();
		}
		if (std::abs(H + J*sigma1_m + K*sigma2_p + L*sigma1_m*sigma2_p) < my_eps) {
			// if (sigmas1_len[i] != 0)
			// std::cerr << "more than one sigma1, sigma2 pair satisfies the 3rd constraint" << std::endl;
			sigmas1[i][sigmas1_len[i]] = sigma1_m.real();
			sigmas2[i][sigmas1_len[i]++] = sigma2_p.real();
		}
		// if (sigmas1_len[i] == 0) // isempty(sigmas1[i])
		// std::cerr << "no sigma1, sigma2 pair satisfies the 3rd constraint" << std::endl;
	}
}

template<typename T>
inline void
invm3x3(const T (&input_m)[3][3], T (&output_m)[3][3])
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
	a = input_m[0][0], b = input_m[0][1], c = input_m[0][2],
	d = input_m[1][0], e = input_m[1][1], f = input_m[1][2],
	g = input_m[2][0], h = input_m[2][1], i = input_m[2][2];

	const T 
	A =  (e*i - f*h), B = -(d*i - f*g), C =  (d*h - e*g),
	D = -(b*i - c*h), E =  (a*i - c*g), F = -(a*h - b*g),
	G =  (b*f - c*e), H = -(a*f - c*d), I =  (a*e - b*d);

	const T invdet_M = 1. / (a*A + b*B + c*C);
	output_m[0][0] = invdet_M * A; output_m[0][1] = invdet_M * D; output_m[0][2] = invdet_M * G;
	output_m[1][0] = invdet_M * B; output_m[1][1] = invdet_M * E; output_m[1][2] = invdet_M * H;
	output_m[2][0] = invdet_M * C; output_m[2][1] = invdet_M * F; output_m[2][2] = invdet_M * I;
}

template<typename T>
inline void
multm3x3(const T (&m1)[3][3], const T (&m2)[3][3], T output_m[][3])
{
	for (unsigned i = 0; i < 3; i++)
		for (unsigned j = 0; j < 3; j++) {
			output_m[i][j] = 0;
			for (unsigned k = 0; k < 3; k++)
				output_m[i][j] += m1[i][k] * m2[k][j];
		}
}

template<typename T>
void pose_poly<T>::
get_r_t_from_rhos(
	const unsigned ts_len,
	const T sigmas1[TS_MAX_LEN][TS_MAX_LEN], const unsigned sigmas1_len[TS_MAX_LEN],
	const T sigmas2[TS_MAX_LEN][TS_MAX_LEN],
	const T rhos1[ROOT_IDS_LEN], const T rhos2[ROOT_IDS_LEN],
	const T gama1[3], const T tgt1[3],
	const T gama2[3], const T tgt2[3],
	const T Gama1[3], const T Tgt1[3],
	const T Gama2[3], const T Tgt2[3],
	T (*output)[RT_MAX_LEN][4][3], unsigned *output_len
)
{
	T lambdas1[TS_MAX_LEN][TS_MAX_LEN]; T lambdas2[TS_MAX_LEN][TS_MAX_LEN];
	const T DGama[3] = {Gama1[0]-Gama2[0], Gama1[1]-Gama2[1], Gama1[2]-Gama2[2]};
  
	for (unsigned i = 0; i < ts_len; i++) {
    const T dgamas_rhos[3] = {
     rhos1[i]*gama1[0] - rhos2[i]*gama2[0],
     rhos1[i]*gama1[1] - rhos2[i]*gama2[1],
     rhos1[i]*gama1[2] - rhos2[i]*gama2[2]};
		for (unsigned j = 0; j < sigmas1_len[i]; j++) {
			lambdas1[i][j] = 
        (DGama[0]*Tgt1[0]+DGama[1]*Tgt1[1] + DGama[2]*Tgt1[2]) / 
        (dgamas_rhos[0]*(rhos1[i]*tgt1[0] + sigmas1[i][j]*gama1[0]) + 
        dgamas_rhos[1]*(rhos1[i]*tgt1[1] + sigmas1[i][j]*gama1[1]) +
        dgamas_rhos[2]*(rhos1[i]*tgt1[2] + sigmas1[i][j]*gama1[2]));
			lambdas2[i][j] = 
        (DGama[0]*Tgt2[0]+DGama[1]*Tgt2[1] + DGama[2]*Tgt2[2]) /
        (dgamas_rhos[0]*(rhos2[i]*tgt2[0] + sigmas2[i][j]*gama2[0]) + 
        dgamas_rhos[1]*(rhos2[i]*tgt2[1] + sigmas2[i][j]*gama2[1]) +
        dgamas_rhos[2]*(rhos2[i]*tgt2[2] + sigmas2[i][j]*gama2[2]));
		}
	}

	//% Rotation:
	const T A[3][3] = {
		DGama[0], Tgt1[0], Tgt2[0],
		DGama[1], Tgt1[1], Tgt2[1],
		DGama[2], Tgt1[2], Tgt2[2]};
	T inv_A[3][3]; invm3x3(A, inv_A);

	// Matrix containing Rotations and Translations
	T (&RT)[RT_MAX_LEN][4][3] = *output;
	unsigned &RT_len               = *output_len; RT_len = 0;
	for (unsigned i = 0; i < ts_len; i++) {
		for (unsigned j = 0; j < sigmas1_len[i]; j++, RT_len++) {
			T (&Rots)[4][3] = RT[RT_len]; T (&Transls)[3] = RT[RT_len][3];

			#define B_row(r) \
				rhos1[i]*gama1[(r)] - rhos2[i]*gama2[(r)], \
				lambdas1[i][j]*(rhos1[i]*tgt1[(r)] + sigmas1[i][j]*gama1[(r)]), \
				lambdas2[i][j]*(rhos2[i]*tgt2[(r)] + sigmas2[i][j]*gama2[(r)])

			const T B[3][3] = { B_row(0), B_row(1), B_row(2) };
			multm3x3(B, inv_A, Rots);
      Transls[0] = rhos1[i]*gama1[0] - Rots[0][0] * Gama1[0] - Rots[0][1] * Gama1[1] - Rots[0][2] * Gama1[2];
      Transls[1] = rhos1[i]*gama1[1] - Rots[1][0] * Gama1[0] - Rots[1][1] * Gama1[1] - Rots[1][2] * Gama1[2];
      Transls[2] = rhos1[i]*gama1[2] - Rots[2][0] * Gama1[0] - Rots[2][1] * Gama1[1] - Rots[2][2] * Gama1[2];
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
  assert(3 == bearing_vectors.rows());
  assert(3 == X.rows());
  assert(bearing_vectors.cols() == X.cols());
  assert(bearing_vectors.cols() == 2);

  unsigned nsols;
  double degen;
	double rotation_translation_solutions[RT_MAX_LEN][4][3];

  p2pt<double>::pose_from_point_tangents(
    bearing_vectors.col(0).data(), tangent_vectors.col(0).data(),
    bearing_vectors.col(1).data(), tangent_vectors.col(1).data(),
    X.col(0).data(), T.col(0).data(),
    X.col(1).data(), T.col(1).data(),
    &rotation_translation_solutions, &nsols, &degen
  );

	for (unsigned i = 0; i < nsols; ++i) {
    Mat34 P;
    for (unsigned j = 0 ; j < 3; ++j)
      for (unsigned k = 0 ; k < 3; ++k)
        P(j,k) = rotation_translation_solutions[i][j][k];

    for (unsigned k = 0 ; k < 3; ++k)
      P(k,3) = rotation_translation_solutions[i][3][k];
    models->push_back(P);
  }
};

} // namespace euclidean_resection
} // namespace openMVG
