/** @file mathop_sse2.c
 ** @brief mathop for SSE2 - Definition
 ** @author Andrea Vedaldi
 **/

/*
Copyright (C) 2007-12 Andrea Vedaldi and Brian Fulkerson.
All rights reserved.

This file is part of the VLFeat library and is made available under
the terms of the BSD license (see the COPYING file).
*/

/* ---------------------------------------------------------------- */
#ifndef VL_MATHOP_SSE2_INSTANTIATING
#define VL_MATHOP_SSE2_INSTANTIATING

#ifndef VL_DISABLE_SSE2
#ifndef __SSE2__
#  error "mathop_sse2.c must be compiled with SSE2 intrinsics enabled"
#endif

#include "generic.h"
#include "mathop.h"
#include "mathop_sse2.h"
#include <emmintrin.h>

#undef FLT
#define FLT VL_TYPE_DOUBLE
#include "mathop_sse2.c"

#undef FLT
#define FLT VL_TYPE_FLOAT
#include "mathop_sse2.c"

/* VL_DISABLE_SSE2 */
#endif

/* ---------------------------------------------------------------- */
/* VL_MATHOP_SSE2_INSTANTIATING */
#else

#include "float.th"

VL_INLINE T
VL_XCAT(_vl_vhsum_sse2_, SFX)(VTYPE x)
{
  T acc ;
#if (VSIZE == 4)
  {
    VTYPE sum ;
    VTYPE shuffle ;
    /* shuffle = [1 0 3 2] */
    /* sum     = [3+1 2+0 1+3 0+2] */
    /* shuffle = [2+0 3+1 0+2 1+3] */
    /* vacc    = [3+1+2+0 3+1+2+0 1+3+0+2 0+2+1+3] */
    shuffle = VSHU (x, x, _MM_SHUFFLE(1, 0, 3, 2)) ;
    sum     = VADD (x, shuffle) ;
    shuffle = VSHU (sum, sum, _MM_SHUFFLE(2, 3, 0, 1)) ;
    x       = VADD (sum, shuffle) ;
  }
#else
  {
    VTYPE shuffle ;
    /* acc     = [1   0  ] */
    /* shuffle = [0   1  ] */
    /* sum     = [1+0 0+1] */
    shuffle = VSHU (x, x, _MM_SHUFFLE2(0, 1)) ;
    x       = VADD (x, shuffle) ;
  }
#endif
  VST1(&acc, x);
  return acc ;
}

VL_EXPORT T
VL_XCAT(_vl_distance_l2_sse2_, SFX)
(vl_size dimension, T const * X, T const * Y)
{
  T const * X_end = X + dimension ;
  T const * X_vec_end = X_end - VSIZE + 1 ;
  T acc ;
  VTYPE vacc = VSTZ() ;
  vl_bool dataAligned = VALIGNED(X) & VALIGNED(Y) ;

  if (dataAligned) {
    while (X < X_vec_end) {
      VTYPE a = *(VTYPE*)X ;
      VTYPE b = *(VTYPE*)Y ;
      VTYPE delta = VSUB(a, b) ;
      VTYPE delta2 = VMUL(delta, delta) ;
      vacc = VADD(vacc, delta2) ;
      X += VSIZE ;
      Y += VSIZE ;
    }
  } else {
    while (X < X_vec_end) {
      VTYPE a = VLDU(X) ;
      VTYPE b = VLDU(Y) ;
      VTYPE delta = VSUB(a, b) ;
      VTYPE delta2 = VMUL(delta, delta) ;
      vacc = VADD(vacc, delta2) ;
      X += VSIZE ;
      Y += VSIZE ;
    }
  }

  acc = VL_XCAT(_vl_vhsum_sse2_, SFX)(vacc) ;

  while (X < X_end) {
    T a = *X++ ;
    T b = *Y++ ;
    T delta = a - b ;
    acc += delta * delta ;
  }

  return acc ;
}

VL_EXPORT T
VL_XCAT(_vl_distance_l1_sse2_, SFX)
(vl_size dimension, T const * X, T const * Y)
{
  T const * X_end = X + dimension ;
  T const * X_vec_end = X + dimension - VSIZE ;
  T acc ;
  VTYPE vacc = VSTZ() ;
  VTYPE vminus = VL_XCAT(_mm_set1_p, VSFX) ((T) -0.0) ; /* sign bit */
  vl_bool dataAligned = VALIGNED(X) & VALIGNED(Y) ;

  if (dataAligned) {
    while (X < X_vec_end) {
      VTYPE a = *(VTYPE*)X ;
      VTYPE b = *(VTYPE*)Y ;
      VTYPE delta = VSUB(a, b) ;
      vacc = VADD(vacc, VANDN(vminus, delta)) ;
      X += VSIZE ;
      Y += VSIZE ;
    }
  } else {
    while (X < X_vec_end) {
      VTYPE a = VLDU(X) ;
      VTYPE b = VLDU(Y) ;
      VTYPE delta = VSUB(a, b) ;
      vacc = VADD(vacc, VANDN(vminus, delta)) ;
      X += VSIZE ;
      Y += VSIZE ;
    }
  }

  acc = VL_XCAT(_vl_vhsum_sse2_, SFX)(vacc) ;

  while (X < X_end) {
    T a = *X++ ;
    T b = *Y++ ;
    T delta = a - b ;
    acc += VL_MAX(delta, - delta) ;
  }

  return acc ;
}

VL_EXPORT T
VL_XCAT(_vl_distance_chi2_sse2_, SFX)
(vl_size dimension, T const * X, T const * Y)
{
  T const * X_end = X + dimension ;
  T const * X_vec_end = X + dimension - VSIZE ;
  T acc ;
  VTYPE vacc = VSTZ() ;
  vl_bool dataAligned = VALIGNED(X) & VALIGNED(Y) ;

  if (dataAligned) {
    while (X < X_vec_end) {
      VTYPE a = *(VTYPE*)X ;
      VTYPE b = *(VTYPE*)Y ;
      VTYPE delta = VSUB(a, b) ;
      VTYPE denom = VADD(a, b) ;
      VTYPE numer = VMUL(delta, delta) ;
      VTYPE ratio = VDIV(numer, denom) ;
      ratio = VAND(ratio, VNEQ(denom, VSTZ())) ;
      vacc = VADD(vacc, ratio) ;
      X += VSIZE ;
      Y += VSIZE ;
    }
  } else {
    while (X < X_vec_end) {
      VTYPE a = VLDU(X) ;
      VTYPE b = VLDU(Y) ;
      VTYPE delta = VSUB(a, b) ;
      VTYPE denom = VADD(a, b) ;
      VTYPE numer = VMUL(delta, delta) ;
      VTYPE ratio = VDIV(numer, denom) ;
      ratio = VAND(ratio, VNEQ(denom, VSTZ())) ;
      vacc = VADD(vacc, ratio) ;
      X += VSIZE ;
      Y += VSIZE ;
    }
  }

  acc = VL_XCAT(_vl_vhsum_sse2_, SFX)(vacc) ;

  while (X < X_end) {
    T a = *X++ ;
    T b = *Y++ ;
    T delta = a - b ;
    T denom = a + b ;
    T numer = delta * delta ;
    if (denom) {
      T ratio = numer / denom ;
      acc += ratio ;
    }
  }
  return acc ;
}


VL_EXPORT T
VL_XCAT(_vl_kernel_l2_sse2_, SFX)
(vl_size dimension, T const * X, T const * Y)
{
  T const * X_end = X + dimension ;
  T const * X_vec_end = X_end - VSIZE + 1 ;
  T acc ;
  VTYPE vacc = VSTZ() ;
  vl_bool dataAligned = VALIGNED(X) & VALIGNED(Y) ;

  if (dataAligned) {
    while (X < X_vec_end) {
      VTYPE a = *(VTYPE*)X ;
      VTYPE b = *(VTYPE*)Y ;
      vacc = VADD(vacc, VMUL(a,b)) ;
      X += VSIZE ;
      Y += VSIZE ;
    }
  } else {
    while (X < X_vec_end) {
      VTYPE a = VLDU(X) ;
      VTYPE b = VLDU(Y) ;
      vacc = VADD(vacc, VMUL(a,b)) ;
      X += VSIZE ;
      Y += VSIZE ;
    }
  }

  acc = VL_XCAT(_vl_vhsum_sse2_, SFX)(vacc) ;

  while (X < X_end) {
    T a = *X++ ;
    T b = *Y++ ;
    acc += a * b ;
  }
  return acc ;
}

VL_EXPORT T
VL_XCAT(_vl_kernel_l1_sse2_, SFX)
(vl_size dimension, T const * X, T const * Y)
{
  T const * X_end = X + dimension ;
  T const * X_vec_end = X_end - VSIZE + 1 ;
  T acc ;
  VTYPE vacc = VSTZ() ;
  VTYPE vminus = VL_XCAT(_mm_set1_p, VSFX) ((T) -0.0) ;
  vl_bool dataAligned = VALIGNED(X) & VALIGNED(Y) ;

  if (dataAligned) {
    while (X < X_vec_end) {
      VTYPE a = *(VTYPE*)X ;
      VTYPE b = *(VTYPE*)Y ;
      VTYPE a_ = VANDN(vminus, a) ;
      VTYPE b_ = VANDN(vminus, b) ;
      VTYPE sum = VADD(a_,b_) ;
      VTYPE diff = VSUB(a, b) ;
      VTYPE diff_ = VANDN(vminus, diff) ;
      vacc = VADD(vacc, VSUB(sum, diff_)) ;
      X += VSIZE ;
      Y += VSIZE ;
    }
  } else {
    while (X < X_vec_end) {
      VTYPE a = VLDU(X) ;
      VTYPE b = VLDU(Y) ;
      VTYPE a_ = VANDN(vminus, a) ;
      VTYPE b_ = VANDN(vminus, b) ;
      VTYPE sum = VADD(a_,b_) ;
      VTYPE diff = VSUB(a, b) ;
      VTYPE diff_ = VANDN(vminus, diff) ;
      vacc = VADD(vacc, VSUB(sum, diff_)) ;
      X += VSIZE ;
      Y += VSIZE ;
    }
  }

  acc = VL_XCAT(_vl_vhsum_sse2_, SFX)(vacc) ;

  while (X < X_end) {
    T a = *X++ ;
    T b = *Y++ ;
    T a_ = VL_XCAT(vl_abs_, SFX) (a) ;
    T b_ = VL_XCAT(vl_abs_, SFX) (b) ;
    acc += a_ + b_ - VL_XCAT(vl_abs_, SFX) (a - b) ;
  }

  return acc / ((T)2) ;
}

VL_EXPORT T
VL_XCAT(_vl_kernel_chi2_sse2_, SFX)
(vl_size dimension, T const * X, T const * Y)
{
  T const * X_end = X + dimension ;
  T const * X_vec_end = X + dimension - VSIZE ;
  T acc ;
  VTYPE vacc = VSTZ() ;
  vl_bool dataAligned = VALIGNED(X) & VALIGNED(Y) ;

  if (dataAligned) {
    while (X < X_vec_end) {
      VTYPE a = *(VTYPE*)X ;
      VTYPE b = *(VTYPE*)Y ;
      VTYPE denom = VADD(a, b) ;
      VTYPE numer = VMUL(a,b) ;
      VTYPE ratio = VDIV(numer, denom) ;
      ratio = VAND(ratio, VNEQ(denom, VSTZ())) ;
      vacc = VADD(vacc, ratio) ;
      X += VSIZE ;
      Y += VSIZE ;
    }
  } else {
    while (X < X_vec_end) {
      VTYPE a = VLDU(X) ;
      VTYPE b = VLDU(Y) ;
      VTYPE denom = VADD(a, b) ;
      VTYPE numer = VMUL(a,b) ;
      VTYPE ratio = VDIV(numer, denom) ;
      ratio = VAND(ratio, VNEQ(denom, VSTZ())) ;
      vacc = VADD(vacc, ratio) ;
      X += VSIZE ;
      Y += VSIZE ;
    }
  }

  acc = VL_XCAT(_vl_vhsum_sse2_, SFX)(vacc) ;

  while (X < X_end) {
    T a = *X++ ;
    T b = *Y++ ;
    T denom = a + b ;
    if (denom) {
      T ratio = a * b / denom ;
      acc += ratio ;
    }
  }
  return ((T)2) * acc ;
}

/* VL_MATHOP_SSE2_INSTANTIATING */
#endif
