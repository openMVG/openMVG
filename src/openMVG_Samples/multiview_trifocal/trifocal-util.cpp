//:\file
//\author Ricardo Fabbri, Brown & Rio de Janeiro State U. (rfabbri.github.io) 
//\date Tue Jun  1 15:03:36 -03 2021
//
// Misc. Utilities used by trifocal solver
#include "trifocal-util.h"

namespace trifocal3pt {

//------------------------------------------------------------------------------
// Utilities
//------------------------------------------------------------------------------

// NoteToSelf(gabriel) For intrinsic parameter transform, see big notes eq. 5.2.13 at beginning of the code.
void
revert_intrinsics(
    const double K[/*3 or 2 ignoring last line*/][3], 
    double px_coords[2], 
    const double normalized_coords[2])
{
  double *px = px_coords;
  const double *nrm = normalized_coords;
  // XXX: usar a inversa da formula exatamente como em invert_intrinsics.
  //      ter certeza que funciona se a entrada e saida forem mesmas posicoes de
  //      memoria
  px[0] = nrm[0]*K[0][0]+nrm[1]*K[0][1]+nrm[2]*K[0][2];
  px[1] = nrm[0]*K[1][0]+nrm[1]*K[1][1]+nrm[2]*K[1][2];
}

void
revert_intrinsics_tgt(
    const double K[/*3 or 2 ignoring last line*/][3], 
    double px_tgt_coords[2], 
    const double normalized_tgt_coords[2])
{
  double *tp = px_tgt_coords;
  const double *t = normalized_tgt_coords;
  tp[0] = t[0]*K[0][0]+t[1]*K[0][1]+t[2]*K[0][2];
  tp[1] = t[0]*K[1][0]+t[1]*K[1][1]+t[2]*K[1][2];
}

void
invert_intrinsics(
    const double K[/*3 or 2 ignoring last line*/][3], 
    const double px_coords[2], 
    double normalized_coords[2])
{
  const double *px = px_coords;
  double *nrm = normalized_coords;
  nrm[1] = (px[1] - K[1][2]) /K[1][1];
  nrm[0] = (px[0] - K[0][1]*nrm[1] - K[0][2])/K[0][0];
}

void
invert_intrinsics_tgt(
    const double K[/*3 or 2 ignoring last line*/][3], 
    const double px_tgt_coords[2], 
    double normalized_tgt_coords[2])
{
  const double *tp = px_tgt_coords;
  double *t = normalized_tgt_coords;
  t[1] = tp[1]/K[1][1];
  t[0] = (tp[0] - K[0][1]*tp[1])/K[0][0];
}
// See big notes eq. 5.2.13 at beginning of the code.
}
