//:\file
//\author Ricardo Fabbri, Brown & Rio de Janeiro State U. (rfabbri.github.io) 
//\date Tue Jun  1 15:03:36 -03 2021
//
// Misc. Utilities used by trifocal solver

namespace trifocal3pt {

//------------------------------------------------------------------------------
// Utilities
//------------------------------------------------------------------------------
void
revert_intrinsics(
    const double K[/*3 or 2 ignoring last line*/][3], 
    double pix_coords[2], 
    const double normalized_coords[2]);

void
revert_intrinsics_tgt(
    const double K[/*3 or 2 ignoring last line*/][3], 
    double pix_tgt_coords[2], 
    const double normalized_tgt_coords[2]);

void
invert_intrinsics(
    const double K[/*3 or 2 ignoring last line*/][3], 
    const double pix_coords[2], 
    double normalized_coords[2]);

void
invert_intrinsics_tgt(
    const double K[/*3 or 2 ignoring last line*/][3], 
    const double pix_tgt_coords[2], 
    double normalized_tgt_coords[2]);

} // namespace trifocal3pt
