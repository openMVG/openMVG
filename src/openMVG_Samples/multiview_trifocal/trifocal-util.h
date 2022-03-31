//:\file
//\author Ricardo Fabbri, Brown & Rio de Janeiro State U. (rfabbri.github.io) 
//\date Tue Jun  1 15:03:36 -03 2021
//
// Misc. Utilities used by trifocal solver

namespace trifocal3pt {

//------------------------------------------------------------------------------
// Utilities
// TODO: most of these can be made inline
//------------------------------------------------------------------------------
void apply_intrinsics(
    const double K[/*3 or 2 ignoring last line*/][3], 
    double px_coords[2], 
    const double normalized_coords[2]);

void apply_intrinsics_tgt(
    const double K[/*3 or 2 ignoring last line*/][3], 
    double px_tgt_coords[2], 
    const double normalized_tgt_coords[2]);

// px_coords and normalized_coords can be the same vector (in-place)
void invert_intrinsics(
    const double K[/*3 or 2 ignoring last line*/][3], 
    const double px_coords[2], 
    double normalized_coords[2]);

// px_coords and normalized_coords can be the same vector (in-place)
void invert_intrinsics_tgt(
    const double K[/*3 or 2 ignoring last line*/][3], 
    const double px_tgt_coords[2], 
    double normalized_tgt_coords[2]);

} // namespace trifocal3pt
