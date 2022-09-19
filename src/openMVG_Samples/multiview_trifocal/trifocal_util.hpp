//:\file
//\author Ricardo Fabbri, Rio de Janeiro State U. (rfabbri.github.io) 
//\date Tue Jun  1 15:03:36 -03 2021
//\author Gabriel Andrade, Rio de Janeiro State U.
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

// Get a reasonable error threshold in normalized coordinates
//
// Take a (threshold,0) vector along the x axis and 
// 
// transform to normalized coordinates
// Currently ignores skew
// 
// TODO(better guess is possible)
inline double threshold_pixel_to_normalized(double threshold, const double K[2][3]) {
  return threshold/K[0][0];
}

inline double threshold_normalized_to_pixel(double threshold, const double K[2][3]) {
  return threshold*K[0][0];
}
} // namespace trifocal3pt
