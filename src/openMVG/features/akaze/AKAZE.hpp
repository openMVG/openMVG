// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2014 Romuald PERROT, Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_FEATURES_AKAZE_HPP
#define OPENMVG_FEATURES_AKAZE_HPP

#ifdef _MSC_VER
#pragma warning(once:4244)
#endif

//------------------
//-- Bibliography --
//------------------
//- [1] "Fast Explicit Diffusion for Accelerated Features in Nonlinear Scale Spaces."
//- Authors: Pablo F. Alcantarilla, Jesús Nuevo and Adrien Bartoli.
//- Date: September 2013.
//- Conference : BMVC.
//------
// Notes:
// This implementation is done for OpenMVG:
//  - code look the same but less memory are allocated,
//  - only Perona and Malik G2 diffusion (PM_G2) is implemented.
//------
// To see the reference implementation (OpenCV based):
// please visit https://github.com/pablofdezalc/akaze
// Authors: Pablo F. Alcantarilla (1), Jesus Nuevo (2)
// Institutions:
//  Toshiba Research Europe Ltd (1)
//  TrueVision Solutions (2)
//------

#include <vector>

#include "openMVG/image/image_container.hpp"

namespace openMVG {
namespace features {

struct AKAZEKeypoint{

  AKAZEKeypoint()
  {
    x = y = size = angle = response = 0.f;
    class_id = octave = 0;
  }

  float x,y;      //!< coordinates of the keypoints
  float size;     //!< diameter of the meaningful keypoint neighborhood
  float angle;    //!< computed orientation of the keypoint (-1 if not applicable);
                  //!< it's in [0,360) degrees and measured relative to
                  //!< image coordinate system, ie in clockwise.
  float response; //!< the response by which the most strong keypoints have been selected. Can be used for the further sorting or subsampling
  unsigned char octave; //!< octave (pyramid layer) from which the keypoint has been extracted
  unsigned char class_id; //!< object class (if the keypoints need to be clustered by an object they belong to)
};

struct TEvolution
{
  image::Image<float>
    cur,    ///< Current gaussian image
    Lx,     ///< Current x derivatives
    Ly,     ///< Current y derivatives
    Lhess;  ///< Current Determinant of Hessian
};

/* ************************************************************************* */
// AKAZE Class Declaration
class AKAZE {

public:
  struct Params
  {
    Params():
      iNbOctave(4),
      iNbSlicePerOctave(4),
      fSigma0(1.6f),
      fThreshold(0.0008f),
      fDesc_factor(1.f)
    {
    }

    template<class Archive>
    void serialize(Archive & ar);

    int iNbOctave; ///< Octave to process
    int iNbSlicePerOctave; ///< Levels per octave
    float fSigma0; ///< Initial sigma offset (used to suppress low level noise)
    float fThreshold;  ///< Hessian determinant threshold
    float fDesc_factor;   ///< Magnifier used to describe an interest point
  };

private:

  Params options_;               ///< Configuration options for AKAZE
  std::vector<TEvolution> evolution_;  ///< Vector of nonlinear diffusion evolution (Scale Space)
  image::Image<float> in_;            ///< Input image

public:

  /// Constructor
  AKAZE(const image::Image<unsigned char> & in, const Params & options);

  /// Compute the AKAZE non linear diffusion scale space per slice
  void Compute_AKAZEScaleSpace();

  /// Detect AKAZE feature in the AKAZE scale space
  void Feature_Detection(std::vector<AKAZEKeypoint>& kpts) const;

  /// Sub pixel refinement of the detected keypoints
  void Do_Subpixel_Refinement(std::vector<AKAZEKeypoint>& kpts) const;

  /// Sub pixel refinement of a keypoint
  bool Do_Subpixel_Refinement(AKAZEKeypoint & kpts, const image::Image<float> & Ldet) const;

  /// Scale Space accessor
  const std::vector<TEvolution> & getSlices() const {return evolution_;}

  /**
   * @brief This method computes the main orientation for a given keypoint
   * @param kpt Input keypoint
   * @note The orientation is computed using a similar approach as described in the
   * original SURF method. See Bay et al., Speeded Up Robust Features, ECCV 2006
  */
  void Compute_Main_Orientation(
    AKAZEKeypoint& kpt,
    const image::Image<float> & Lx,
    const image::Image<float> & Ly) const;

private:

  /// Compute an AKAZE slice
  static
  void ComputeAKAZESlice(
    const image::Image<float> & src, // Input image for the given octave
    const int p , // octave index
    const int q , // slice index
    const int nbSlice , // slices per octave
    const float sigma0 , // first octave initial scale
    const float contrast_factor ,
    image::Image<float> & Li, // Diffusion image
    image::Image<float> & Lx, // X derivatives
    image::Image<float> & Ly, // Y derivatives
    image::Image<float> & Lhess // Det(Hessian)
    );

  /// Compute Contrast Factor
  static float ComputeAutomaticContrastFactor(
    const image::Image<float> & src,
    const float percentile );
};
/* ************************************************************************* */

} // namespace features
} // namespace openMVG

#endif // OPENMVG_FEATURES_AKAZE_HPP
