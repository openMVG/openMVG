// Copyright (c) 2014 Romuald PERROT, Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_FEATURES_AKAZE_H
#define OPENMVG_FEATURES_AKAZE_H

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

#include "openMVG/image/image.hpp"
#include "openMVG/numeric/numeric.h"
#include "openMVG/numeric/math_trait.hpp"

#include "openMVG/features/feature.hpp"
#include "openMVG/features/descriptor.hpp"

#include "openMVG/features/akaze/mldb_descriptor.hpp"
#include "openMVG/features/akaze/msurf_descriptor.hpp"

namespace openMVG {

struct AKAZEConfig
{
  AKAZEConfig():
    iNbOctave(4),
    iNbSlicePerOctave(4),
    fSigma0(1.6f),
    fThreshold(0.0008f),
    fDesc_factor(1.f)
  {
  }

  int iNbOctave; ///< Octave to process 
  int iNbSlicePerOctave; ///< Levels per octave 
  float fSigma0; ///< Initial sigma offset (used to suppress low level noise)
  float fThreshold;  ///< Hessian determinant threshold
  float fDesc_factor;   ///< Magnifier used to describe an interest point
};

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
  Image<float>
    cur,    ///< Current gaussian image
    Lx,     ///< Current x derivatives
    Ly,     ///< Current y derivatives
    Lhess;  ///< Current Determinant of Hessian
};

/* ************************************************************************* */
// AKAZE Class Declaration
class AKAZE {

private:

  AKAZEConfig options_;               ///< Configuration options for AKAZE
  std::vector<TEvolution> evolution_;	///< Vector of nonlinear diffusion evolution (Scale Space)
  Image<float> in_;                  ///< Input image

public:

  /// Constructor
  AKAZE(const Image<unsigned char> & in, const AKAZEConfig & options);

  /// Compute the AKAZE non linear diffusion scale space per slice
  void Compute_AKAZEScaleSpace(void);

  /// Detect AKAZE feature in the AKAZE scale space
  void Feature_Detection(std::vector<AKAZEKeypoint>& kpts) const;

  /// Sub pixel refinement of the detected keypoints
  void Do_Subpixel_Refinement(std::vector<AKAZEKeypoint>& kpts) const;
  
  /// Sub pixel refinement of a keypoint
  bool Do_Subpixel_Refinement( AKAZEKeypoint & kpts, const Image<float> & Ldet) const;

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
    const Image<float> & Lx,
    const Image<float> & Ly) const;
  
  /// Compute an AKAZE slice    
  static
  void ComputeAKAZESlice(
    const Image<float> & src , // Input image for the given octave
    const int p , // octave index
    const int q , // slice index
    const int nbSlice , // slices per octave
    const float sigma0 , // first octave initial scale
    const float contrast_factor ,
    Image<float> & Li , // Diffusion image
    Image<float> & Lx , // X derivatives
    Image<float> & Ly , // Y derivatives
    Image<float> & Lhess // Det(Hessian)
    );

  /// Compute Contrast Factor
  static float ComputeAutomaticContrastFactor(
    const Image<float> & src,
    const float percentile );
};
/* ************************************************************************* */

/*
 * AKAZE keypoint extraction and description:
 * - MSURF: type = float; N = 64
 * - MLDB:  type = bool;  N = 486
 */
template<typename type, int N>
static bool AKAZEDetector(
  const Image<unsigned char>& I,
  std::vector<SIOPointFeature>& vec_feat,
  std::vector<Descriptor<type,N> >& vec_desc,
  AKAZEConfig options = AKAZEConfig())
{
  std::cerr << "You cannot compute the provided descriptor on AKAZE keypoints." << std::endl;
  return false;
}

/// Compute AKAZE keypoints and their associated MSURF descriptors
template<>
#if (WIN32)
static 
#endif //WIN32
bool AKAZEDetector<float,64>(
	const Image<unsigned char>& I,
	std::vector<SIOPointFeature>& vec_feat,
	std::vector<Descriptor<float,64> >& vec_desc,
	AKAZEConfig options)
{
	options.fDesc_factor = 10.f*sqrtf(2.f);

	AKAZE akaze(I, options);
	akaze.Compute_AKAZEScaleSpace();
	std::vector<AKAZEKeypoint> kpts;
	kpts.reserve(5000);
	akaze.Feature_Detection(kpts);
	akaze.Do_Subpixel_Refinement(kpts);

	vec_feat.resize(kpts.size());
	vec_desc.resize(kpts.size());
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
	for (int i = 0; i < static_cast<int>(kpts.size()); ++i)
	{
		AKAZEKeypoint ptAkaze = kpts[i];
		const TEvolution & cur_slice = akaze.getSlices()[ ptAkaze.class_id ] ;

		// Compute point orientation
		akaze.Compute_Main_Orientation(ptAkaze, cur_slice.Lx, cur_slice.Ly);

		vec_feat[i] = SIOPointFeature(ptAkaze.x, ptAkaze.y, ptAkaze.size, ptAkaze.angle);

		// Compute descriptor (MSURF)
		ComputeMSURFDescriptor( cur_slice.Lx , cur_slice.Ly , ptAkaze.octave , vec_feat[ i ] , vec_desc[ i ] ) ;
	}
	return true;
}

/// Compute AKAZE keypoints and their associated MLDB descriptors
template<>
#if (WIN32)
static 
#endif //WIN32
bool AKAZEDetector<bool,486>(
	const Image<unsigned char>& I,
	std::vector<SIOPointFeature>& vec_feat,
	std::vector<Descriptor<bool,486> >& vec_desc,
	AKAZEConfig options)
{
	options.fDesc_factor = 10.f*sqrtf(2.f);

	AKAZE akaze(I, options);
	akaze.Compute_AKAZEScaleSpace();
	std::vector<AKAZEKeypoint> kpts;
	kpts.reserve(5000);
	akaze.Feature_Detection(kpts);
	akaze.Do_Subpixel_Refinement(kpts);

	vec_feat.resize(kpts.size());
	vec_desc.resize(kpts.size());
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
	for (int i = 0; i < static_cast<int>(kpts.size()); ++i)
	{
		AKAZEKeypoint ptAkaze = kpts[i];
		const TEvolution & cur_slice = akaze.getSlices()[ ptAkaze.class_id ] ;

		// Compute point orientation
		akaze.Compute_Main_Orientation(ptAkaze, cur_slice.Lx, cur_slice.Ly);

		vec_feat[i] = SIOPointFeature(ptAkaze.x, ptAkaze.y, ptAkaze.size, ptAkaze.angle);

		// Compute descriptor (FULL MLDB)
		ComputeMLDBDescriptor( cur_slice.cur , cur_slice.Lx , cur_slice.Ly , ptAkaze.octave , vec_feat[i] , vec_desc[i] ) ;
	}
	return true;
}

}; // namespace openMVG

#endif //OPENMVG_FEATURES_AKAZE_H
