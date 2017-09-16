// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/*

== Patent Warning and License =================================================

The SIFT method is patented

    [2] "Method and apparatus for identifying scale invariant features
      in an image."
        David G. Lowe
        Patent number: 6711293
        Filing date: Mar 6, 2000
        Issue date: Mar 23, 2004
        Application number: 09/519,89

 These source codes are made available for the exclusive aim of serving as
 scientific tool to verify the soundness and completeness of the algorithm
 description. Compilation, execution and redistribution of this file may
 violate patents rights in certain countries. The situation being different
 for every country and changing over time, it is your responsibility to
 determine which patent rights restrictions apply to you before you compile,
 use, modify, or redistribute this file. A patent lawyer is qualified to make
 this determination. If and only if they don't conflict with any patent terms,
 you can benefit from the following license terms attached to this file.

The implementation is based on

    [1] "Anatomy of the SIFT Method."
        I. Rey Otero  and  M. Delbracio
        Image Processing Online, 2013.
        http://www.ipol.im/pub/algo/rd_anatomy_sift/
*/

#ifndef OPENMVG_FEATURES_SIFT_SIFT_DESCRIPTOR_EXTRACTOR_HPP
#define OPENMVG_FEATURES_SIFT_SIFT_DESCRIPTOR_EXTRACTOR_HPP

#include <algorithm>
#include <limits>
#include <vector>

#include "openMVG/features/feature.hpp"
#include "openMVG/features/sift/hierarchical_gaussian_scale_space.hpp"
#include "openMVG/features/sift/sift_keypoint.hpp"
#include "openMVG/image/image_container.hpp"

namespace openMVG{
namespace features{
namespace sift{

/**
 ** @brief Class managing SIFT descriptor computation (orientation and description)
 **/
class Sift_DescriptorExtractor
{
public:
  /**
   ** @brief Constructor
   ** @param nb_orientation_histogram_bin Number of orientation histogram bin to compute orientations of a keypoint
   ** @param m_orientation_scale Sets the width of the orientation patch
   ** @param nb_smoothing_step Number of smoothing step applied to orientation histogram before computing orientations of a keypoint
   ** @param percentage_keep Percentage of highest orientation histogram above which a new orientation is computed
   ** @param b_root_sift Compute the root sift normalized histogram
   ** @param descriptor_scale Scale of descriptor box computation
   ** @param nb_split2d Number of split (in x and in y) applied to descriptor histogram
   ** @param nb_split_angle Number of split (in z) applied to descriptor histogram
   ** @param clip_value Value above this threshold are clamped before computing final descriptor
   **/
  Sift_DescriptorExtractor
  (
    // Orientation computation parameters
    const int nb_orientation_histogram_bin = 36,
    const float orientation_scale = 1.5f,
    const int nb_smoothing_step = 6,
    const float percentage_keep = 0.8f,
    // Descriptor computation parameters
    const bool b_root_sift = true,
    const float descriptor_scale = 6.0f,
    const int nb_split2d = 4,
    const int nb_split_angle = 8,
    const float clip_value = 0.2f
  ):
  m_nb_orientation_histogram_bin(nb_orientation_histogram_bin),
  m_orientation_scale(orientation_scale),
  m_nb_smoothing_step(nb_smoothing_step),
  m_percentage_keep(percentage_keep),
  m_b_root_sift(b_root_sift),
  m_descriptor_scale(descriptor_scale),
  m_nb_split2d(nb_split2d),
  m_nb_split_angle(nb_split_angle),
  m_clip_value(clip_value)
  {
  }

  /**
  * @brief Compute the local orientation and descriptor of a list of given Keypoints
  * @param[in] octave A Gaussian Octave
  * @param[in,out] keypoints The list of found extrema as Keypoints with updated orientation and descriptor
  */
  void operator()
  (
    const Octave & octave ,
    std::vector< Keypoint > & keypoints
  )
  {
    Compute_Gradients(octave);
    Keypoints_orientations(keypoints);
    Keypoints_description(keypoints);
  }

  void Compute_Orientations
  (
    const Octave & octave,
    std::vector< Keypoint > & keypoints
  )
  {
    Compute_Gradients(octave);
    Keypoints_orientations(keypoints);
  }

protected:

  /**
  * @brief Compute the gradient of a given Gaussian Octave
  * @param[in] octave A gaussian Octave
  */
  void Compute_Gradients(const Octave & octave)
  {
    const int nSca = octave.slices.size();
    m_xgradient.slices.resize(nSca);
    m_ygradient.slices.resize(nSca);
    m_xgradient.delta = octave.delta;
    m_ygradient.delta = octave.delta;
    m_xgradient.octave_level = octave.octave_level;
    m_ygradient.octave_level = octave.octave_level;
    for (int s = 1; s < nSca-1; ++s)
    {
      // only in range [1; n-1] (since first and last images were only used for non max suppression)
      image::ImageXDerivative( octave.slices[s], m_xgradient.slices[s]);
      image::ImageYDerivative( octave.slices[s], m_ygradient.slices[s]);
    }
  }

  // Compute the x modulus y
  // Modulo than handle negative values (for circular way)
  static inline float modulus(float x, float y)
  {
    float z = x;
    if (z < 0)
    {
      z += ((int)( -z / y) + 1) * y;
    }
    z -= (int)(z / y) * y;
    return z;
  }

  static inline int ori_to_bin(float ori, int nbins)
  {
    if (ori < 0)
      ori += 2 * M_PI;
    return static_cast<int>((ori/(2*M_PI)*nbins+0.5))%nbins;
  }

  static inline float bin_to_ori(float bin, int nbins)
  {
    float ori = (bin+0.5)*2*M_PI/(float)nbins;
    if (ori > M_PI)
      ori -= 2*M_PI;
    return ori;
  }

  /**
  * @brief Accumulate gradient orientation histogram around a keypoint
  * @param key The keypoint for which orientation histogram
  * @param[out] hist The orientation histogram
  */
  void Keypoint_orientation_histogram
  (
    const Keypoint & key,
    openMVG::Vecf & hist) const
  {
    hist.resize(m_nb_orientation_histogram_bin);
    hist.fill(0.0f);

    // conversion to the octave's coordinates (sampling)
    const float delta  = m_xgradient.delta;
    const float x = key.x / delta;
    const float y = key.y / delta;
    const float sigma = key.sigma / delta;

    const int w = m_xgradient.slices[1].Width(); // 1 since 0 is empty
    const int h = m_xgradient.slices[1].Height();

    const image::Image<float> & xgradient = m_xgradient.slices[key.s];
    const image::Image<float> & ygradient = m_ygradient.slices[key.s];

    // Contributing pixels are inside a patch [siMin;siMax] X [sjMin;sjMax]
    // of width w 6*lambda_ori*sigma_key (=9*sigma_key)
    const float R = 3 * m_orientation_scale * sigma;
    const int siMin = std::max(0, (int)(x - R + 0.5));
    const int sjMin = std::max(0, (int)(y - R + 0.5));
    const int siMax = std::min((int)(x + R + 0.5), w-1);
    const int sjMax = std::min((int)(y + R + 0.5), h-1);

    // For each pixel neighbourhood patch
    for (int si = siMin; si <= siMax; ++si)
    {
      for (int sj = sjMin; sj <= sjMax; ++sj)
      {
        // gradient orientation (theta)
        const float dx = xgradient(sj, si);
        const float dy = ygradient(sj, si);
        const float ori = modulus(atan2(dy, dx), 2*M_PI);

        // gradient magnitude with Gaussian weighting
        const float sX = (si - x) / sigma;
        const float sY = (sj - y) / sigma;
        const float r2 = Square(sX) + Square(sY);
        const float M = hypot(dx,dy) * exp(-r2/(2*Square(m_orientation_scale)));

        /// Determine the bin index in the circular histogram
        const int gamma = ori_to_bin(ori, m_nb_orientation_histogram_bin);

        /// Add the contribution to the orientation histogram
        hist[gamma] += M;
      }
    }
  }

  /**
  * @brief Iterative box filter based of a triangular length 3 kernel
  * @param[in] iter The number of filtering operation
  * @param[in,out] orientation_histogram The histogram to smooth
  */
  static void smooth_circular_histogram
  (
    int niter,
    openMVG::Vecf & orientation_histogram
  )
  {
    int i,i_prev,i_next;
    openMVG::Vecf tmp = orientation_histogram;
    const int nbins = tmp.rows();
    // Convolution with box filters
    for (; niter > 0; --niter){
      tmp = orientation_histogram;
      for (i = 0; i < nbins; ++i){
        i_prev = (i-1+nbins)%nbins;
        i_next = (i+1)%nbins;
        orientation_histogram[i] = (tmp[i_prev]/4+tmp[i]/2+tmp[i_next]/4);
      }
    }
  }

  /**
  * @brief Quadratic interpolation
  * @retval The quadratic interpolated value
  */
  static float interpolate_peak(float h1, float h2, float h3)
  {
    return (h1-h3)/(2*(h1+h3-2*h2));
  }

  /**
  * @brief Extract the principal orientations from a gradient orientation histogram
  *        Local maxima are considered secondary
  *        principal orientation if they exceed
  *        threshold times the absolute maximum (expressed in percentage of the max value)
  * @param[in] orientation_histogram The orientation histogram
  * @param[out] principal_orientations the principal orientation value (as bin index)
  * @retval the number of found principal orientation
  */
  int Extract_principal_orientations
  (
    openMVG::Vecf & orientation_histogram,
    std::vector<float> & principal_orientations
  ) const
  {
    // number of principal orientations (the return value).
    int o = 0;

    // Smooth histogram : iterated box filters
    smooth_circular_histogram(m_nb_smoothing_step, orientation_histogram);

    // What is the value of the global maximum
    const float max_value = orientation_histogram.maxCoeff();
    // Search for local extrema in the histogram
    const int nbins = orientation_histogram.rows();
    for (int i = 0; i < nbins; ++i){
      const int i_prev = (i-1+nbins)%nbins;
      const int i_next = (i+1)%nbins;
      if ( (orientation_histogram[i] > m_percentage_keep * max_value)
        && (orientation_histogram[i]>orientation_histogram[i_prev])
        && (orientation_histogram[i]>orientation_histogram[i_next]))
      {
        // Quadratic interpolation of the position of each local maximum
        const float offset = interpolate_peak(orientation_histogram[i_prev], orientation_histogram[i], orientation_histogram[i_next]);
        // Add to vector of principal orientations (expressed in [0,2pi]
        principal_orientations[o] = bin_to_ori((float)i + offset, nbins);
        ++o;
      }
    }
    // return the number of principal orientations
    return o;
  }

  /**
  * @brief Compute the orientations of a list of Keypoints
  * @param[in,out] keypoints The list of found Keypoints
  */
  void Keypoints_orientations
  (
    std::vector<Keypoint> & keypoints
  ) const
  {
    std::vector<Keypoint> kps;
    kps.reserve(keypoints.size());
#ifdef OPENMVG_USE_OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for (int i_key = 0; i_key < keypoints.size(); ++i_key)
    {
      const Keypoint & key = keypoints[i_key];

      openMVG::Vecf orientation_histogram;
      Keypoint_orientation_histogram(key, orientation_histogram);

      // Compute principal orientation(s)
      std::vector<float> principal_orientations(m_nb_orientation_histogram_bin);
      const int n_prOri = Extract_principal_orientations(orientation_histogram, principal_orientations);

      // Updating keypoints and save it in the new list
#ifdef OPENMVG_USE_OPENMP
      #pragma omp critical
#endif
      for (int n = 0; n < n_prOri; ++n)
      {
        Keypoint kp = key;
        kp.theta = principal_orientations[n];
        kps.emplace_back(kp);
      }
    }
    keypoints = std::move(kps);
    keypoints.shrink_to_fit();
  }

  /// Fast approximation of atan2
  /// http://www.dspguru.com/dsp/tricks/fixed-point-atan2-with-self-normalization
  static inline float atan2_fast
  (
    float y,
    float x
  )
  {
    float angle, r;
    const float c3 = 0.1821f;
    const float c1 = 0.9675f;
    const float abs_y = std::abs(y) + std::numeric_limits<float>::min();

    if (x >= 0) {
      r = (x - abs_y) / (x + abs_y);
      angle = (float)(M_PI / 4);
    }
    else {
      r = (x + abs_y) / (abs_y - x);
      angle = (float)(3 * M_PI / 4);
    }
    angle += (c3 * r * r - c1) * r;
    return (y < 0) ? -angle : angle;
  }

  static inline float mod_2pi_fast
  (
    float x
  )
  {
    while (x >(float)(2 * M_PI)) x -= (float)(2 * M_PI);
    while (x < 0.0f) x += (float)(2 * M_PI);
    return x;
  }


  /**
  * @brief Extract a SIFT descriptor for a given Keypoint
  * @param[in] k The keypoint for which the descriptor is required
  * @param[out] descr The found SIFT descriptor
  */
  void Extract_sift_feature_vector
  (
    const Keypoint &k,
    openMVG::Vecf & descr
  ) const
  {
    descr.fill(0.0f);
    const float delta = m_xgradient.delta;
    const int w = m_xgradient.slices[1].Width();
    const int h = m_xgradient.slices[1].Height();
    const float x = k.x / delta;
    const float y = k.y / delta;
    const float sigma = k.sigma / delta;

    const float c = cos(-k.theta);
    const float s = sin(-k.theta);

    const image::Image<float> & xgradient = m_xgradient.slices[k.s];
    const image::Image<float> & ygradient = m_ygradient.slices[k.s];

    // Local patch size is [siMin;siMax]X[sjMin;sjMax] which width W
    // W = 2*lambda_descr*sigma_key*(m_nb_split2d+1)/m_nb_split2d
    const float R =  (1+1/(float)m_nb_split2d)* m_descriptor_scale * sigma;
    const float Rp =  1.41421356237309504880 * R;
    const int siMin = std::max(0, (int)(x - Rp + 0.5f));
    const int sjMin = std::max(0, (int)(y - Rp + 0.5f));
    const int siMax = std::min((int)(x + Rp + 0.5f), w-1);
    const int sjMax = std::min((int)(y + Rp + 0.5f), h-1);
    /// For each pixel inside the patch.
    for (int si = siMin; si < siMax; ++si)
    {
      for (int sj = sjMin; sj < sjMax; ++sj)
      {
        // Compute pixel coordinates (sX,sY) on keypoint's invariant referential.
        const float Xref = si - x;
        const float Yref = sj - y;
        const float X = c * Xref - s * Yref;
        const float Y = s * Xref + c * Yref;
        // Does this sample fall inside the descriptor area ?
        if (std::max(std::abs(X), std::abs(Y)) < R)
        {
          // Compute the gradient orientation (theta) on keypoint referential.
          const float dx = xgradient(sj, si);
          const float dy = ygradient(sj, si);
          const float atan2val = atan2_fast(dy, dx);
          const float ori = mod_2pi_fast(atan2val - k.theta);

          // Compute the gradient magnitude and apply a Gaussian weighting to give less emphasis to distant sample
          const float t = m_descriptor_scale * sigma;
          const float M = hypot(dx, dy) * exp(-(Square(X)+Square(Y))/(2*Square(t)));

          // bin indices, Compute the (tri)linear weightings ...
          const float alpha = X/(2 * m_descriptor_scale * sigma / m_nb_split2d) + (m_nb_split2d-1.0)/2.0;
          const float beta  = Y/(2 * m_descriptor_scale * sigma / m_nb_split2d) + (m_nb_split2d-1.0)/2.0;
          const float gamma = ori / (2*M_PI) * m_nb_split_angle;
          //    ...and add contributions to respective bins in different histograms.
          // a loop with 1 or two elements
          const int i0 = std::floor(alpha);
          const int j0 = std::floor(beta);
          for (int i = std::max(0,i0); i <= std::min(i0+1,m_nb_split2d-1); ++i)
          {
            for (int j = std::max(0,j0); j <= std::min(j0+1,m_nb_split2d-1); ++j)
            { // looping through all surrounding histograms.

              const int index = i*m_nb_split2d*m_nb_split_angle+j*m_nb_split_angle;
              // Contribution to left bin.
              int k = ((int)gamma + m_nb_split_angle) % m_nb_split_angle;
              descr[index+k]
                += (1.0f-(gamma-floor(gamma)))
                   *(1.0f-std::abs((float)i-alpha))
                   *(1.0f-std::abs((float)j-beta))
                   *M;

              // Contribution to right bin.
              k = ((int)gamma + 1 + m_nb_split_angle) % m_nb_split_angle;
              descr[index+k]
                += (1.0f-(floor(gamma)+1.f-gamma))
                  *(1.0f-std::abs((float)i-alpha))
                  *(1.0f-std::abs((float)j-beta))
                  *M;
            }
          }
        }
      }
    }
  }

  /**
  * @brief Compute the SIFT descriptor of a list of Keypoints, apply the requested normalization and discretization [0;255]
  * @param[in,out] keypoints The list of keypoints that must have their SIFT descriptor computed
  */
  void Keypoints_description
  (
    std::vector<Keypoint> & keypoints
  ) const
  {
    const int n_descr = Square(m_nb_split2d) * m_nb_split_angle;
#ifdef OPENMVG_USE_OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for (int i_key = 0; i_key < keypoints.size(); ++i_key)
    {
      Keypoint & key = keypoints[i_key];
      // Compute the SIFT descriptor
      key.descr.resize(n_descr);
      Extract_sift_feature_vector(key, key.descr);

      // Threshold bins
      float norm = (m_b_root_sift ? key.descr.lpNorm<1>() : key.descr.lpNorm<2>())
        + std::numeric_limits<float>::epsilon();
      for (int i = 0; i < key.descr.rows(); ++i){
        key.descr[i] = std::min(key.descr[i], m_clip_value * norm);
      }
      // Quantization
      if (m_b_root_sift)
      {
        // scaling(rootsift) = sqrt( sift / sum(sift) );
        norm = key.descr.lpNorm<1>() + std::numeric_limits<float>::epsilon();
        key.descr = 512.f * (key.descr / norm).array().sqrt();
      }
      else
      {
        // scaling(normalized descriptor)
        norm = key.descr.lpNorm<2>() + std::numeric_limits<float>::epsilon();
        key.descr = 512.f * (key.descr / norm).array();
      }
    }
  }
protected:

  Octave m_xgradient; // gradient along x direction
  Octave m_ygradient; // gradient along y direction

  // Orientation histogram params
  int m_nb_orientation_histogram_bin;   // Nb bin for computing orientation
  float m_orientation_scale;            // Sets the width of the orientation patch
  int m_nb_smoothing_step;              // Smoothing step before orientation assignment
  float m_percentage_keep;              // Percentage to keep an orientation (minimal value in the histogram)

  // Orientation & Descriptor extraction parameters
  bool m_b_root_sift;       // Apply or not root sift normalized histogram
  float m_descriptor_scale; // Scale of descriptor box computation (sets the width of the descriptor patch)
  int m_nb_split2d;         // Number of split (in x and in y) applied to descriptor histogram
  int m_nb_split_angle;     // Number of split (in z) applied to descriptor histogram
  float m_clip_value;       // Threshold value to clamp large descriptor peak
};

} // namespace sift
} // namespace features
} // namespace openMVG

#endif // OPENMVG_FEATURES_SIFT_SIFT_DESCRIPTOR_EXTRACTOR_HPP
