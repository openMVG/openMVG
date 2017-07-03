// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
#ifndef OPENMVG_FEATURES_SIFT_SIFT_KEYPOINT_EXTRACTOR_HPP
#define OPENMVG_FEATURES_SIFT_SIFT_KEYPOINT_EXTRACTOR_HPP


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

#include <vector>

#include "openMVG/features/feature.hpp"
#include "openMVG/features/sift/hierarchical_gaussian_scale_space.hpp"
#include "openMVG/features/sift/sift_keypoint.hpp"
#include "openMVG/image/image_container.hpp"

namespace openMVG{
namespace features{
namespace sift{

struct SIFT_KeypointExtractor
{
  /**
  * @brief SIFT_KeypointExtractor constructor
  * @param peak_threshold Threshold on DoG operator
  * @param edge_threshold Threshold on the ratio of principal curvatures
  * @param nb_refinement_step Maximum number of refinement step to find exact location of interest point
  */
  SIFT_KeypointExtractor
  (
    float peak_threshold,     // i.e => 0.04 / slices per octave
    float edge_threshold,     // i.e => 10
    int nb_refinement_step = 5
  ):
  m_peak_threshold(peak_threshold),
  m_edge_threshold(edge_threshold),
  m_nb_refinement_step(nb_refinement_step)
  {
  }

  /**
  * @brief Detect Scale Invariant points using Difference of Gaussians
  * @param octave A Gaussian octave
  * @param[out] keypoints The found Scale Invariant keypoint
  */
  void operator()( const Octave & octave , std::vector< Keypoint > & keypoints )
  {
    if (!ComputeDogs(octave))
      return;
    Find_3d_discrete_extrema(keypoints, 0.8f);
    Keypoints_refine_position(keypoints);
  }

protected:
  /**
  * @brief Compute the Difference of Gaussians (Dogs) for a Gaussian octave
  * @param Octave The input Gaussian octave
  * @retval true If Dogs have been computed
  * @retval false if Dogs cannot be computed (less than 2 images)
  */
  bool ComputeDogs(const Octave & octave)
  {
    const int n = octave.slices.size();
    if (n < 2)
      return false;

    m_Dogs.slices.resize(n-1);
    m_Dogs.octave_level = octave.octave_level;
    m_Dogs.delta = octave.delta;
    m_Dogs.sigmas = octave.sigmas;
    for (int s = 0; s < m_Dogs.slices.size(); ++s)
    {
      const image::Image<float> &P = octave.slices[s+1];
      const image::Image<float> &M = octave.slices[s];
      m_Dogs.slices[s] = P - M;
    }
    return true;
  }

  /**
  * @brief Tell if a point is local maximum/minimum
  * @param slices A Dog octave
  * @param id_slice The "middle" index, the slice id
  * @param id_row The discrete y point position
  * @param id_col The discrete x point position
  * @retval true If the point is local maximum/minimum
  * @retval false If the point is not a local maximum/minimum
  */
  static inline bool is_local_min_max
  (
    const std::vector< image::Image<float> > & slices,
    const size_t id_slice,
    const size_t id_row,
    const size_t id_col
  )
  {
    const float pix_val = std::abs(slices[id_slice](id_row, id_col));

    const bool is_min_or_max =
      // Current slice
      (pix_val > std::abs(slices[id_slice](id_row - 1, id_col - 1))) &&
      (pix_val > std::abs(slices[id_slice](id_row - 1, id_col))) &&
      (pix_val > std::abs(slices[id_slice](id_row - 1, id_col + 1))) &&
      (pix_val > std::abs(slices[id_slice](id_row, id_col - 1))) &&
      (pix_val > std::abs(slices[id_slice](id_row, id_col + 1))) &&
      (pix_val > std::abs(slices[id_slice](id_row + 1, id_col - 1))) &&
      (pix_val > std::abs(slices[id_slice](id_row + 1, id_col))) &&
      (pix_val > std::abs(slices[id_slice](id_row + 1, id_col + 1))) &&
      // Above slice
      (pix_val > std::abs(slices[id_slice - 1](id_row - 1, id_col - 1))) &&
      (pix_val > std::abs(slices[id_slice - 1](id_row - 1, id_col))) &&
      (pix_val > std::abs(slices[id_slice - 1](id_row - 1, id_col + 1))) &&
      (pix_val > std::abs(slices[id_slice - 1](id_row, id_col - 1))) &&
      (pix_val > std::abs(slices[id_slice - 1](id_row, id_col))) &&
      (pix_val > std::abs(slices[id_slice - 1](id_row, id_col + 1))) &&
      (pix_val > std::abs(slices[id_slice - 1](id_row + 1, id_col - 1))) &&
      (pix_val > std::abs(slices[id_slice - 1](id_row + 1, id_col))) &&
      (pix_val > std::abs(slices[id_slice - 1](id_row + 1, id_col + 1))) &&
      // Bottom slice
      (pix_val > std::abs(slices[id_slice + 1](id_row - 1, id_col - 1))) &&
      (pix_val > std::abs(slices[id_slice + 1](id_row - 1, id_col))) &&
      (pix_val > std::abs(slices[id_slice + 1](id_row - 1, id_col + 1))) &&
      (pix_val > std::abs(slices[id_slice + 1](id_row, id_col - 1))) &&
      (pix_val > std::abs(slices[id_slice + 1](id_row, id_col))) &&
      (pix_val > std::abs(slices[id_slice + 1](id_row, id_col + 1))) &&
      (pix_val > std::abs(slices[id_slice + 1](id_row + 1, id_col - 1))) &&
      (pix_val > std::abs(slices[id_slice + 1](id_row + 1, id_col))) &&
      (pix_val > std::abs(slices[id_slice + 1](id_row + 1, id_col + 1)));

    return is_min_or_max;
  }


  /**
  * @brief Compute the 2D Hessian response of the DoG operator is computed via finite difference schemes
  * @param key A Keypoint (the field edgeResp will be updated)
  * @retval the Harris and Stephen Edge response value
  */
  float Compute_edge_response
  (
    Keypoint & key
  ) const
  {
    const int s = key.s;
    const int i = key.i; // i = id_row
    const int j = key.j; // j = id_col
    const Octave & octave = m_Dogs;
    const image::Image<float> & im = octave.slices[s];
    // Compute the 2d Hessian at pixel (i,j)
    const float hXX = im(j,i-1) + im(j,i+1) - 2.f * im(j,i);
    const float hYY = im(j+1,i) + im(j-1,i) - 2.f * im(j,i);
    const float hXY = ((im(j+1,i+1) - im(j-1,i+1)) - ((im(j+1,i-1) - im(j-1,i-1))))/4.f;
    // Harris and Stephen Edge response
    key.edgeResp = Square(hXX + hYY)/(hXX*hYY - hXY*hXY);
    return key.edgeResp;
  }

  /**
  * @brief Find discrete extrema position (position, scale) in the Dog domain
  * @param[out] keypoints The list of found extrema as Keypoints
  * @param percent Percentage applied on of the internal Edge threshold value
  */
  void Find_3d_discrete_extrema
  (
    std::vector< Keypoint > & keypoints,
    float percent = 1.0f
  ) const
  {
    const int ns = m_Dogs.slices.size();
    const float delta = m_Dogs.delta;
    const int h = m_Dogs.slices[0].Height();
    const int w = m_Dogs.slices[0].Width();

    // Loop through the slices of the image stack (one octave)
    for (int s = 1; s < ns-1; ++s)
    {
      for (int id_row = 1; id_row < h-1; ++id_row )
      {
        for (int id_col = 1; id_col < w-1; ++id_col )
        {
          const float pix_val = m_Dogs.slices[s](id_row, id_col);
          if (std::abs(pix_val) > m_peak_threshold * percent)
          if (is_local_min_max(m_Dogs.slices, s, id_row, id_col))
          {
            // if 3d discrete extrema, save a candidate keypoint
            Keypoint key;
            key.i = id_col;
            key.j = id_row;
            key.s = s;
            key.o = m_Dogs.octave_level;
            key.x = delta * id_col;
            key.y = delta * id_row;
            key.sigma = m_Dogs.sigmas[s];
            key.val = pix_val;
            keypoints.emplace_back(key);
          }
        }
      }
    }
    keypoints.shrink_to_fit();
  }


  /**
  * @brief Refine the 3D location of a Keypoint using the local Hessian value (discrete to subpixel)
  * @param stack The list of found extrema as Keypoints
  * @param i Input discrete x location of the keypoint
  * @param j Input discrete y location of the keypoint
  * @param s Input scale of the keypoint
  * @param di Refined subpixel x location
  * @param dj Refined subpixel y location
  * @param ds Refined scale value
  * @param val Refined local intensity of the keypoint (Peak value)
  * @param ofstMax Upper limit of variation of the refined parameters
  * @retval true if the value have been refined
  * @retval false if the value cannot be refined (inside the provide range)
  */
  static bool Inverse_3D_Taylor_second_order_expansion
  (
    const Octave & stack, // the dog stack
    int i, int j, int s,
    float *di, float *dj, float *ds, float *val,
    const float ofstMax
  )
  {
    float hXX,hXY,hXS,hYY,hYS,hSS;
    float gX,gY,gS;
    float ofstX, ofstY, ofstS, ofstVal;

    const image::Image<float> & slice = stack.slices[s];
    const image::Image<float> & sliceU = stack.slices[s+1];
    const image::Image<float> & sliceD = stack.slices[s-1];

    // Compute the 3d Hessian at pixel (i,j,s)  Finite difference scheme
    hXX = slice(j,i-1) + slice(j,i+1) - 2.f*slice(j,i);
    hYY = slice(j+1,i) + slice(j-1,i) - 2.f*slice(j,i);
    hSS = sliceU(j,i)  + sliceD(j,i)  - 2.f*slice(j,i);
    hXY = (  (slice(j+1,i+1) - slice(j-1,i+1))
           - (slice(j+1,i-1) - slice(j-1,i-1)) )/4.f;
    hXS = (  (sliceU(j,i+1)  - sliceU(j,i-1))
           - (sliceD(j,i+1)  - sliceD(j,i-1)) )/4.f;
    hYS = (  (sliceU(j+1,i)  - sliceU(j-1,i))
           - (sliceD(j+1,i)  - sliceD(j-1,i)) )/4.f;

    // Compute the 3d gradient at pixel (i,j,s)
    gX = ( slice(j,i+1) - slice(j,i-1) ) / 2.f;
    gY = ( slice(j+1,i) - slice(j-1,i) ) / 2.f;
    gS = ( sliceU(j,i)  - sliceD(j,i)  ) / 2.f;

    // Inverse the Hessian - Fitting a quadratic function
    Eigen::Matrix<float, 3, 3> A;
    Vec3f b;
    A << hXX, hXY, hXS, hXY, hYY, hYS, hXS, hYS, hSS;
    b << -gX, -gY, -gS;

    // solve for offset
    Eigen::FullPivLU<Eigen::Matrix<float, 3, 3> > lu_decomp(A);
    if (!lu_decomp.isInvertible())
      return false;

    const Vec3f dst = lu_decomp.solve(b);
    ofstX = dst(0);
    ofstY = dst(1);
    ofstS = dst(2);
    // Compute the DoG value offset
    ofstVal = (gX * ofstX + gY * ofstY + gS * ofstS) / 2.f;

    // output
    *di = ofstX;
    *dj = ofstY;
    *ds = ofstS;
    *val = slice(j,i) + ofstVal;

    // return true is the quadratic model is consistent (to the given range)
    return std::abs(ofstX) < ofstMax && std::abs(ofstY) < ofstMax && std::abs(ofstS) < ofstMax;
  }

  /**
  * @brief Tell if a keypoint is close to the border according it's scale
  * @param key The list of found extrema as Keypoints
  * @param w Slice image width
  * @param h Slice image height
  * @param lambda Scaling factor apply on the keypoint scale
  * @retval true if the keypoint is in the [0+lambda*key.sigma; (w or h) - lambda*key.sigma] range
  * @retval false if the keypoint is too close the image border
  */
  inline bool Border_Check(const Keypoint & key, int w, int h, float lambda = 1.0f) const
  {
    const float x = key.x;
    const float y = key.y;
    const int octave = key.o;
    const float ratio = 1 << octave; //pow(2,p);
    const float sigma = key.sigma;
    const bool isIn = (x - lambda * sigma > 0.0 )&&( x + lambda * sigma < (float)w * ratio)
                   && (y - lambda * sigma > 0.0 )&&( y + lambda * sigma < (float)h * ratio);
    return isIn;
  }

  /**
  * @brief Refine the keypoint position (location in space and scale), discard keypoints that cannot be refined.
  * @param[in,out] key The list of refined keypoints
  */
  void Keypoints_refine_position
  (
    std::vector< Keypoint > & keypoints
  ) const
  {
    std::vector< Keypoint > kps;
    kps.reserve(keypoints.size());

    const float ofstMax = 0.6f;

    // Ratio between two consecutive scales in the slice
    const float sigma_ratio = m_Dogs.sigmas[1] / m_Dogs.sigmas[0];
    const float edge_thres = Square(m_edge_threshold + 1) / m_edge_threshold;

    const Octave & octave = m_Dogs;
    const int w = octave.slices[0].Width();
    const int h = octave.slices[0].Height();
    const float delta  = octave.delta;

    for (const auto & key : keypoints)
    {
      float val = key.val;

      int ic = key.i; // current discrete value of x coordinate - at each interpolation
      int jc = key.j; // current discrete value of y coordinate - at each interpolation
      int sc = key.s; // current discrete value of s coordinate - at each interpolation

      float ofstX = 0.0f;
      float ofstY = 0.0f;
      float ofstS = 0.0f;

      bool isConv = false;
      // While position cannot be refined and the number of refinement is not too much
      for ( int nIntrp = 0; !isConv && nIntrp < m_nb_refinement_step; ++nIntrp)
      {
        // Extrema interpolation via a quadratic function
        //   only if the detection is not too close to the border (so the discrete 3D Hessian is well defined)
        if ( 0 < ic &&  ic < (w-1) && 0 < jc && jc < (h-1) )
        {
          if (Inverse_3D_Taylor_second_order_expansion(octave, ic, jc, sc, &ofstX, &ofstY, &ofstS, &val, ofstMax))
            isConv = true;
        }
        else
        {
          isConv = false;
        }
        if (!isConv)
        {
          // let's explore the neighbourhood in
          // space...
          if (ofstX > +ofstMax && (ic+1) < (w-1) ) {++ic;}
          if (ofstX < -ofstMax && (ic-1) >  0    ) {--ic;}
          if (ofstY > +ofstMax && (jc+1) < (h-1) ) {++jc;}
          if (ofstY < -ofstMax && (jc-1) >  0    ) {--jc;}
          // ... and scale.
          /*
          if (ofstS > +ofstMax && (sc+1) < (ns-1)) {++sc;}
          if (ofstS < -ofstMax && (sc-1) >    0  ) {--sc;}
          */
        }
      }

      if (isConv)
      {
        // Peak threshold check
        if ( std::abs(val) > m_peak_threshold )
        {
          Keypoint kp = key;
          kp.x = (ic + ofstX) * delta;
          kp.y = (jc + ofstY) * delta;
          kp.i = ic;
          kp.j = jc;
          kp.s = sc;
          kp.sigma = octave.sigmas[sc] * pow(sigma_ratio, ofstS); // logarithmic scale
          kp.val = val;
          // Edge check
          if (Compute_edge_response(kp) >=0 && std::abs(kp.edgeResp) <= edge_thres)
          {
            // Border check
            if (Border_Check(kp, w, h))
            {
              kps.emplace_back(kp);
            }
          }
        }
      }
    }
    keypoints = std::move(kps);
    keypoints.shrink_to_fit();
  }

protected:
  Octave m_Dogs;

  // Keypoint detection parameters
  float m_peak_threshold;     // threshold on DoG operator
  float m_edge_threshold;    // threshold on the ratio of principal curvatures
  int m_nb_refinement_step; // Maximum number of refinement step to find exact location of interest point
};

} // namespace sift
} // namespace features
} // namespace openMVG

#endif // OPENMVG_FEATURES_SIFT_SIFT_KEYPOINT_EXTRACTOR_HPP
