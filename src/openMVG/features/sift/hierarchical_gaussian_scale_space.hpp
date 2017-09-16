// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON, Romuald PERROT.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_FEATURES_SIFT_HIERARCHICAL_GAUSSIAN_SCALE_SPACE_HPP
#define OPENMVG_FEATURES_SIFT_HIERARCHICAL_GAUSSIAN_SCALE_SPACE_HPP

#include <algorithm>
#include <vector>

#include "openMVG/features/sift/octaver.hpp"
#include "openMVG/image/image_filtering.hpp"
#include "openMVG/image/image_resampling.hpp"
#include "openMVG/numeric/numeric.h"

namespace openMVG{
namespace features{

struct Octave{
  int octave_level;                   // the octave level
  float delta;                        // sampling rate in this octave
  std::vector<float> sigmas;          // sigma values
  std::vector<image::Image<float> > slices;  // octave slice (from fine to coarse)

};

struct GaussianScaleSpaceParams
{
  float sigma_min; // minimal level of blur featured in the scale-space
  float delta_min; // minimal inter-pixel distance featured in the scale-space
  float sigma_in;  // assumed level of blur in the input image
  int supplementary_levels; //necessary for post processing (i.e DOG)

  /*

  Explanation for the supplementary_levels parameters:
   For some computation you need supplementary levels,
    to ensure manipulation accross the gaussian pyramid.
   i.e. SIFT: when maxima/minima are searched over three levels


                       Slice 0 => ex. 1 slice per octave was asked
                     /
               DOG 0
             /       \
           /           Slice 1
         /           /
 MAX 0   ----  DOG 1
         \           \
           \           Slice 2
             \       /
               DOG 2
                     \
                       Slice 3 => 3 supplementary level are required to produce 3 DoGs

    to resume : Max computation required n+2 Dog so it implies n+3 Slice
  */

  GaussianScaleSpaceParams(
    const float sigma_min_ = 1.6f,
    const float delta_min_ = 1.0f,
    const float sigma_in_  = 0.5f,
    const int supplementary_levels_ = 0
    ):
    sigma_min(sigma_min_),
    delta_min(delta_min_),
    sigma_in(sigma_in_),
    supplementary_levels(supplementary_levels_)
  {  }
};

struct HierarchicalGaussianScaleSpace: public Octaver<Octave>
{
  /**
  * @brief HierarchicalGaussianScaleSpace constructor
  * @param nb_octave Number of Octave
  * @param nb_slice Number of slice per octave
  * @param params Parameters of the Gaussian scale space
  */
  HierarchicalGaussianScaleSpace(
    const int nb_octave = 6,
    const int nb_slice = 3,
    const GaussianScaleSpaceParams & params =
      std::move(GaussianScaleSpaceParams())
  ) :Octaver<Octave>(nb_octave, nb_slice),
    m_params(params),
    m_cur_octave_id(0)
  {
  }

  /**
  * @brief Set Initial image and update nb_octave if necessary
  * @param img Input image
  */
  virtual void SetImage(const image::Image<float> & img)
  {
    const double sigma_extra =
      sqrt(Square(m_params.sigma_min) - Square(m_params.sigma_in)) / m_params.delta_min;
    if (m_params.delta_min == 1.0f)
    {
      image::ImageGaussianFilter(img, sigma_extra, m_cur_base_octave_image);
    }
    else  // delta_min == 1
    {
      if (m_params.delta_min == 0.5f)
      {
        image::Image<float> tmp;
        ImageUpsample(img, tmp);
        image::ImageGaussianFilter(tmp, sigma_extra, m_cur_base_octave_image);
      }
      else
      {
        std::cerr
          << "Upsampling or downsampling with delta equal to: "
          << m_params.delta_min << " is not yet implemented" << std::endl;
      }
    }
    //-- Limit the size of the last octave to be at least 32x32 pixels
    const int nbOctaveMax = std::ceil(std::log2( std::min(m_cur_base_octave_image.Width(), m_cur_base_octave_image.Height())/32));
    m_nb_octave = std::min(m_nb_octave, nbOctaveMax);
  }

  /**
  * @brief Compute a full octave
  * @param[out] oct Computed octave
  * @retval true If an octave have been computed
  * @retval false if no octave have been computed (process ended)
  */
  virtual bool NextOctave(Octave & octave)
  {
    if (m_cur_octave_id >= m_nb_octave)
    {
      return false;
    }
    else
    {
      octave.octave_level = m_cur_octave_id;
      if (m_cur_octave_id == 0)
      {
        octave.delta = m_params.delta_min;
      }
      else
      {
        octave.delta *= 2.0f;
      }

      // init the "blur"/sigma scale spaces values
      octave.slices.resize(m_nb_slice + m_params.supplementary_levels);
      octave.sigmas.resize(m_nb_slice + m_params.supplementary_levels);
      for (int s = 0; s < m_nb_slice  + m_params.supplementary_levels; ++s)
      {
        octave.sigmas[s] =
          octave.delta / m_params.delta_min * m_params.sigma_min * pow(2.0,(float)s/(float)m_nb_slice);
      }

      // Build the octave iteratively
      octave.slices[0] = m_cur_base_octave_image;
      for (int s = 1; s < octave.sigmas.size(); ++s)
      {
        // Iterative blurring the previous image
        const image::Image<float> & im_prev = octave.slices[s-1];
        image::Image<float> & im_next = octave.slices[s];
        const double sig_prev = octave.sigmas[s-1];
        const double sig_next = octave.sigmas[s];
        const double sigma_extra = sqrt(Square(sig_next) - Square(sig_prev)) / octave.delta;

        image::ImageGaussianFilter(im_prev, sigma_extra, im_next);
      }
      /*
      // Debug: Export DoG scale space on disk
      for (int s = 0; s < octave.sigmas.size(); ++s)
      {
        std::stringstream os;
        os << "DoG_out_00" << m_cur_octave_id << "_s" << "00" << s << ".png";
        image::WriteImage(os.str().c_str(), octave.slices[s]);
      }
      */

      // Prepare for next octave computation -> Decimation
      ++m_cur_octave_id;
      if (m_cur_octave_id < m_nb_octave)
      {
        // Decimate => sigma * 2 for the next iteration
        const int index = (m_params.supplementary_levels == 0) ? 1 : m_params.supplementary_levels;
        ImageDecimate(octave.slices[octave.sigmas.size()-index], m_cur_base_octave_image);
      }
      return true;
    }
  }

protected:
  GaussianScaleSpaceParams m_params;  // The Gaussian scale space parameters
  image::Image<float> m_cur_base_octave_image; // The image that will be used to generate the next octave
  int m_cur_octave_id; // The current Octave id [0 -> Octaver::m_nb_octave]
};

} // namespace features
} // namespace openMVG

#endif // OPENMVG_FEATURES_SIFT_HIERARCHICAL_GAUSSIAN_SCALE_SPACE_HPP
