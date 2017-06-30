// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_FEATURES_SIFT_OCTAVER_HPP
#define OPENMVG_FEATURES_SIFT_OCTAVER_HPP

namespace openMVG { namespace image { template <typename T> class Image; } }

namespace openMVG{
namespace features{

template <typename OctaveT>
struct Octaver
{
  public:
  /**
  * @brief Octaver constructor
  * @param nb_octave Number of Octave
  * @param nb_slice Number of slice per octave
  */
  Octaver(const int nb_octave = 4, const int nb_slice = 3)
  :m_nb_octave(nb_octave), m_nb_slice(nb_slice)
  {};

  /**
  * @brief Set Initial image and update nb_octave if necessary
  * @param img Input image
  */
  virtual void SetImage(const image::Image<float> & img) = 0;

  /**
  * @brief Compute a full octave
  * @param[out] oct Computed octave
  * @retval true If an octave have been computed
  * @retval false if no octave have been computed (process ended)
  */
  virtual bool NextOctave(OctaveT & oct) = 0;

  /**
  * @brief Get number of slice in an octave
  */
  int NbSlice() const { return m_nb_slice;}

  /**
  * @brief Get number of octave computed
  * @note maybe less than specified in constructor because of input image size
  */
  int NbOctave() const {return m_nb_octave;}

protected:
  // Constant parameters
  int m_nb_octave;
  int m_nb_slice;
};

} // namespace features
} // namespace openMVG

#endif // OPENMVG_FEATURES_SIFT_OCTAVER_HPP
