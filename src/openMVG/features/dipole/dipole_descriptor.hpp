// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_FEATURES_DIPOLE_DIPOLE_DESCRIPTOR_HPP
#define OPENMVG_FEATURES_DIPOLE_DIPOLE_DESCRIPTOR_HPP


#include "openMVG/features/feature.hpp"
#include "openMVG/features/descriptor.hpp"
#include "openMVG/image/image_container.hpp"
#include "openMVG/image/sample.hpp"

//------------------
//-- Bibliography --
//------------------
//- [1] "New local descriptor based on dissociated dipoles"
//- Authors: Alexis Joly.
//- Date: December 2007.
//- Conference: CIVR.

namespace openMVG
{
namespace features
{
  // Create a DIPOLE describer [1].
  //
  // Note :
  // - Angle is in radians.
  // - data the output array (must be allocated to 20 values).
  template<typename Real>
  void PickNaiveDipole
  (
    const image::Image<Real> & image,
    float x,
    float y,
    float scale,
    float angle,
    float * data
  )
  {
    // Use bilinear sampling
    const image::Sampler2d<image::SamplerLinear> sampler;
    // Setup the rotation center.
    const float & cx = x, & cy = y;

    const float lambda1 = scale;
    const float lambda2 = lambda1 / 2.0f;
    static const float angleSubdiv = 2.0f * M_PI / 12.0f;
    const float memoizeCos[12] =
    {
      std::cos(angle + 0.f * angleSubdiv),
      std::cos(angle + 1.f * angleSubdiv),
      std::cos(angle + 2.f * angleSubdiv),
      std::cos(angle + 3.f * angleSubdiv),
      std::cos(angle + 4.f * angleSubdiv),
      std::cos(angle + 5.f * angleSubdiv),
      std::cos(angle + 6.f * angleSubdiv),
      std::cos(angle + 7.f * angleSubdiv),
      std::cos(angle + 8.f * angleSubdiv),
      std::cos(angle + 9.f * angleSubdiv),
      std::cos(angle + 10.f * angleSubdiv),
      std::cos(angle + 11.f * angleSubdiv)
    };
    const float memoizeSin[12] =
    {
      std::sin(angle + 0.f * angleSubdiv),
      std::sin(angle + 1.f * angleSubdiv),
      std::sin(angle + 2.f * angleSubdiv),
      std::sin(angle + 3.f * angleSubdiv),
      std::sin(angle + 4.f * angleSubdiv),
      std::sin(angle + 5.f * angleSubdiv),
      std::sin(angle + 6.f * angleSubdiv),
      std::sin(angle + 7.f * angleSubdiv),
      std::sin(angle + 8.f * angleSubdiv),
      std::sin(angle + 9.f * angleSubdiv),
      std::sin(angle + 10.f * angleSubdiv),
      std::sin(angle + 11.f * angleSubdiv)
    };

    Vecf dipoleF1(12);
    for (int i = 0; i < 12; ++i)
    {
      const float xi = cx + lambda1 * memoizeCos[i];
      const float yi = cy + lambda1 * memoizeSin[i];
      dipoleF1(i) = sampler(image, yi, xi);
    }
    Matf A(8,12);
    A <<  0, 0, 0, 1, 0, 0, 0, 0, 0,-1, 0, 0,
          0,-1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
          0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0, 1,
          0, 0, 0, 0, 1, 0, 0,-1, 0, 0, 0, 0,
          0, 0, 0, 0, 0, 0, 1, 0, 0,-1, 0, 0,
          0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,-1,
          0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
          1, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0;

    // Add the second order F2 dipole
    Vecf dipoleF2(12);
    for (int i = 0; i < 12; ++i)
    {
      const float xi = cx + (lambda1 + lambda2) * memoizeCos[i];
      const float yi = cy + (lambda1 + lambda2) * memoizeSin[i];

      const float xii = cx + (lambda1 - lambda2) * memoizeCos[i];
      const float yii = cy + (lambda1 - lambda2) * memoizeSin[i];

      // Bilinear interpolation
      dipoleF2(i) = sampler(image, yi, xi) - sampler(image, yii, xii);
    }
    // Normalize to be affine luminance invariant (a*I(x,y)+b).
    Map<Vecf> dataMap( data, 20);
    dataMap.block<8,1>(0,0) = (A * dipoleF1).normalized();
    dataMap.block<12,1>(8,0) = dipoleF2.normalized();
  }

  // Pick an angular smoothed dipole
  template<typename Real>
  void PickASDipole
  (
    const image::Image<Real> & image,
    float x,
    float y,
    float scale,
    float angle,
    float * data)
  {
    const image::Sampler2d<image::SamplerLinear> sampler;
    // Setup the rotation center.
    const float & cx = x, & cy = y;

    const float lambda1 = scale;
    const float lambda2 = lambda1 / 2.0f;
    const float angleSubdiv = 2.0f * M_PI / 12.0f;

    //-- First order dipole:
    Vecf dipoleF1(12);
    for (int i = 0; i < 12; ++i)
    {
      const float xi = cx + lambda1 * std::cos(angle + i * angleSubdiv);
      const float yi = cy + lambda1 * std::sin(angle + i * angleSubdiv);
      const float xi0 = cx + lambda1 * std::cos(angle + i * angleSubdiv - angleSubdiv/2.0);
      const float yi0 = cy + lambda1 * std::sin(angle + i * angleSubdiv - angleSubdiv/2.0);
      const float xi3 = cx + lambda1 * std::cos(angle + i * angleSubdiv + angleSubdiv/2.0);
      const float yi3 = cy + lambda1 * std::sin(angle + i * angleSubdiv + angleSubdiv/2.0);
      // Bilinear interpolation
      dipoleF1(i) =
        (sampler(image, yi, xi) +
         sampler(image, yi0, xi0) +
         sampler(image, yi3, xi3))/3.0f;
    }
    Matf A(8,12);
    A <<  0, 0, 0, 1, 0, 0, 0, 0, 0,-1, 0, 0,
          0,-1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
          0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0, 1,
          0, 0, 0, 0, 1, 0, 0,-1, 0, 0, 0, 0,
          0, 0, 0, 0, 0, 0, 1, 0, 0,-1, 0, 0,
          0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,-1,
          0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
          1, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0;

    //-- Second order dipole:
    Vecf dipoleF2(12);
    for (int i = 0; i < 12; ++i)
    {
      const double angleSample = i * angleSubdiv;
      const float xi = cx + (lambda1 + lambda2) * std::cos(angle + angleSample);
      const float yi = cy + (lambda1 + lambda2) * std::sin(angle + angleSample);
      const float xi0 = cx + (lambda1 + lambda2) * std::cos(angle + angleSample - angleSubdiv/2.0);
      const float yi0 = cy + (lambda1 + lambda2) * std::sin(angle + angleSample - angleSubdiv/2.0);
      const float xi3 = cx + (lambda1 + lambda2) * std::cos(angle + angleSample + angleSubdiv/2.0);
      const float yi3 = cy + (lambda1 + lambda2) * std::sin(angle + angleSample + angleSubdiv/2.0);

      const float xii = cx + (lambda1 - lambda2) * std::cos(angle + angleSample);
      const float yii = cy + (lambda1 - lambda2) * std::sin(angle + angleSample);
      const float xii0 = cx + (lambda1 - lambda2) * std::cos(angle + angleSample - angleSubdiv/2.0);
      const float yii0 = cy + (lambda1 - lambda2) * std::sin(angle + angleSample - angleSubdiv/2.0);
      const float xii3 = cx + (lambda1 - lambda2) * std::cos(angle + angleSample + angleSubdiv/2.0);
      const float yii3 = cy + (lambda1 - lambda2) * std::sin(angle + angleSample + angleSubdiv/2.0);

      // Bilinear interpolation
      dipoleF2(i) =
        (sampler(image, yi, xi)   +
         sampler(image, yi0, xi0) +
         sampler(image, yi3, xi3)) /3.0f
        -
       (sampler(image, yii, xii)   +
        sampler(image, yii0, xii0) +
        sampler(image, yii3, xii3)) /3.0f;
    }
    // Normalize to be affine luminance invariant (a*I(x,y)+b).
    Map<Vecf> dataMap( data, 20);
    dataMap.block<8,1>(0,0) = (A * dipoleF1).normalized();
    dataMap.block<12,1>(8,0) = dipoleF2.normalized();
  }

   /**
    ** @brief Compute DIPOLE descriptor for a given interest point
    ** @param Li Input image
    ** @param ipt Input interest point
    ** @param desc output descriptor (floating point descriptor)
    ** @param bAngularSmoothedDipole Tell if we must extract an upright or an
    **  angular smoothed dipole
    ** @param magnif_factor Scaling factor used to rescale the dipole sampling
    **/
  template<typename Real>
  void ComputeDipoleDescriptor
  (
    const image::Image<Real> & Li,
    const SIOPointFeature & ipt,
    Descriptor<float, 20> & desc,
    bool bAngularSmoothedDipole = true,
    const float magnif_factor = 3.5f
  )
  {
    if (bAngularSmoothedDipole)
      PickASDipole(Li, ipt.x(), ipt.y(), ipt.scale() * magnif_factor, ipt.orientation(), &desc[0]);
    else
      PickNaiveDipole(Li, ipt.x(), ipt.y(), ipt.scale() * magnif_factor, ipt.orientation(), &desc[0]);
  }

} // namespace features
} // namespace openMVG

#endif // OPENMVG_FEATURES_DIPOLE_DIPOLE_DESCRIPTOR_HPP
