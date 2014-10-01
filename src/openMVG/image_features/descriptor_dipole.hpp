// Copyright (c) 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_IMAGE_DESCRIPTION_DIPOLE_DESCRIPTOR_H
#define OPENMVG_IMAGE_DESCRIPTION_DIPOLE_DESCRIPTOR_H

#include "openMVG/features/feature.hpp"
#include "openMVG/numeric/math_trait.hpp"

//------------------
//-- Bibliography --
//------------------
//- [1] "New local descriptor based on dissociated dipoles"
//- Authors: Alexis Joly.
//- Date: December 2007.
//- Conference: CIVR.

namespace openMVG {
namespace image_features{

  // Creates a naive dissociated dipole describer [1].
  //
  // Note :
  // - Angle is in radian.
  // - data the output array (must be allocated to 20 values).
  template<typename Real, typename T>
  void PickNaiveDipole(
    const Image<Real> & image,
    float x, float y, float scale, float angle,
    T * data)
  {
    // Setup the rotation center.
    const float & cx = x, & cy = y;

    const float lambda1 = scale;
    const float lambda2 = lambda1 / 2.0f;
    static const float angleSubdiv = 2.0f * M_PI / 12.0f;
    static const float memoizeCos[12] =
    {
      MathTrait<float>::cos(angle + 0 * angleSubdiv),
      MathTrait<float>::cos(angle + 1 * angleSubdiv),
      MathTrait<float>::cos(angle + 2 * angleSubdiv),
      MathTrait<float>::cos(angle + 3 * angleSubdiv),
      MathTrait<float>::cos(angle + 4 * angleSubdiv),
      MathTrait<float>::cos(angle + 5 * angleSubdiv),
      MathTrait<float>::cos(angle + 6 * angleSubdiv),
      MathTrait<float>::cos(angle + 7 * angleSubdiv),
      MathTrait<float>::cos(angle + 8 * angleSubdiv),
      MathTrait<float>::cos(angle + 9 * angleSubdiv),
      MathTrait<float>::cos(angle + 10 * angleSubdiv),
      MathTrait<float>::cos(angle + 11 * angleSubdiv)
    };
    static const float memoizeSin[12] =
    {
      MathTrait<float>::sin(angle + 0 * angleSubdiv),
      MathTrait<float>::sin(angle + 1 * angleSubdiv),
      MathTrait<float>::sin(angle + 2 * angleSubdiv),
      MathTrait<float>::sin(angle + 3 * angleSubdiv),
      MathTrait<float>::sin(angle + 4 * angleSubdiv),
      MathTrait<float>::sin(angle + 5 * angleSubdiv),
      MathTrait<float>::sin(angle + 6 * angleSubdiv),
      MathTrait<float>::sin(angle + 7 * angleSubdiv),
      MathTrait<float>::sin(angle + 8 * angleSubdiv),
      MathTrait<float>::sin(angle + 9 * angleSubdiv),
      MathTrait<float>::sin(angle + 10 * angleSubdiv),
      MathTrait<float>::sin(angle + 11 * angleSubdiv)
    };

    Vecf dipoleF1(12);
    for (int i = 0; i < 12; ++i)  {
      const float xi = cx + lambda1 * memoizeCos[i];
      const float yi = cy + lambda1 * memoizeSin[i];
      // Bilinear interpolation
      dipoleF1(i) = SampleLinear(image, yi, xi);
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
    for (int i = 0; i < 12; ++i)  {
      const float xi = cx + (lambda1 + lambda2) * memoizeCos[i];
      const float yi = cy + (lambda1 + lambda2) * memoizeSin[i];

      const float xii = cx + (lambda1 - lambda2) * memoizeCos[i];
      const float yii = cy + (lambda1 - lambda2) * memoizeSin[i];

      // Bilinear interpolation
      dipoleF2(i) = SampleLinear(image, yi, xi) - SampleLinear(image, yii, xii);
    }
    // Normalize to be affine luminance invariant (a*I(x,y)+b).
    Map< Vecf > dataMap( data, 20);
    dataMap.block<8,1>(0,0) = (A * dipoleF1).normalized();
    dataMap.block<12,1>(8,0) = dipoleF2.normalized();
  }

  /// Pick an angular smoothed dissociated dipole
  template<typename Real, typename T>
  void PickASDipole(
    const Image<Real> & image,
    float x, float y, float scale, float angle,
    T * data)
  {
    // Setup the rotation center.
    const float & cx = x, & cy = y;

    const float lambda1 = scale;
    const float lambda2 = lambda1 / 2.0f;
    const float angleSubdiv = 2.0f * M_PI / 12.0f;

    //-- First order dipole :
    Vecf dipoleF1(12);
    for (int i = 0; i < 12; ++i)  {
      const float xi = cx + lambda1 * MathTrait<float>::cos(angle + i * angleSubdiv);
      const float yi = cy + lambda1 * MathTrait<float>::sin(angle + i * angleSubdiv);
      const float xi0 = cx + lambda1 * MathTrait<float>::cos(angle + i * angleSubdiv - angleSubdiv/2.0);
      const float yi0 = cy + lambda1 * MathTrait<float>::sin(angle + i * angleSubdiv - angleSubdiv/2.0);
      const float xi3 = cx + lambda1 * MathTrait<float>::cos(angle + i * angleSubdiv + angleSubdiv/2.0);
      const float yi3 = cy + lambda1 * MathTrait<float>::sin(angle + i * angleSubdiv + angleSubdiv/2.0);
      // Bilinear interpolation
      dipoleF1(i) =
        (SampleLinear(image, yi, xi) +
         SampleLinear(image, yi0, xi0) +
         SampleLinear(image, yi3, xi3))/3.0f;
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

    //-- Second order dipole :
    Vecf dipoleF2(12);
    for (int i = 0; i < 12; ++i)  {
      const double angleSample = i * angleSubdiv;
      const float xi = cx + (lambda1 + lambda2) * MathTrait<float>::cos(angle + angleSample);
      const float yi = cy + (lambda1 + lambda2) * MathTrait<float>::sin(angle + angleSample);
      const float xi0 = cx + (lambda1 + lambda2) * MathTrait<float>::cos(angle + angleSample - angleSubdiv/2.0);
      const float yi0 = cy + (lambda1 + lambda2) * MathTrait<float>::sin(angle + angleSample - angleSubdiv/2.0);
      const float xi3 = cx + (lambda1 + lambda2) * MathTrait<float>::cos(angle + angleSample + angleSubdiv/2.0);
      const float yi3 = cy + (lambda1 + lambda2) * MathTrait<float>::sin(angle + angleSample + angleSubdiv/2.0);

      const float xii = cx + (lambda1 - lambda2) * MathTrait<float>::cos(angle + angleSample);
      const float yii = cy + (lambda1 - lambda2) * MathTrait<float>::sin(angle + angleSample);
      const float xii0 = cx + (lambda1 - lambda2) * MathTrait<float>::cos(angle + angleSample - angleSubdiv/2.0);
      const float yii0 = cy + (lambda1 - lambda2) * MathTrait<float>::sin(angle + angleSample - angleSubdiv/2.0);
      const float xii3 = cx + (lambda1 - lambda2) * MathTrait<float>::cos(angle + angleSample + angleSubdiv/2.0);
      const float yii3 = cy + (lambda1 - lambda2) * MathTrait<float>::sin(angle + angleSample + angleSubdiv/2.0);

      // Bilinear interpolation
      dipoleF2(i) =
        (SampleLinear(image, yi, xi)   +
         SampleLinear(image, yi0, xi0) +
         SampleLinear(image, yi3, xi3) )/3.0f
        -
       (SampleLinear(image, yii, xii)   +
        SampleLinear(image, yii0, xii0) +
        SampleLinear(image, yii3, xii3)) /3.0f;
    }
    // Normalize to be affine luminance invariant (a*I(x,y)+b).
    Map< Vecf > dataMap( data, 20);
    dataMap.block<8,1>(0,0) = (A * dipoleF1).normalized();
    dataMap.block<12,1>(8,0) = dipoleF2.normalized();
  }

   /**
    ** @brief Compute a dissociated DIPOLE descriptor for a given interest point
    ** @param Li Input Octave slice
    ** @param Lx Input X-derivative
    ** @param Ly Input Y-derivative
    ** @param id_octave Id of current octave
    ** @param ipt Input interest point
    ** @param desc output descriptor (floating point descriptor)
    ** @param bAngularSmoothedDipole Extract or not a angular smoothed dipole
    **/
  template< typename Real>
  void ComputeDipoleDescriptor(
    const Image<Real> & Li , const Image<Real> &Lx , const Image<Real> &Ly ,
    const int id_octave , const SIOPointFeature & ipt , Descriptor<float, 20> & desc,
    bool bAngularSmoothedDipole = true)
  {
    const float factor_size = 3.5f;

    // Sampling size according to the scale value
    const float inv_octave_scale = static_cast<float>( 1.f ) /  static_cast<float>( 1 << id_octave );
    const float sigma_scale = MathTrait<float>::round( ipt.scale() * inv_octave_scale );

    const float x = ipt.x() * inv_octave_scale;
    const float y = ipt.y() * inv_octave_scale;

    if (bAngularSmoothedDipole)
      PickASDipole(Li, x, y, sigma_scale * factor_size, ipt.orientation(), &desc[0]);
    else
      PickNaiveDipole(Li, x, y, sigma_scale * factor_size, ipt.orientation(), &desc[0]);
  }
} // namespace image_features
} // namespace openMVG

#endif
