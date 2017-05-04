// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_GEOMETRY_RIGID_TRANSFORMATION_3D_SRT_HPP
#define OPENMVG_GEOMETRY_RIGID_TRANSFORMATION_3D_SRT_HPP

#include "openMVG/numeric/lm.hpp"
#include "openMVG/numeric/eigen_alias_definition.hpp"

namespace openMVG
{
namespace geometry
{

/** 3D rigid transformation estimation (7 dof)
 * Compute a Scale Rotation and Translation rigid transformation.
 * This transformation provide a distortion-free transformation
 * using the following formulation Xb = S * R * Xa + t.
 * "Least-squares estimation of transformation parameters between two point patterns",
 * Shinji Umeyama, PAMI 1991, DOI: 10.1109/34.88573
 *
 * \param[in] x1 The first 3xN matrix of euclidean points
 * \param[in] x2 The second 3xN matrix of euclidean points
 * \param[out] S The scale factor
 * \param[out] t The 3x1 translation
 * \param[out] R The 3x3 rotation
 *
 * \return true if the transformation estimation has succeeded
 *
 * \note Need at least 3 points
 */
  bool FindRTS
  (
    const Mat &x1,
    const Mat &x2,
    double * S,
    Vec3 * t,
    Mat3 * R
  );

/**
* @brief Eigen Levemberg-Marquardt functor to refine translation, Rotation and Scale parameter.
*/
struct lm_SRTRefine_functor : Functor<double>
{
  /**
  * @brief Constructor
  * @param inputs Number of inputs (nb elements to refine)
  * @param values Number of samples
  * @param x1 Input samples for first dataset
  * @param x2 Input samples for second dataset
  * @param S Scale
  * @param R Rotation
  * @param t Translation
  */
  lm_SRTRefine_functor( int inputs, int values,
                        const Mat &x1, const Mat &x2,
                        const double &S, const Mat3 & R, const Vec &t );

  /**
  * @brief Computes error given a sample
  * @param x a Sample
  * @param[out] fvec Error for each values
  */
  int operator()( const Vec &x, Vec &fvec ) const;

  Mat x1_, x2_;
  Vec3 t_;
  Mat3 R_;
  double S_;
};


/**
* @brief Eigen LM functor to refine Rotation.
*/
struct lm_RRefine_functor : Functor<double>
{
  /**
  * @brief Constructor
  * @param inputs Number of inputs (elements to refine)
  * @param values Number of samples
  * @param x1 Input samples for first dataset
  * @param x2 Input samples for second dataset
  * @param S Scale
  * @param R Rotation
  * @param t Translation
  */
  lm_RRefine_functor( int inputs, int values,
                      const Mat &x1, const Mat &x2,
                      const double &S, const Mat3 & R, const Vec &t );

  /**
   * @brief Computes error given a sample
   * @param x a Sample
   * @param[out] fvec Error for each values
   */
  int operator()( const Vec &x, Vec &fvec ) const;

  Mat x1_, x2_;
  Vec3 t_;
  Mat3 R_;
  double S_;
};

/** 3D rigid transformation refinement using LM
 * Refine the Scale, Rotation and translation rigid transformation
 * using a Levenberg-Marquardt opimization.
 *
 * \param[in] x1 The first 3xN matrix of euclidean points
 * \param[in] x2 The second 3xN matrix of euclidean points
 * \param[out] S The initial scale factor
 * \param[out] t The initial 3x1 translation
 * \param[out] R The initial 3x3 rotation
 *
 * \return none
 */
void Refine_RTS
(
  const Mat &x1,
  const Mat &x2,
  double * S,
  Vec3 * t,
  Mat3 * R
);

} // namespace geometry
} // namespace openMVG

#endif  // OPENMVG_GEOMETRY_RIGID_TRANSFORMATION_3D_SRT_HPP
