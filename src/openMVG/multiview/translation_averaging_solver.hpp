//  Copyright (c) 2014, Kyle Wilson
//  All rights reserved.
//
//  Redistribution and use in source and binary forms, with or without
//  modification, are permitted provided that the following conditions are met:
//
//  1. Redistributions of source code must retain the above copyright notice, this
//  list of conditions and the following disclaimer.
//  2. Redistributions in binary form must reproduce the above copyright notice,
//  this list of conditions and the following disclaimer in the documentation
//  and/or other materials provided with the distribution.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
//  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
//  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
//  ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
//  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
//  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
//  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
//  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
//  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
//  The views and conclusions contained in the software and documentation are those
//  of the authors and should not be interpreted as representing official policies,
//  either expressed or implied, of the FreeBSD Project.

// Copyright (c) 2014 Pierre MOULON (Updated for openMVG).

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef __TRANS_SOLVER_H__
#define __TRANS_SOLVER_H__

//------------------
//-- Bibliography --
//------------------
//- [1] "Robust Global Translations with 1DSfM."
//- Authors: Kyle Wilson and Noah Snavely.
//- Date: September 2014.
//- Conference: ECCV.
//
//- [2] "Global Fusion of Relative Motions for Robust, Accurate and Scalable Structure from Motion."
//- Authors: Pierre Moulon, Pascal Monasse and Renaud Marlet.
//- Date: December 2013.
//- Conference: ICCV.

namespace openMVG {

/// Implementation of [1] : "5. Solving the Translations Problem" equation (3)
/// Compute camera center positions from relative camera translations (translation directions).
bool
solve_translations_problem_l2_chordal
(
  const int* edges,
  const double* poses,
  const double* weights,
  int num_edges,
  double loss_width,
  double* X,
  double function_tolerance,
  double parameter_tolerance,
  int max_iterations
);

/**
* @brief Registration of relative translations to global translations. It implements LInf minimization of  [2]
*  as a SoftL1 minimization. It can use 2/3 views relative translation vectors (bearing, or triplets of translations).
*
* @param[in] vec_initial_estimates relative motion information
* @param[in] b_translation_triplets tell if relative motion comes 3 or 2 views
*   false: 2-view estimates -> 1 relativeInfo per 2 view estimates,
*   true:  3-view estimates -> triplet of translations: 3 relativeInfo per triplet.
* @param[in] nb_poses the number of camera nodes in the relative motion graph
* @param[out] translations found global camera translations
* @param[in] d_l1_loss_threshold optionnal threshold for SoftL1 loss (-1: no loss function)
* @return True if the registration can be solved
*/
bool
solve_translations_problem_softl1
(
  const std::vector<openMVG::relativeInfo > & vec_initial_estimates,
  const bool b_translation_triplets,
  const int nb_poses,
  std::vector<Eigen::Vector3d> & translations,
  const double d_l1_loss_threshold = 0.01
);

} // namespace openMVG

#endif /* __TRANS_SOLVER_H__ */
