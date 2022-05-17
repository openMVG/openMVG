// This file is part of a pose-only algorithm of Linear Global Translation (LiGT)

// Copyright (c) 2022, Qi Cai and Yuanxin Wu

// This Source Code Form is subject to the license terms of
// Creative Commons Attribution Share Alike 4.0 International.
// Details are available at https://choosealicense.com/licenses/cc-by-sa-4.0/

#ifndef LIGT_ALGORITHM
#define LIGT_ALGORITHM

#pragma once
#include "LiGT_types.hpp"
#include <string>

namespace LiGT {

// ============== The LiGT Algorithm (Version 1.1) =============
// [Version History]
// v1.0: first release; parallelism by Pierre Moulon.
// v1.1: Spectra replaces Eigen; Block manipulation to implement LTL matrix.
//
// Coded by: Drs. Qi Cai and Xinrui Li
// Refined by: Pierre Moulon
// Email: sipangyiyou@sjtu.edu.cn, yuanxin.wu@sjtu.edu.cn
//
// [Conditions of Use]: the LiGT algorithm is distributed under
// the License of Attribution-ShareAlike 4.0 International
// (https://creativecommons.org/licenses/by-sa/4.0/).
//
//------------------
//-- Bibliography --
//------------------
// If you use it for a publication, please cite the following paper:
//- [1] "A Pose-only Solution to Visual Reconstruction and Navigation".
//- Authors: Qi Cai, Lilian Zhang, Yuanxin Wu, Wenxian Yu and Dewen Hu.
//- Date: December 2021.
//- Journal: IEEE T-PAMI.
//
//- [2] "Equivalent constraints for two-view geometry: pose solution/pure rotation identification and 3D reconstruction".
//- Authors: Qi Cai; Yuanxin Wu; Lilian Zhang; Peike Zhang.
//- Date: December 2019.
//- Journal: IJCV.
//
// This is a dedicated version of the LiGT algorithm for openMVG, which is supported by
// the Inertial and Visual Fusion (VINF) research group in Shanghai Jiao Tong University
// @ https://www.researchgate.net/lab/Inertial-visual-Fusion-VINF-Yuanxin-Wu
//
// Note:
// 1. The LiGT algorithm is subject to the University proprietary right with patent protection.
// The current version of the LiGT algorithm is generally consistent with
// the above-mentioned T-PAMI paper in terms of accuracy.
//
// 2. It does not consider the rank condition in Proposition 6 of the T-PAMI paper.
//

#include <Eigen/Core>

class LiGTProblem {
public:
  LiGTProblem();

  virtual ~LiGTProblem() = default;

  // [Step.2 in Pose-only Algorithm]: select the left/right-base views
  void SelectBaseViews(const Track& track,
             ViewId& lbase_view_id,
             ViewId& rbase_view_id,
             ObsId& id_lbase,
             ObsId& id_rbase);

  // [Step.3 in Pose-only algorithm]: calculate local L matrix, update LTL and A_lr matrix
  void BuildLTL(Eigen::MatrixXd& LTL,
          Eigen::MatrixXd& A_lr);

  //[Step.4 in Pose-only Algorithm]: obtain the translation solution by using SVD
  bool SolveLiGT(const Eigen::MatrixXd& LTL,
           Eigen::VectorXd &evectors);

  // [Step.5 in Pose-only Algorithm]: identify the correct sign of the translation solution after using SVD
  void IdentifySign(const Eigen::MatrixXd& A_lr,
            Eigen::VectorXd& evectors);

  // LiGT solution
  bool Solution();

  // get tracks pointer (help to set track inputs)
  Tracks GetTracks();

  // get rotations (help to set rotation inputs)
  Rotations GetRotations();

  // get poses
  Poses GetPoses();

  // check information in tracks (such as num_view/num_pts/num_obs)
  // and build estimated information EstInfo
  void CheckTracks();

  // recover the view id into original view ids
  void RecoverViewIds();

  // print copyright claim
  void PrintCopyright() const;

protected:
  unsigned int num_view_;
  unsigned int num_pts_;
  unsigned int num_obs_;
  // set the minimal number of image observations in a track
  // recommend value: 2~3 (default -> 2)
  unsigned int min_track_length_;
  double time_use_;

  // tracks
  LiGT::Tracks tracks_;

  // estimated view information
  LiGT::EstInfo est_info_;

  // [Note]: global rotations are sorted by estimated view ids
  LiGT::Rotations global_rotations_;

  // [Note]: global translations are sorted by estimated view ids
  LiGT::Translations global_translations_;

  // [Note]: global poses are sorted by original view ids
  LiGT::Poses poses_;

  // reference view id
  ViewId fixed_id_;
};

}
#endif
