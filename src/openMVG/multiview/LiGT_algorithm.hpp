// This file is part of a pose-only algorithm of Linear Global Translation (LiGT)

// Copyright (c) 2022, Qi Cai and Yuanxin Wu

// This Source Code Form is subject to the license terms of
// Creative Commons Attribution Share Alike 4.0 International.
// Details are available at https://choosealicense.com/licenses/cc-by-sa-4.0/

#ifndef LIGT_ALGORITHM
#define LIGT_ALGORITHM

#pragma once
#include "LiGT_types.hpp"
#include <iostream>

namespace LiGT {

// ============== The LiGT Algorithm (Version 1.0) =============
// Coded by: Drs. Qi Cai and Xinrui Li
// Email: sipangyiyou@sjtu.edu.cn, yuanxin.wu@sjtu.edu.cn
//
// Conditions of Use: the LiGT algorithm is distributed under
// the License of Attribution-ShareAlike 4.0 International
// (https://creativecommons.org/licenses/by-sa/4.0/).
//
// If you use it for a publication, please cite the following paper:
//- [1] "A Pose-only Solution to Visual Reconstruction and Navigation".
//- Authors: Qi Cai, Lilian Zhang, Yuanxin Wu, Wenxian Yu and Dewen Hu.
//- Date: December 2022.
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
// The Spectra library was not utilized in solving the SVD problem.


class LiGTProblem {
public:

    LiGTProblem();

    explicit LiGTProblem(LiGTProblem& problem);

    LiGTProblem(Tracks tracks,
                         Attitudes global_R);

    LiGTProblem(const std::string& globalR_file,
                         const std::string& track_file,
                         const std::string& output_file,
                         const std::string& time_file,
                         const int& fixed_id);

    virtual ~LiGTProblem() = default;

    // initialize the LiGT problem
    void Init(const int& fixed_id = 0);

    // set global rotation input
    void SetupOrientation(const std::string& globalR_file);

    // load tracks from a bal-format track file
    // [web details in bal format: https://grail.cs.washington.edu/projects/bal/]
    void LoadTracks(const std::string& track_file);

    // obtain the antisymmetric matrix from a vector
    Matrix3d hat(const Vector3d& v);

    // write translation result
    void WriteTranslation(const std::string output_file);

    // write time result
    void WriteTime(const std::string& time_file);

    // [Step.3 in Pose-only algorithm]: calculate local L matrix, update LTL and A_lr matrix
    // the update way of A_lr is based on the point id,
    // where the dimension of A_lr is [num_pts] x [3*num_view]
    void UpdateLiGTMatrix(const ViewId &lbase_view_id,
                            const ViewId &rbase_view_id,
                            const ObsId &id_lbase,
                            const ObsId &id_rbase,
                            const Track &track,
                            const PtsId &pts_id,
                            MatrixXd& A_lr,
                            Eigen::MatrixXd& LTL);

    // [Step.5 in Pose-only Algorithm]: identify the correct sign of the translation solution after using SVD
    void IdentifySign(const MatrixXd &A_lr,
                        VectorXd& evectors);

    // LiGT solution
    void Solution();

    // get tracks pointer (help to set track inputs)
    Tracks GetTracks();

    // get rotations (help to set rotation inputs)
    Attitudes GetRotations();

    // get poses
    Poses GetPoses();

    // check information in tracks (such as num_view/num_pts/num_obs)
    // and build estimated information EstInfo
    void CheckTracks();

    // recover the view id into original view ids
    void RecoverViewIds();

protected:

    unsigned int num_view_;
    unsigned int num_pts_;
    unsigned int num_obs_;
    double time_use_;

    // output files
    std::string output_file_;
    std::string time_file_;

    // tracks
    LiGT::Tracks tracks_;

    // estimated view information
    LiGT::EstInfo est_info_;

    // [Note]: global rotations are sorted by estimated view ids
    LiGT::Attitudes global_rotations_;

    // [Note]: global translations are sorted by estimated view ids
    LiGT::Translations global_translations_;

    // [Note]: global poses are sorted by original view ids
    LiGT::Poses poses_;

    // reference view id
    ViewId fixed_id_;
};


}
#endif
