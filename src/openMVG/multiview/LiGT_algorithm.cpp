// This file is part of a pose-only algorithm of Linear Global Translation (LiGT)

// Copyright (c) 2022, Qi Cai and Yuanxin Wu

// This Source Code Form is subject to the license terms of
// Creative Commons Attribution Share Alike 4.0 International.
// Details are available at https://choosealicense.com/licenses/cc-by-sa-4.0/

#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SVD>
#include <chrono>
#include <fstream>
#include <iomanip>

#include "LiGT_algorithm.hpp"

// set the minimal number of image observations in a track
#define MIN_TRACKING_LENGTH 1000000

using namespace std;
using namespace chrono;
using namespace LiGT;

namespace LiGT {

// ====================   LiGT Problem ========================

LiGTProblem::LiGTProblem() {
    num_view_ = 0;
    num_pts_ = 0;
    num_obs_ = 0;
    time_use_ = 0;
    fixed_id_ = 0;

    output_file_ = "";
    time_file_ = "";

}


LiGTProblem::LiGTProblem(const std::string& globalR_file,
                         const std::string& track_file,
                         const std::string& output_file,
                         const std::string& time_file,
                         const int& fixed_id) {
    num_view_ = 0;
    num_pts_ = 0;
    num_obs_ = 0;
    time_use_ = 0;
    fixed_id_ = fixed_id;

    output_file_ = output_file;
    time_file_ = time_file;
    cout<<"====================================================="<<endl;
    cout<<"##  input info in LiGT algorithm (basic version)  ##"<<endl;
    cout<<"Rotation file:="<<globalR_file<<endl;
    cout<<"Track file:="<<track_file<<endl;
    cout<<"fixed id:="<<fixed_id<<endl;
    cout<<"-----------------------------------------------------"<<endl;
    SetupOrientation(globalR_file);
    LoadTracks(track_file);
}

LiGTProblem::LiGTProblem(LiGTProblem& problem){
    tracks_ = problem.tracks_;

    time_use_ = problem.time_use_;
    fixed_id_ = problem.fixed_id_;
    num_view_ = problem.num_view_;
    num_pts_ = problem.num_pts_;
    num_obs_ = problem.num_obs_;
    output_file_ = problem.output_file_;
    time_file_ = problem.time_file_;
    est_info_.origin_view_ids = problem.est_info_.origin_view_ids;
    est_info_.estimated_view_ids = problem.est_info_.estimated_view_ids;
    est_info_.est2origin_view_ids = problem.est_info_.est2origin_view_ids;
    est_info_.origin2est_view_ids = problem.est_info_.origin2est_view_ids;

    global_rotations_ = problem.global_rotations_;
    global_translations_ = problem.global_translations_;
    poses_ = problem.poses_;

}
LiGTProblem::LiGTProblem(Tracks tracks,
                         Attitudes global_R){

    tracks_.swap(tracks);
    global_rotations_.swap(global_R);

    // check tracks and build estimated information [EstInfo]
    if (!tracks_.size()){
        cout<<"error: wrong track inputs to run LiGT" << endl;
        return;
    }

    if (!global_rotations_.size()){
        cout<<"error: wrong rotation inputs to run LiGT" << endl;
        return;
    }

    time_use_ = 0;
    fixed_id_ = 0;

    CheckTracks();
}

void LiGTProblem::Init(const int& fixed_id){

    // check tracks and build estimated information [EstInfo]
    if (!tracks_.size()){
        cout<<"error: wrong track inputs to run LiGT" << endl;
        return;
    }

    if (!global_rotations_.size()){
        cout<<"error: wrong rotation inputs to run LiGT" << endl;
        return;
    }

    time_use_ = 0;
    fixed_id_ = fixed_id;

    CheckTracks();
}

void LiGTProblem::SetupOrientation(const std::string& globalR_file) {

    auto start_time = steady_clock::now();

    // load global rotation file
    fstream infile;
    infile.open(globalR_file, std::ios::in);
    if (!infile.is_open()) {
        cout << "gR_file cannot open" << endl;
        return;
    }

    infile >> num_view_;

    global_rotations_.clear();
    global_rotations_.resize(num_view_);

    for (ViewId i = 0; i < num_view_; i++) {
        Eigen::Matrix3d tmp_R;
        infile >> tmp_R(0, 0) >> tmp_R(1, 0) >> tmp_R(2, 0);
        infile >> tmp_R(0, 1) >> tmp_R(1, 1) >> tmp_R(2, 1);
        infile >> tmp_R(0, 2) >> tmp_R(1, 2) >> tmp_R(2, 2);

        global_rotations_[i] = (tmp_R);
    }

    infile.close();

    // print time cost information
    auto end_time = steady_clock::now();
    auto duration = duration_cast<microseconds>(end_time - start_time);
    cout << ">> time for loading global rotation file: "
         << double(duration.count()) * microseconds::period::num / microseconds::period::den
         << "s" << endl;

}

void LiGTProblem::LoadTracks(const std::string& track_file) {

    auto start_time = steady_clock::now();

    // load tracks file
    std::fstream tracks_file(track_file, ios::in);
    if (!tracks_file.is_open()) {
        std::cout << "tracks file cannot load" << std::endl;
        return;
    }

    // load header info
    tracks_file >> num_view_ >> num_pts_ >> num_obs_;

    // load tracks
    int record_pts_id = -1;

    TrackInfo track_info;
    ObsInfo tmp_obs;

    tracks_.reserve(num_pts_);

    for (ObsId i = 0; i < num_obs_; i++) {

        tracks_file >> tmp_obs.view_id
                >> tmp_obs.pts_id
                >> tmp_obs.coord(0)
                >> tmp_obs.coord(1);

        tmp_obs.coord(0) = -tmp_obs.coord(0);
        tmp_obs.coord(1) = -tmp_obs.coord(1);
        tmp_obs.coord(2) = 1;

        if (tmp_obs.pts_id != record_pts_id) {
            if (i > 0) {
                if (track_info.track.size() < MIN_TRACKING_LENGTH) {
                    //only need use observation of tracking length >= MIN_TRACKING_LENGTH
                    track_info.is_used = 0;
                } else {
                    track_info.is_used = 1;
                }
                tracks_.emplace_back(track_info);
            }
            track_info.track.clear();
        }
        track_info.track.emplace_back(tmp_obs);
        record_pts_id = tmp_obs.pts_id;
    }
    tracks_file.close();
    tracks_.emplace_back(track_info);

    // check tracks and build estimated information [EstInfo]
    CheckTracks();

    // print time cost information
    auto end_time = steady_clock::now();
    auto duration = duration_cast<microseconds>(end_time - start_time);
    cout << ">> time for loading tracks file: "
         << double(duration.count()) * microseconds::period::num / microseconds::period::den
         << "s" << endl;

}

Matrix3d LiGTProblem::hat(const Vector3d& vector) {

    Matrix3d hat_vector;
    hat_vector <<  0, -vector(2),  vector(1)
            ,  vector(2),     0, -vector(0)
            , -vector(1),  vector(0),     0;
    return hat_vector;

}


void LiGTProblem::CheckTracks(){

    cout << "checking tracks information..." << endl;

    Size tmp_num_obs = 0;

    std::set<ViewId> est_view;
    for ( PtsId i = 0; i < tracks_.size(); ++i){
        Track& track = tracks_[i].track;
        for (ObsId j = 0; j < track.size(); ++j){
            est_view.insert(track[j].view_id);
        }
        tmp_num_obs += track.size();
    }

    // check information
    if (num_view_ != est_view.size())
    {
        cout<<"warning: the number of input views might be wrong!"<<endl;
        num_view_ = est_view.size();
    }

    if (num_pts_ != tracks_.size())
    {
        cout<<"warning: the number of input pts might be wrong!"<<endl;
        num_pts_ = tracks_.size();
    }

    if (num_obs_ != tmp_num_obs)
    {
        cout<<"warning: the number of input obs might be wrong!"<<endl;
        num_obs_ = tmp_num_obs;
    }

    // build estimated view information
    ViewId id = 0;
    for ( auto& view_id : est_view){
        est_info_.origin_view_ids.emplace_back(view_id);
        est_info_.estimated_view_ids.emplace_back(id);
        id ++;
    }

    // select the reference view id by the fixed camera id
    auto& origin_view_id = est_info_.origin_view_ids[0];
    for ( ViewId i = 0; i < est_info_.origin_view_ids.size(); ++i){
        if ( est_info_.origin_view_ids[i] == fixed_id_){
            est_info_.origin_view_ids[0] = est_info_.origin_view_ids[i];
            est_info_.origin_view_ids[i] = origin_view_id;
        }
    }

    est_info_.BuildMap();

    // update view id in track
    for ( PtsId i = 0; i < tracks_.size(); ++i){
        Track& track = tracks_[i].track;
        for ( ObsId j = 0; j < track.size(); ++j){
            track[j].view_id = est_info_.origin2est_view_ids[track[j].view_id];
        }
    }
}

void LiGTProblem::RecoverViewIds(){

    cout << "recover the estimated view ids into original view ids" << endl;

    for (ViewId i = 0; i < num_view_; ++i){
        Pose tmp_pose(global_rotations_[i],global_translations_[i]);
        poses_.insert({est_info_.est2origin_view_ids[i], tmp_pose});
    }

}

void LiGTProblem::WriteTranslation(const std::string output_file) {

    fstream output(output_file, std::ios::out | ios::trunc);
    if (!output.is_open()) {
        cout << "output file cannot create, please check path" << endl;
        return;
    }

    Size t_size = num_view_;
    Vector3d* tmp_position = nullptr;
    Matrix3d* tmp_R = nullptr;
    unsigned int i = 0;

    output << t_size << endl;

    for (i = 0; i < t_size; ++i) {
        tmp_position = &global_translations_[i];
        tmp_R = &global_rotations_[i];

        output << setprecision(16) << endl << (*tmp_R)(0, 0) << " " << (*tmp_R)(1, 0) << " " << (*tmp_R)(2, 0) << " ";
        output << setprecision(16) << (*tmp_R)(0, 1) << " " << (*tmp_R)(1, 1) << " " << (*tmp_R)(2, 1) << " ";
        output << setprecision(16) << (*tmp_R)(0, 2) << " " << (*tmp_R)(1, 2) << " " << (*tmp_R)(2, 2) << " ";

        output << setprecision(16) << (*tmp_position)(0) << " " << (*tmp_position)(1) << " " << (*tmp_position)(2);
    }
    output.close();

}

void LiGTProblem::WriteTime(const std::string& time_file) {

    fstream output(time_file, std::ios::out | ios::trunc);
    if (!output.is_open()) {
        cout << "time file cannot create, please check path" << endl;
        return;
    }
    output << setprecision(16) << time_use_ << ' ' << 1 << endl;
    output.close();

}

void LiGTProblem::UpdateLiGTMatrix(const ViewId& lbase_view_id,
                                     const ViewId& rbase_view_id,
                                     const ObsId& id_lbase,
                                     const ObsId& id_rbase,
                                     const Track& track,
                                     const PtsId& pts_id,
                                     Eigen::MatrixXd& A_lr,
                                     Eigen::MatrixXd& LTL){

    for (ObsId i = 0; i < track.size(); i++) {
        ViewId i_view_id = track[i].view_id; // the current view id

        if (i_view_id != lbase_view_id) {
            Eigen::Matrix3d xi_cross = hat(track[i].coord);
            Eigen::Matrix3d R_li = global_rotations_[i_view_id] * global_rotations_[lbase_view_id].transpose();
            Eigen::Matrix3d R_lr = global_rotations_[rbase_view_id] * global_rotations_[lbase_view_id].transpose();

            Eigen::Vector3d tmp_a_lr = hat(R_lr * track[id_lbase].coord)
                    * track[id_rbase].coord;

            // a_lr (Row) vector in a_lr * t > 0
            Eigen::RowVector3d a_lr = tmp_a_lr.transpose() * hat(track[id_rbase].coord);

            // combine all a_lr vectors into a matrix form A, i.e., At > 0
            A_lr.row(pts_id).block<1, 3>(0, lbase_view_id * 3) = a_lr * global_rotations_[rbase_view_id];
            A_lr.row(pts_id).block<1, 3>(0, rbase_view_id * 3) = -a_lr * global_rotations_[rbase_view_id];

            // theta_lr
            Eigen::Vector3d theta_lr_vector = hat(track[id_rbase].coord)
                    * R_lr
                    * track[id_lbase].coord;
            double theta_lr = theta_lr_vector.squaredNorm();

            // caculate matrix B
            Eigen::Matrix3d Coefficient_B =
                    xi_cross * R_li * track[id_lbase].coord * a_lr * global_rotations_[rbase_view_id];

            // caculate matrix C
            Eigen::Matrix3d Coefficient_C =
                    theta_lr * hat(track[i].coord) * global_rotations_[i_view_id];

            // caculate matrix D
            Eigen::Matrix3d Coefficient_D = -(Coefficient_B + Coefficient_C);

            // calculate temp matrix L for a single 3D matrix
            Eigen::MatrixXd tmp_LiGT_vec = Eigen::MatrixXd::Zero(3, num_view_ * 3);// L matrix for a single 3D point
            tmp_LiGT_vec.block<3, 3>(0, rbase_view_id * 3) += Coefficient_B;
            tmp_LiGT_vec.block<3, 3>(0, i_view_id * 3) += Coefficient_C;
            tmp_LiGT_vec.block<3, 3>(0, lbase_view_id * 3) += Coefficient_D;

            // delete the reference view column
            Eigen::MatrixXd tmp_vec = tmp_LiGT_vec.rightCols(tmp_LiGT_vec.cols() - 3);

            // calculate matrix LTL for SVD
            LTL += tmp_vec.transpose() * tmp_vec;
        }
    }

}

void LiGTProblem::IdentifySign(const MatrixXd& A_lr,
                                 VectorXd& evectors) {

    int positive_count = 0;
    int nagative_count = 0;

    VectorXd judgeValue = A_lr * evectors;

    for (int i = 0; i < judgeValue.size(); ++i) {
        if (judgeValue[i] > 0){
            positive_count++;
        }
        if (judgeValue[i] < 0){
            nagative_count++;
        }
    }

    if (positive_count < nagative_count) {
        evectors = -evectors;
    }

}


void LiGTProblem::Solution() {
    cout<<"************  LiGT Solve Summary  **************"<<endl;
    cout << "num_view = " << num_view_ << "; num_pts = " << num_pts_ << "; num_obs = " << num_obs_ << endl;

    // start time clock
    auto start_time = steady_clock::now();

    // allocate memory for LTL matrix where Lt=0
    Eigen::MatrixXd LTL = Eigen::MatrixXd::Zero(num_view_ * 3-3, num_view_ * 3-3);

    // use A_lr * t > 0 to identify the correct sign of the translation result
    Eigen::MatrixXd A_lr = Eigen::MatrixXd::Zero(tracks_.size(), 3 * num_view_);

    // construct LTL and A_lr matrix from 3D points
    for (PtsId track_id = 0; track_id < tracks_.size(); track_id++) {

        TrackInfo& track_info = tracks_[track_id];
        Track& track = track_info.track;

        double best_criterion_value = 0;

        ViewId lbase_view_id = 0;
        ViewId rbase_view_id = 0;

        ObsId id_lbase = 0;
        ObsId id_rbase = 0;

        Size track_size = track.size();

        // [Step.2 in Pose-only Algorithm]: select the left/right-base views
        for (ObsId i = 0; i < track_size - 1; i++) {
            for (ObsId j = i + 1; j < track_size; j++) {

                ViewId& i_view_id = track[i].view_id;
                ViewId& j_view_id = track[j].view_id;

                Vector3d& i_coord = track[i].coord;
                Vector3d& j_coord = track[j].coord;

                Matrix3d& R_i = global_rotations_[i_view_id];
                Matrix3d& R_j = global_rotations_[j_view_id];

                Matrix3d R_ij = R_j * R_i.transpose();
                Vector3d theta_ij = j_coord.cross(R_ij * i_coord);

                double criterion_value = theta_ij.norm();

                if (criterion_value > best_criterion_value) {

                    best_criterion_value = criterion_value;

                    lbase_view_id = track[i].view_id;
                    rbase_view_id = track[j].view_id;

                    id_lbase = i;
                    id_rbase = j;
                }
            }
        }

        // [Step.3 in Pose-only algorithm]: calculate local L matrix, update LTL and A_lr matrix
        UpdateLiGTMatrix( lbase_view_id, rbase_view_id,
                            id_lbase, id_rbase,
                            track, track_id,
                            A_lr, LTL);

    }

    //[Step.4 in Pose-only Algorithm]: obtain the translation solution by using SVD
    JacobiSVD<Eigen::MatrixXd> svd(LTL, ComputeFullU | ComputeFullV);
    MatrixXd V = svd.matrixV();
    VectorXd evectors = VectorXd::Zero(V.rows() + 3);
    evectors.bottomRows(V.rows()) = V.col(V.cols() - 1);

    //[Step.5 in Pose-only Algorithm]: identify the right global translation solution
    IdentifySign(A_lr, evectors);

    // algorithm time cost
    auto end_time = steady_clock::now();
    auto duration = duration_cast<microseconds>(end_time - start_time);
    cout <<  ">> time for LiGT algorithm:"
          << double(duration.count()) * microseconds::period::num / microseconds::period::den
          << "s" << endl;
    time_use_ = double(duration.count()) * microseconds::period::num / microseconds::period::den;

    // (optional) transform solution evectors into global_translations_
    global_translations_.clear();
    global_translations_.resize(num_view_);

    for (ViewId i = 0; i < num_view_; ++i) {
        global_translations_[i] = evectors.middleRows<3>(3*i);
    }

    // (optional) translation recovery
    RecoverViewIds();
}

Tracks LiGTProblem::GetTracks(){

    return tracks_;

}

Poses LiGTProblem::GetPoses(){

    return poses_;

}

Attitudes LiGTProblem::GetRotations(){

    return global_rotations_;
}


}
