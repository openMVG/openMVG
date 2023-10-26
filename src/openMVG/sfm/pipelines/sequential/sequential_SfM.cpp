// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON, Ricardo Fabbri and Gabriel Andrade

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <array>
#include <functional>
#include <iostream>
#include <utility>

//#include "openMVG/geometry/pose3.hpp"
#include "openMVG/multiview/triangulation.hpp"
#include "openMVG/multiview/triangulation_nview.hpp"
#include "openMVG/multiview/trifocal/solver_trifocal_metrics.hpp"
#include "openMVG/sfm/base/sfm_features_provider.hpp"
#include "openMVG/sfm/base/sfm_matches_provider.hpp"
#include "openMVG/sfm/base/SfM_Localizer.hpp"
#include "openMVG/sfm/base/sfm_robust_model_estimation.hpp"
#include "openMVG/sfm/base/sfm_data.hpp"
#include "openMVG/sfm/base/sfm_data_BA.hpp"
#include "openMVG/sfm/base/sfm_data_BA_ceres.hpp"
#include "openMVG/sfm/base/sfm_data_filters.hpp"
#include "openMVG/sfm/base/sfm_data_io.hpp"
#include "openMVG/stl/stl.hpp" 
#include "openMVG/system/logger.hpp" 
#include "openMVG/system/loggerprogress.hpp"

#include "third_party/histogram/histogram_raw.hpp"
#include "third_party/htmlDoc/htmlDoc.hpp"

#include <ceres/types.h>
#include "openMVG/sfm/pipelines/sequential/sequential_SfM.hpp"
#include "openMVG/sfm/pipelines/sequential/sequential_SfM_util.hpp"
#include "openMVG/sfm/pipelines/sequential/sfm_robust_model_estimation_trifocal.hpp"


#ifdef _MSC_VER
#pragma warning( once : 4267 ) //warning C4267: 'argument' : conversion from 'size_t' to 'const int', possible loss of data
#endif

namespace openMVG {
namespace sfm {

using namespace openMVG::cameras;
using namespace openMVG::geometry;
using namespace openMVG::matching;
using namespace histogramming;

//-------------------
//-- Incremental reconstruction
//-------------------
bool SequentialSfMReconstructionEngine::Process() 
{
  if (!InitLandmarkTracks()) return false;
  if (!MakeInitialSeedReconstruction()) return false;
  if (!ConsistencyCheck()) {
    OPENMVG_LOG_INFO << "Internal SfM Consistency check failure";
    return false;
  }

  if (!ResectOneByOneTilDone()) return false;
  FinalStatistics();
  return true;
}

inline static void TriangulateTangent2View
(
  const Mat3 &R0,
  const Vec3 &bearing0,
  const Vec3 &tangent0,
  const Mat3 &R1,
  const Vec3 &bearing1,
  const Vec3 &tangent1,
  Vec3 &Tangent
)
{
  Tangent = (R0.transpose()*tangent0.cross(bearing0)).cross(R1.transpose()*tangent1.cross(bearing1));
}

// Go through sfm_data and reconstructs all tangents
// At this stage, for robustness we use no TriangulateDLT for tangents.
// For each pair of views in the reconstruction
//    - reconstruct 3D tangents
//
// Input:
//  - Points and cameras already reconstructed
//    - Possibly, some points have not yet been reconstructed 
//    (not inliers / other filters removed them)
//  - Assume tracks connect all images 
//  (all reconstructed features visible in all views)
//  - 3D tangents will be computed only for those points that have already been
//  reconstructed
//  - 2D point observations assumed not to have orientation/tangent yet
//  (ObservationsInfo is not filled yet)
//
// Output:
// - 3D tangents
// - 2D tangents from feature provider stored in ObservationInfo
//
// See also
// bool track_triangulation in sfm_data_triangulation.cpp
// 
void SequentialSfMReconstructionEngine::ReconstructAllTangents()
{
  unsigned const nviews = sfm_data_.num_views() - set_remaining_view_id_.size();
  std::vector<const Observation *> ob(nviews); // use only nviews of this;
  std::vector<ObservationInfo *>obi(nviews);
  std::vector<IndexT> vi(nviews);
  std::vector<const Pose3 *> pose(nviews);
  std::vector<std::shared_ptr<cameras::Pinhole_Intrinsic> > intrinsics(nviews);

  for (auto &lit : sfm_data_.GetStructure()) {
    const Landmark &l = lit.second;  LandmarkInfo &li = sfm_data_.info[lit.first]; // creates info
    const Observations &obs = l.obs; ObservationsInfo &iobs = li.obs_info;

    //Observations::const_iterator iterObs[nviews];
    unsigned v = 0;
    for (auto &o : obs) {
      vi[v] = o.first;
      ob[v]  = &o.second;
      obi[v] = &li.obs_info[vi[v]];  // creates obs
      const features::SIOPointFeature *feature = &(features_provider_->sio_feats_per_view[vi[v]][ob[v]->id_feat]);
      double theta = feature->orientation();
      obi[v]->t = Vec2(std::cos(theta), std::sin(theta));
      pose[v] = &sfm_data_.poses.at(sfm_data_.GetViews().at(vi[v])->id_pose);

      intrinsics[v] = std::dynamic_pointer_cast<cameras::Pinhole_Intrinsic>(
         sfm_data_.GetIntrinsics().at(sfm_data_.GetViews().at(vi[v])->id_intrinsic));
      v++;
    }

    // determine the best two views to triangulate
    double best_angle = 0;
    IndexT best_v0 = 0;
    IndexT best_v1 = 0;
    for (IndexT v0= 0; v0 + 1 < nviews; ++v0)
      for (IndexT v1 = v0 + 1; v1 < nviews; ++v1) {
        const double angle = AngleBetweenRay(
          *pose[v0], intrinsics[v0].get(), *pose[v1], intrinsics[v1].get(), ob[v0]->x, ob[v1]->x);
        // We are chosing the biggest angle.
        // TODO choose the closest to 90 degrees
        // we could also require the tangent backprojection planes be as
        // orthogonal as possible
        if (angle > best_angle) {
          best_angle = angle;
          best_v0 = v0;
          best_v1 = v1;
        }
      }

   // reconstruct T from best_v0 and best_v1

   //- bearing: invert intrinsic
   
   Vec3 bearing0 = ((*intrinsics[best_v0])(ob[best_v0]->x));
   Vec3 bearing1 = ((*intrinsics[best_v0])(ob[best_v1]->x));

   Vec3 tangent0, tangent1;
   Pinhole_Intrinsic::invert_intrinsics_tgt(intrinsics[best_v0]->K(), obi[best_v0]->t.data(), tangent0.data());
   Pinhole_Intrinsic::invert_intrinsics_tgt(intrinsics[best_v1]->K(), obi[best_v1]->t.data(), tangent1.data());

   TriangulateTangent2View (
     pose[best_v0]->rotation(),
     bearing0,
     tangent0,
     pose[best_v1]->rotation(),
     bearing1,
     tangent1,
     li.T
   );
   li.T.normalize();
  } // end for each landmark
}

// checks sfm_data_ internal consistency and consistency with external sources
// such as provider
bool SequentialSfMReconstructionEngine::ConsistencyCheck() const
{
  OPENMVG_LOG_INFO << "Running consistency check";
  // For all landmarks
  //    - check X !=0

  //  for each observation
  //  Check the view of the observation hash matches any view id in the vieset
  //  Check id_feat points to a real feature with same .x OK
  unsigned const nviews_assumed = sfm_data_.num_views() - set_remaining_view_id_.size();
  for (const auto &lit : sfm_data_.GetStructure()) {
    OPENMVG_LOG_INFO << "A";
    const Landmark       &l = lit.second;
    const Observations &obs = l.obs;
    unsigned nviews = obs.size(); 
    assert(nviews == nviews_assumed);
    assert(l.X[0]);

    for (const auto &o : obs) {
      OPENMVG_LOG_INFO << "B";
      unsigned vi = o.first;
      const Observation &ob = o.second;
      Vec2 xf;
      const features::PointFeature *feature = &(features_provider_->feats_per_view[vi][ob.id_feat]);
      OPENMVG_LOG_INFO << "C";
      xf = feature->coords().cast<double>();

      if (features_provider_->has_sio_features()) {
        const features::PointFeature *siofeature = &(features_provider_->sio_feats_per_view[vi][ob.id_feat]);
        OPENMVG_LOG_INFO << "C";
        Vec2 xfs = siofeature->coords().cast<double>();
        assert((xfs - xf).norm() < 1e-9);
      }

      OPENMVG_LOG_INFO << "D";
      // OPENMVG_LOG_INFO << xf.transpose() << " ob: " << ob->x.transpose();
      assert((xf - ob.x).norm() < 1e-9);
      OPENMVG_LOG_INFO << "D";
    }
  }

  for (const auto &vit : sfm_data_.GetViews()) {
   OPENMVG_LOG_INFO << "E";
    assert(vit.first == vit.second->id_view);
  }

  // TODO check residuals are not huge

  return true;
}

bool SequentialSfMReconstructionEngine::ConsistencyCheckOriented() const
{
  ConsistencyCheck();

  OPENMVG_LOG_INFO << "Running oriented consitency check";
  assert(sfm_data_.is_oriented());
  for (const auto &lit : sfm_data_.GetStructure()) {
    const Landmark       &l = lit.second;

    OPENMVG_LOG_INFO << "A";
    assert(sfm_data_.info.count(lit.first));
    const LandmarkInfo &li = sfm_data_.GetInfo().at(lit.first); // [lit.first] but const
    const Observations &obs = l.obs;
    const ObservationsInfo &iobs = li.obs_info;
    unsigned nviews = obs.size(); 
    unsigned inviews = iobs.size();
    assert(inviews == nviews);
    assert(l.X.norm());

    for (const auto &o : obs) {
      OPENMVG_LOG_INFO << "B";
      unsigned vi = o.first;
      const Observation *ob = &o.second;
      const features::SIOPointFeature *feature = &(features_provider_->sio_feats_per_view[vi][ob->id_feat]);
      double theta = feature->orientation();
      OPENMVG_LOG_INFO << "C";
      assert(theta);
      Vec2 orient(std::cos(theta),std::sin(theta));
      assert(iobs.count(vi));
      assert((orient - iobs.at(vi).t).norm() < 1e-9);
      OPENMVG_LOG_INFO << "D";
    }
  }
  return true;
}

// Compute robust Resection of remaining images
// - group of images will be selected and resection + scene completion will be tried
bool SequentialSfMReconstructionEngine::ResectOneByOneTilDone()
{
  OPENMVG_LOG_INFO << "-------------------------------------------------------";
  OPENMVG_LOG_INFO << "Robust resection";
  size_t resectionGroupIndex = 0;
  std::vector<uint32_t> vec_possible_resection_indexes;
  while (FindImagesWithPossibleResection(vec_possible_resection_indexes)) {
    if (resection_method_ == resection::SolverType::P2Pt_FABBRI_ECCV12) {
      ReconstructAllTangents();
      if (!SequentialSfMReconstructionEngine::ConsistencyCheckOriented()) {
        OPENMVG_LOG_INFO << "Resection: Fail Oriented Test after 3d rec";
        return false;
      } else {
        OPENMVG_LOG_INFO << "Resection: Pass Oriented Test after 3d rec";
      }
    }

    // ------------------------------------------------------------------------
    // Add images to the 3D reconstruction
    bool bImageAdded = false;
    for (const auto & iter : vec_possible_resection_indexes) {
      bImageAdded |= Resection(iter); // <<------------------------------------
      set_remaining_view_id_.erase(iter);
    }
    // ------------------------------------------------------------------------
    if (bImageAdded) {
      // Scene logging as ply for visual debug
      std::ostringstream os;
      os << std::setw(8) << std::setfill('0') << resectionGroupIndex << "_Resection";
      Save(sfm_data_, stlplus::create_filespec(sOut_directory_, os.str(), ".ply"), ESfM_Data(ALL));

      // Perform BA until all point are under the given precision
      do BundleAdjustment(); while (badTrackRejector(4.0, 50));
      eraseUnstablePosesAndObservations(sfm_data_, 4); // XXX we are allowing 4
                                                       // points per pose as we
                                                       // are working with more
                                                       // radical pipeline
                                                       // situations
    }
    ++resectionGroupIndex;
  }
  // Ensure there is no remaining outliers
  //   if (badTrackRejector(4.0, 0))
  //     eraseUnstablePosesAndObservations(sfm_data_, 4);
  return true;
}

// Compute the initial 3D seed (First camera t=0; R=Id, second and third by
// Fabbri etal CVPR20 trifocal algorithm ). Computes the robust calibrated trifocal
// tensor / relative pose for ImageId [t[0],t[1],t[2]]
//
// Output
//  - sfm_data_
//  - map_ACThreshold_: if ACRansac used
//  - set_remaining_view_id_: remaining views to reconstruct
bool SequentialSfMReconstructionEngine::
MakeInitialTriplet3D(const Triplet &current_triplet)
{
  static constexpr unsigned nviews = 3;
  std::array<IndexT, nviews> t{ {std::get<0>(current_triplet),
                                 std::get<1>(current_triplet),
                                 std::get<2>(current_triplet)} };
  std::sort(t.begin(),t.end());

  for (unsigned v = 0; v < nviews; ++v)
    if (sfm_data_.GetViews().count(t[v]) == 0) {
      OPENMVG_LOG_ERROR << "Cannot find the view corresponding to the view id: " << t[v];
      return false;
    }

  // ---------------------------------------------------------------------------
  // a. Assert we have valid cameras
  const View *view[nviews] = {
    sfm_data_.GetViews().at(t[0]).get(),
    sfm_data_.GetViews().at(t[1]).get(),
    sfm_data_.GetViews().at(t[2]).get()
  };
    
  const Intrinsics::const_iterator iterIntrinsic[nviews] = {
    sfm_data_.GetIntrinsics().find(view[0]->id_intrinsic),
    sfm_data_.GetIntrinsics().find(view[1]->id_intrinsic),
    sfm_data_.GetIntrinsics().find(view[2]->id_intrinsic)
  };

  for (unsigned v=0; v < nviews; ++v)
    if (iterIntrinsic[v] == sfm_data_.GetIntrinsics().end()) {
      OPENMVG_LOG_ERROR << "Views with valid intrinsic data are required but this failed for view " << v;
      return false;
    }
  
  const cameras::IntrinsicBase *cam[nviews] = {
    iterIntrinsic[0]->second.get(),
    iterIntrinsic[1]->second.get(),
    iterIntrinsic[2]->second.get()
  };
  
  for (unsigned v=0; v < nviews; ++v) {
    if (!cam[v]) {
      OPENMVG_LOG_ERROR << "Cannot get back the camera intrinsic model for the triplet.";
      return false;
    }
    if (!dynamic_cast<const Pinhole_Intrinsic *>(cam[v])) {
      OPENMVG_LOG_ERROR << "Trifocal initialization only works for pinhole intrinsics K matrix.";
      return false;
    }
    OPENMVG_LOG_INFO << "K for v " << v << std::endl << dynamic_cast<const Pinhole_Intrinsic *>(cam[v])->K();
  }
  
  OPENMVG_LOG_INFO << "Putative starting triplet info\n\tindex:";
  for (unsigned v = 0; v < nviews; ++v)
    OPENMVG_LOG_INFO << t[v] << " ";
  OPENMVG_LOG_INFO << "\tview basename:";
  for (unsigned v = 0; v < nviews; ++v)
    OPENMVG_LOG_INFO << stlplus::basename_part(view[v]->s_Img_path) << " ";

  // ---------------------------------------------------------------------------
  // b. Get common features between the three views
  // use the track to have a more dense match correspondence set
  OPENMVG_LOG_INFO << "Geting common features between the three views";
  if (!features_provider_->has_sio_features()) {
    OPENMVG_LOG_ERROR << "Trifocal initialization only works for oriented features";
    return false;
  }
  openMVG::tracks::STLMAPTracks map_tracksCommon;
  shared_track_visibility_helper_->GetTracksInImages({t[0], t[1], t[2]}, map_tracksCommon);
  const size_t n = map_tracksCommon.size();
  OPENMVG_LOG_INFO << "number of tracks showing up in the three views = " << n ;
  std::array<Mat, nviews> pxdatum; // x,y,orientation across 3 views datum[view](coord,point)
  Mat scdatum;
  scdatum.resize(3,n);
  for (unsigned v = 0; v < nviews; ++v)
    pxdatum[v].resize(4,n);
  
  uint32_t cptIndex = 0;
  for (const auto &track_iter : map_tracksCommon) {
    auto iter = track_iter.second.cbegin();
    uint32_t i = iter->second;
    for (unsigned v = 0; v < nviews; ++v) {
      assert(features_provider_->sio_feats_per_view[t[v]].size());
      const features::SIOPointFeature *feature = &(features_provider_->sio_feats_per_view[t[v]][i]);
      pxdatum[v].col(cptIndex) << feature->x(), feature->y(), 
                                 cos(feature->orientation()), sin(feature->orientation());
      scdatum(v,cptIndex) = double(feature->scale());
      // TODO(trifocal future): provide undistortion for tangents for models that need it (get_ud_pixel)
      i=(++iter)->second;
    }
    ++cptIndex;
  }
  // ---------------------------------------------------------------------------
  // c. Robust estimation of the relative pose
  OPENMVG_LOG_INFO << "---------------------------------------------------------";
  OPENMVG_LOG_INFO << "Starting Trifocal robust estimation of the relative pose";
  OPENMVG_LOG_INFO << "---------------------------------------------------------";
  RelativePoseTrifocal_Info relativePose_info; // TODO(trifocal future): include image size
  if (!
      robustRelativePoseTrifocal(cam, pxdatum, relativePose_info, 4.0, maximum_trifocal_ransac_iterations_
     )) {
    OPENMVG_LOG_ERROR 
      << " /!\\ Robust estimation failed to compute calibrated trifocal tensor for this triplet: "
      << "{"<< t[0] << "," << t[1] << "," << t[2] << "}";
    return false;
  }
  //  OPENMVG_LOG_INFO 
  //    << "Trifocal Relative Pose residual from all inliers is: "
  //    << relativePose_info.found_residual_precision;
  // Bound min precision at 1 pix.
  relativePose_info.found_residual_precision = std::max(relativePose_info.found_residual_precision, 1.0);

  OPENMVG_LOG_INFO << "---------------------------------------------------------";
  OPENMVG_LOG_INFO << "Starting Bundle Adjustment for initial triplet";
  OPENMVG_LOG_INFO << "---------------------------------------------------------";

  SfM_Data tiny_scene;
  std::vector<Mat34> P;
  P.reserve(nviews);
  // tiny_scene.poses[view_J->id_pose] = relativePose_info.relativePose;
  for (unsigned v = 0; v < nviews; ++v) {
    // Init views and intrincics
    tiny_scene.views.insert(*sfm_data_.GetViews().find(view[v]->id_view));
    tiny_scene.intrinsics.insert(*iterIntrinsic[v]);
    OPENMVG_LOG_INFO << "Relative pose in _info \n" << relativePose_info.relativePoseTrifocal[v];
    if (v==0) 
      tiny_scene.poses[view[v]->id_pose] = Pose3(Mat3::Identity(), Vec3::Zero());
    else
      tiny_scene.poses[view[v]->id_pose] = Pose3(relativePose_info.relativePoseTrifocal[v].block<3,3>(0,0), 
                                                -relativePose_info.relativePoseTrifocal[v].block<3,3>(0,0).transpose()*relativePose_info.relativePoseTrifocal[v].block<3,1>(0,3));

    // Init projection matrices
    P.push_back(dynamic_cast<const Pinhole_Intrinsic *>(cam[v])->K()*(relativePose_info.relativePoseTrifocal[v]));
  }

  // OPENMVG_LOG_INFO << "\tscale[0]" << scdatum(0);

  // Init structure
  Landmarks &landmarks = tiny_scene.structure;
  { // initial structure ---------------------------------------------------
    Mat3X x(3,3);
    for (const auto &track_iterator : map_tracksCommon) {
      // Get corresponding points
      auto iter = track_iterator.second.cbegin();
      uint32_t ifeat = iter->second;
      for (unsigned v = 0; v < nviews; ++v) {
        x.col(v) =
          features_provider_->sio_feats_per_view[t[v]][ifeat].coords().cast<double>().homogeneous();
        landmarks[track_iterator.first].obs[view[v]->id_view] = Observation(x.col(v).hnormalized(), ifeat);
        // TODO(trifocal future) get_ud_pixel
        ifeat=(++iter)->second;
      }
      
      // triangulate 3 views
      Vec4 X;
      TriangulateNView(x, P, &X);
      landmarks[track_iterator.first].X = X.hnormalized();

      Vec2 residual = cam[0]->residual( tiny_scene.poses[view[0]->id_pose](landmarks[track_iterator.first].X), 
          landmarks[track_iterator.first].obs[view[0]->id_view].x );
      OPENMVG_LOG_INFO << "Residual from reconstructed point after robust-estimation " << residual.transpose();
      OPENMVG_LOG_INFO << "Residual from error()";
      { // For debug
      std::array<Mat, nviews> datum;
      for (unsigned v = 0; v < nviews; ++v) {
        datum[v].resize(4,1);
        // datum[v].col(0).head(2) = (*cam[v])(pxdatum[v].col(0).head<2>()).colwise().hnormalized(); // OK
        // datum[v].col(0).head(2) = (*cam[v])(x.col(v).hnormalized()).colwise().hnormalized();       // OK
        datum[v].col(0).head(2) = 
          (*cam[v])(landmarks[track_iterator.first].obs[view[v]->id_view].x).colwise().hnormalized(); // OK
      }
      OPENMVG_LOG_INFO << trifocal::NormalizedSquaredPointReprojectionOntoOneViewError::Error(
        relativePose_info.relativePoseTrifocal,
        datum[0].col(0), datum[1].col(0), datum[2].col(0));
      }
    }
    Save(tiny_scene, stlplus::create_filespec(sOut_directory_, "initialTriplet.ply"), ESfM_Data(ALL));
  } // !initial structure

  // -----------------------------------------------------------------------
  // - refine only Structure and Rotations & translations (keep intrinsic constant)
  constexpr bool bRefine_using_BA = false;
  if (bRefine_using_BA) { // badj
    Bundle_Adjustment_Ceres::BA_Ceres_options options(true, true);
    options.linear_solver_type_ = ceres::DENSE_SCHUR;
    Bundle_Adjustment_Ceres bundle_adjustment_obj(options);
    if (!bundle_adjustment_obj.Adjust(tiny_scene, Optimize_Options
        ( Intrinsic_Parameter_Type::NONE,           // Keep intrinsic constant
          Extrinsic_Parameter_Type::ADJUST_ALL,     // Adjust camera motion
          Structure_Parameter_Type::ADJUST_ALL))) { // Adjust structure
      return false;
    }
    std::cout << "passed BA\n";
  } // !badj
  // Save computed data
  Pose3 *pose[nviews];
  for (unsigned v = 0; v < nviews; ++v)  {
    sfm_data_.poses[view[v]->id_pose] = tiny_scene.poses[view[v]->id_pose];
    pose[v] = &tiny_scene.poses[view[v]->id_pose];
    map_ACThreshold_.insert({t[v], relativePose_info.found_residual_precision});
    set_remaining_view_id_.erase(view[v]->id_view);
    OPENMVG_LOG_INFO << "pose rc\n" << pose[v]->rotation() << "\n" <<  pose[v]->center().transpose();
  }

  OPENMVG_LOG_INFO << "After triplet BA, recompute inliers and save them";
  // Recompute inliers and save them
  // TODO: this is currently too strict,
  // every 2-view must pass
  for (const auto & landmark_entry : tiny_scene.GetLandmarks()) {
    const IndexT trackId = landmark_entry.first;
    const Landmark &landmark = landmark_entry.second;
    const Observations &obs = landmark.obs;
    OPENMVG_LOG_INFO << "\tTrack id " << trackId;

    Observations::const_iterator iterObs_x[nviews];
    const Observation *ob_x[nviews];
    Vec2 ob_x_ud[nviews];
    for (unsigned v = 0; v < nviews; ++v) {
      iterObs_x[v] = obs.find(view[v]->id_view);
      ob_x[v] = &iterObs_x[v]->second;
      ob_x_ud[v] = cam[v]->get_ud_pixel(ob_x[v]->x);

      // OPENMVG_LOG_INFO << "\t\tPoint in view " << v 
      // << " view id " << view[v]->id_view << " " << ob_x[v]->x << " = " << ob_x_ud[v];
    }
    bool include_landmark = true;
    for (unsigned v0 = 0; v0 + 1 < nviews; ++v0)
      for (unsigned v1 = v0 + 1; v1 < nviews; ++v1) {
        const double angle = AngleBetweenRay(
          *pose[v0], cam[v0], *pose[v1], cam[v1], ob_x_ud[v0], ob_x_ud[v1]);
        
        const Vec2 residual_0 = cam[v0]->residual((*pose[v0])(landmark.X), ob_x[v0]->x);
        const Vec2 residual_1 = cam[v1]->residual((*pose[v1])(landmark.X), ob_x[v1]->x);

        OPENMVG_LOG_INFO << "\t\tv0, v1 = " << v0 << ", " << v1 
          <<  "t\tresiduals norm " << residual_0.norm() << " " << residual_1.norm();
        if (angle <= 2.0) {
          OPENMVG_LOG_INFO << "\t\tFAIL angle test with angle " << angle;
          include_landmark = false;
        } else if (!CheiralityTest((*cam[v0])(ob_x_ud[v0]), *pose[v0],
                            (*cam[v1])(ob_x_ud[v1]), *pose[v1], landmark.X)) {
          OPENMVG_LOG_INFO << "\t\tFAIL Cheirality test ";
          include_landmark = false;
        } else if (residual_0.norm() >= relativePose_info.found_residual_precision ||
                   residual_1.norm() >= relativePose_info.found_residual_precision) {
            OPENMVG_LOG_INFO << "\t\tFAIL residual test: " << residual_0.norm() << " " 
              << residual_1.norm() << " both greater than " << relativePose_info.found_residual_precision;
          include_landmark = false;
        } // else if (UsingOrientedConstraint()) { } TODO 
      }
    if (include_landmark)
      sfm_data_.structure[trackId] = landmarks[trackId];
  }

  OPENMVG_LOG_INFO << "Final initial triplet residual histogram ---------------";
  // Save outlier residual information
  Histogram<double> histoResiduals;
  ComputeResidualsHistogram(&histoResiduals);
  if (!sLogging_file_.empty()) {
    using namespace htmlDocument;
    html_doc_stream_->pushInfo(htmlMarkup("h1","Trifocal tensor."));
    std::ostringstream os;
    os
      << "-------------------------------" << "<br>"
      << "-- Robust calibrated Trifocal tensor: <"  << t[0] << "," << t[1] << "," << t[2] << "> images: "
      << view[0]->s_Img_path << ","
      << view[1]->s_Img_path << ","
      << view[2]->s_Img_path << "<br>"
      << "-- Threshold: " << relativePose_info.found_residual_precision << "<br>"
      << "-- Resection status: " << "OK" << "<br>"
      << "-- Nb points used for robust Trifocal tensor estimation: "
      << pxdatum[0].cols() << "<br>"
      << "-- Nb points validated by robust estimation: "
      << sfm_data_.structure.size() << "<br>"
      << "-- % points validated: "
      << sfm_data_.structure.size()/static_cast<float>(pxdatum[0].cols())
      << "<br>"
      << "-------------------------------" << "<br>";
    html_doc_stream_->pushInfo(os.str());

    html_doc_stream_->pushInfo(htmlMarkup("h2",
      "Residual of the robust estimation (Initial triangulation). Thresholded at: "
      + toString(relativePose_info.found_residual_precision)));

    html_doc_stream_->pushInfo(htmlMarkup("h2","Histogram of residuals"));

    const std::vector<double> xBin = histoResiduals.GetXbinsValue();
    const auto range = autoJSXGraphViewport<double>(xBin, histoResiduals.GetHist());

    htmlDocument::JSXGraphWrapper jsxGraph;
    jsxGraph.init("InitialPairTriangulationKeptInfo",600,300);
    jsxGraph.addXYChart(xBin, histoResiduals.GetHist(), "line,point");
    jsxGraph.addLine(relativePose_info.found_residual_precision, 0,
      relativePose_info.found_residual_precision, histoResiduals.GetHist().front());
    jsxGraph.UnsuspendUpdate();
    jsxGraph.setViewport(range);
    jsxGraph.close();
    html_doc_stream_->pushInfo(jsxGraph.toStr());

    html_doc_stream_->pushInfo("<hr>");

    std::ofstream htmlFileStream( std::string(stlplus::folder_append_separator(sOut_directory_) +
      "Reconstruction_Report_Trifocal.html"));
    htmlFileStream << html_doc_stream_->getDoc();
  }

  return !sfm_data_.structure.empty();
}

/**
 * @brief Add one image to the 3D reconstruction. To the resectioning of
 * the camera and triangulate all the new possible tracks.
 * @param[in] viewIndex: image index to add to the reconstruction.
 *
 * A. Compute 2D/3D matches
 * B. Look if intrinsic data is known or not
 * C. Do the resectioning: compute the camera pose.
 * D. Refine the pose of the found camera
 * E. Update the global scene with the new camera
 * F. Update the observations into the global scene structure
 * G. Triangulate new possible 2D tracks
 */
bool SequentialSfMReconstructionEngine::Resection(const uint32_t viewIndex)
{
  using namespace tracks;

  // A. Compute 2D/3D matches
  // A1. list tracks ids used by the view
  openMVG::tracks::STLMAPTracks map_tracksCommon;
  shared_track_visibility_helper_->GetTracksInImages({viewIndex}, map_tracksCommon);
  std::set<uint32_t> set_tracksIds;
  TracksUtilsMap::GetTracksIdVector(map_tracksCommon, &set_tracksIds);

  // A2. intersects the track list with the reconstructed
  std::set<uint32_t> reconstructed_trackId;
  std::transform(sfm_data_.GetLandmarks().cbegin(), sfm_data_.GetLandmarks().cend(),
    std::inserter(reconstructed_trackId, reconstructed_trackId.begin()),
    stl::RetrieveKey());

  // Get the ids of the already reconstructed tracks
  std::set<uint32_t> set_trackIdForResection;
  std::set_intersection(set_tracksIds.cbegin(), set_tracksIds.cend(),
    reconstructed_trackId.cbegin(), reconstructed_trackId.cend(),
    std::inserter(set_trackIdForResection, set_trackIdForResection.begin()));

  if (set_trackIdForResection.empty()) {
    // No match. The image has no connection with already reconstructed points.
    OPENMVG_LOG_WARNING << "-- Failed to find the pose of the camera index: " << viewIndex
    << " ( the view have no connection with the already reconstructed 3d points)";
    return false;
  }

  // Get back featId associated to a tracksID already reconstructed.
  // These 2D/3D associations will be used for the resection.
  std::vector<uint32_t> vec_featIdForResection;
  TracksUtilsMap::GetFeatIndexPerViewAndTrackId(map_tracksCommon,
    set_trackIdForResection,
    viewIndex,
    &vec_featIdForResection);

  // Localize the image inside the SfM reconstruction
  Image_Localizer_Match_Data resection_data;
  resection_data.min_consensus_ratio = 1; // TODO make this a parameter
  resection_data.pt2D.resize(2, set_trackIdForResection.size());
  resection_data.pt3D.resize(3, set_trackIdForResection.size());
  resection_data.tgt2D.resize(2, set_trackIdForResection.size());
  resection_data.tgt3D.resize(3, set_trackIdForResection.size());
  if (resection_method_ == resection::SolverType::P2Pt_FABBRI_ECCV12) {
    resection_data.max_iteration = 1024; // TODO pose from only 2 points, can reduce a lot more
    assert(sfm_data_.is_oriented());
  }

  // B. Look if the intrinsic data is known or not
  const View *view_I = sfm_data_.GetViews().at(viewIndex).get();
  std::shared_ptr<cameras::IntrinsicBase> optional_intrinsic;
  if (sfm_data_.GetIntrinsics().count(view_I->id_intrinsic))
    optional_intrinsic = sfm_data_.GetIntrinsics().at(view_I->id_intrinsic);

  // Setup the track 2d observation for this new view
  Mat2X pt2D_original(2, set_trackIdForResection.size());
  Mat2X tgt2D_original(2, set_trackIdForResection.size());
  std::set<uint32_t>::const_iterator iterTrackId = set_trackIdForResection.begin();
  std::vector<uint32_t>::const_iterator iterfeatId = vec_featIdForResection.begin();


  for (size_t cpt = 0; cpt < vec_featIdForResection.size(); ++cpt, ++iterTrackId, ++iterfeatId) {
    resection_data.pt3D.col(cpt) = sfm_data_.GetLandmarks().at(*iterTrackId).X;
    if (features_provider_->has_sio_features()) {
      const features::SIOPointFeature *feature = 
        &(features_provider_->sio_feats_per_view.at(viewIndex)[*iterfeatId]);
      resection_data.pt2D.col(cpt) = pt2D_original.col(cpt) = feature->coords().cast<double>();
      double theta = feature->orientation();
      resection_data.tgt2D.col(cpt) = tgt2D_original.col(cpt) = Vec2(std::cos(theta), std::sin(theta));
    } else
      resection_data.pt2D.col(cpt) = pt2D_original.col(cpt) =
          features_provider_->feats_per_view.at(viewIndex)[*iterfeatId].coords().cast<double>();

    // Handle image distortion if intrinsic is known (to ease the resection)
    if (optional_intrinsic && optional_intrinsic->have_disto()) {
      resection_data.pt2D.col(cpt) = optional_intrinsic->get_ud_pixel(resection_data.pt2D.col(cpt));
      if (resection_method_ == resection::SolverType::P2Pt_FABBRI_ECCV12)
        OPENMVG_LOG_ERROR << "Not yet supported";
    }

    if (UsingOrientedConstraint() || resection_method_ == resection::SolverType::P2Pt_FABBRI_ECCV12)
      resection_data.tgt3D.col(cpt) = sfm_data_.GetInfo().at(*iterTrackId).T;

  }
  // C. Do the resectioning: compute the camera pose
  // TODO(p2pt)
  OPENMVG_LOG_INFO << "-- Trying robust Resection of view: " << viewIndex;

  geometry::Pose3 pose;
  const bool bResection = sfm::SfM_Localizer::Localize
  (
    optional_intrinsic ? resection_method_ : resection::SolverType::DLT_6POINTS,
    {view_I->ui_width, view_I->ui_height},
    optional_intrinsic.get(),
    resection_data,
    pose
  );
  resection_data.pt2D = std::move(pt2D_original); // restore original image domain points

  if (!sLogging_file_.empty())
  {
    using namespace htmlDocument;
    std::ostringstream os;
    os << "Resection of Image index: <" << viewIndex << "> image: "
      << view_I->s_Img_path <<"<br> \n";
    html_doc_stream_->pushInfo(htmlMarkup("h1",os.str()));

    os.str("");
    os
      << "-------------------------------" << "<br>"
      << "-- Robust Resection of camera index: <" << viewIndex << "> image: "
      <<  view_I->s_Img_path <<"<br>"
      << "-- Threshold: " << resection_data.error_max << "<br>"
      << "-- Resection status: " << (bResection ? "OK" : "FAILED") << "<br>"
      << "-- Nb points used for Resection: " << vec_featIdForResection.size() << "<br>"
      << "-- Nb points validated by robust estimation: " << resection_data.vec_inliers.size() << "<br>"
      << "-- % points validated: "
      << resection_data.vec_inliers.size()/static_cast<float>(vec_featIdForResection.size()) << "<br>"
      << "-------------------------------" << "<br>";
    html_doc_stream_->pushInfo(os.str());
  }

  if (!bResection)
    return false;

  // D. Refine the pose of the found camera.
  // We use a local scene with only the 3D points and the new camera.
  // TODO: are we just ignoring the inliers info and recomputing it?
  {
    const bool b_new_intrinsic = (optional_intrinsic == nullptr);
    // A valid pose has been found (try to refine it):
    // If no valid intrinsic as input:
    //  init a new one from the projection matrix decomposition
    // Else use the existing one and consider it as constant.
    if (b_new_intrinsic) {
      // setup a default camera model from the found projection matrix
      Mat3 K, R;
      Vec3 t;
      KRt_From_P(resection_data.projection_matrix, &K, &R, &t);

      const double focal = (K(0,0) + K(1,1))/2.0;
      const Vec2 principal_point(K(0,2), K(1,2));

      // Create the new camera intrinsic group
      switch (cam_type_) {
        case PINHOLE_CAMERA:
          optional_intrinsic =
            std::make_shared<Pinhole_Intrinsic>
            (view_I->ui_width, view_I->ui_height, focal, principal_point(0), principal_point(1));
        break;
        case PINHOLE_CAMERA_RADIAL1:
          optional_intrinsic =
            std::make_shared<Pinhole_Intrinsic_Radial_K1>
            (view_I->ui_width, view_I->ui_height, focal, principal_point(0), principal_point(1));
        break;
        case PINHOLE_CAMERA_RADIAL3:
          optional_intrinsic =
            std::make_shared<Pinhole_Intrinsic_Radial_K3>
            (view_I->ui_width, view_I->ui_height, focal, principal_point(0), principal_point(1));
        break;
        case PINHOLE_CAMERA_BROWN:
          optional_intrinsic =
            std::make_shared<Pinhole_Intrinsic_Brown_T2>
            (view_I->ui_width, view_I->ui_height, focal, principal_point(0), principal_point(1));
        break;
        case PINHOLE_CAMERA_FISHEYE:
            optional_intrinsic =
                std::make_shared<Pinhole_Intrinsic_Fisheye>
            (view_I->ui_width, view_I->ui_height, focal, principal_point(0), principal_point(1));
        break;
        default:
          OPENMVG_LOG_ERROR << "Try to create an unknown camera type (id):" << cam_type_;
          return false;
      }
    }
    const bool b_refine_pose = true;
    const bool b_refine_intrinsics = false;
    if (!
        sfm::SfM_Localizer::RefinePose(
        optional_intrinsic.get(), pose,
        resection_data, b_refine_pose, b_refine_intrinsics)
        ){
      OPENMVG_LOG_ERROR << "Unable to refine the pose of the view id: " << viewIndex;
      return false;
    }

    // E. Update the global scene with:
    // - the new found camera pose
    sfm_data_.poses[view_I->id_pose] = pose;
    // - track the view's AContrario robust estimation found threshold
    map_ACThreshold_.insert({viewIndex, resection_data.error_max});
    // - intrinsic parameters (if the view has no intrinsic group add a new one)
    if (b_new_intrinsic)
    {
      // Since the view have not yet an intrinsic group before, create a new one
      IndexT new_intrinsic_id = 0;
      if (!sfm_data_.GetIntrinsics().empty())
      {
        // Since some intrinsic Id already exists,
        //  we have to create a new unique identifier following the existing one
        std::set<IndexT> existing_intrinsicId;
          std::transform(sfm_data_.GetIntrinsics().cbegin(), sfm_data_.GetIntrinsics().cend(),
          std::inserter(existing_intrinsicId, existing_intrinsicId.begin()),
          stl::RetrieveKey());
        new_intrinsic_id = (*existing_intrinsicId.rbegin())+1;
      }
      sfm_data_.views.at(viewIndex)->id_intrinsic = new_intrinsic_id;
      sfm_data_.intrinsics[new_intrinsic_id] = optional_intrinsic;
    }
  }

  // F. List tracks that share content with this view and add observations and new 3D track if required.
  //    - If the track already exists (look if the new view tracks observation are valid)
  //    - If the track does not exists, try robust triangulation & add the new valid view track observation
  //    TODO: filter by tangent orientation
  {
    // Get information of new view
    const IndexT I = viewIndex;
    const View * view_I = sfm_data_.GetViews().at(I).get();
    const IntrinsicBase * cam_I = sfm_data_.GetIntrinsics().at(view_I->id_intrinsic).get();
    const Pose3 pose_I = sfm_data_.GetPoseOrDie(view_I);

    // Vector of all already reconstructed views
    const std::set<IndexT> valid_views = Get_Valid_Views(sfm_data_);

    // Go through each track and look if we must add new view observations or new 3D points
    for (const std::pair<uint32_t, tracks::submapTrack>& trackIt : map_tracksCommon) {
      const uint32_t trackId = trackIt.first;
      const tracks::submapTrack & track = trackIt.second;

      // List the potential view observations of the track
      const tracks::submapTrack & allViews_of_track = map_tracks_[trackId];

      // List to save the new view observations that must be added to the track
      std::set<IndexT> new_track_observations_valid_views;

      // If the track was already reconstructed
      if (sfm_data_.structure.count(trackId) != 0)
        // Since the 3D point was triangulated before we add the new the Inth view observation
        new_track_observations_valid_views.insert(I);
      else
        // Go through the views that observe this track & look if a successful triangulation can be done
        for (const std::pair<IndexT, IndexT>& trackViewIt : allViews_of_track) {
          const IndexT & J = trackViewIt.first;
          // If view is valid try triangulation
          if (J != I && valid_views.count(J) != 0) {
            // If successfully triangulated add the observation from J view
            if (sfm_data_.structure.count(trackId) != 0)
              new_track_observations_valid_views.insert(J);
            else {
              const View * view_J = sfm_data_.GetViews().at(J).get();
              const IntrinsicBase * cam_J = sfm_data_.GetIntrinsics().at(view_J->id_intrinsic).get();
              const Pose3 pose_J = sfm_data_.GetPoseOrDie(view_J);
              Vec2
              xJ = features_provider_->feats_per_view.at(J)[allViews_of_track.at(J)].coords().cast<double>(),
              xI = features_provider_->feats_per_view.at(I)[track.at(I)].coords().cast<double>();

              // Try to triangulate a 3D point from J view
              // A new 3D point must be added
              // Triangulate it
              const Vec2
                xI_ud = cam_I->get_ud_pixel(xI),
                xJ_ud = cam_J->get_ud_pixel(xJ);
              Vec3 X = Vec3::Zero();

              if (Triangulate2View(
                    pose_I.rotation(),
                    pose_I.translation(),
                    (*cam_I)(xI_ud),
                    pose_J.rotation(),
                    pose_J.translation(),
                    (*cam_J)(xJ_ud),
                    X,
                    triangulation_method_)) {
                // Check triangulation result
                const double angle = AngleBetweenRay(
                  pose_I, cam_I, pose_J, cam_J, xI_ud, xJ_ud);
                const Vec2 residual_I = cam_I->residual(pose_I(X), xI);
                const Vec2 residual_J = cam_J->residual(pose_J(X), xJ);
                if ( // TODO: track onrientation constraints
                    //  - Check angle (small angle leads to imprecise triangulation)
                    angle > 2.0 &&
                    //  - Check residual values (must be inferior to the found view's AContrario threshold)
                    residual_I.norm() < std::max(4.0, map_ACThreshold_.at(I)) &&
                    residual_J.norm() < std::max(4.0, map_ACThreshold_.at(J))
                    // Cheirality as been tested already in Triangulate2View
                   ) {
                  // Add a new track
                  Landmark & landmark = sfm_data_.structure[trackId];
                  landmark.X = X;
                  new_track_observations_valid_views.insert(I);
                  new_track_observations_valid_views.insert(J);
                } else { // 3D point is valid
                  OPENMVG_LOG_INFO << "XXX " << "not adding track in views " << I << ", " << J << std::endl;
                }
              } else {
                // We mark the view to add the observations once the point is triangulated
                new_track_observations_valid_views.insert(J);
              } // 3D point is invalid
            }
          }
        }// Go through all the views

      // If successfully triangulated, add the valid view observations
      if (sfm_data_.structure.count(trackId) != 0 && !new_track_observations_valid_views.empty()) {
        Landmark & landmark = sfm_data_.structure[trackId];
        // Check if view feature point observations of the track are valid (residual, depth) or not
        for (const IndexT & J: new_track_observations_valid_views) {
          const View * view_J = sfm_data_.GetViews().at(J).get();
          const IntrinsicBase * cam_J = sfm_data_.GetIntrinsics().at(view_J->id_intrinsic).get();
          const Pose3 pose_J = sfm_data_.GetPoseOrDie(view_J);
          const Vec2 xJ = features_provider_->feats_per_view.at(J)[allViews_of_track.at(J)].coords().cast<double>();
          const Vec2 xJ_ud = cam_J->get_ud_pixel(xJ);

          //    TODO: filter by tangent orientation
          const Vec2 residual = cam_J->residual(pose_J(landmark.X), xJ);
          if (CheiralityTest((*cam_J)(xJ_ud), pose_J, landmark.X)
              && residual.norm() < std::max(4.0, map_ACThreshold_.at(J)))
            landmark.obs[J] = Observation(xJ, allViews_of_track.at(J));
        }
      } // if
    }// All the tracks in the view
  }
  return true;
}

bool SequentialSfMReconstructionEngine::MakeInitialSeedReconstruction()
{
  // Initial pair Essential Matrix and [R|t] estimation.
  // or initial triplet relative pose
  //
  // Initial images choice
  // 
  // TODO(trifocal future) Refine this: If both initial triplet and initial pair are specified,
  // then first try intial pair then try initial triplet if that fails
  // OR try initial triplet first and initial pair if fails, which might be more
  // reliable anyways
  if (!hasInitialPair()) {
    if (!hasInitialTriplet()) {
      if (!AutomaticInitialPairChoice(initial_pair_)) {
        // Cannot find a valid initial pair with the defined settings:
        // TODO(trifocal) Trifocal - try to set it automatically
        // - try to initialize a pair with less strict constraint
        //    testing only X pairs with most matches.
        const auto sorted_pairwise_matches_iterators =
          GetPairWithMostMatches(sfm_data_, matches_provider_->pairWise_matches_, 20);
        for (const auto & it : sorted_pairwise_matches_iterators)
          if (MakeInitialPair3D({it->first.first, it->first.second})) {
            initial_pair_ = {it->first.first, it->first.second};
            break;
          }
        if (sorted_pairwise_matches_iterators.empty() || initial_pair_ == Pair(0,0)) {
          OPENMVG_LOG_INFO << "Cannot find a valid initial pair - stop reconstruction.";
          return false;
        }
      } else
          MakeInitialPair3D(initial_pair_);
    } else { // triplet but no pair
      OPENMVG_LOG_INFO << "Trying 3-view initialization from the provided one.";
      if (!MakeInitialTriplet3D(initial_triplet_)) {
          OPENMVG_LOG_INFO << "Tried 3-view initialization from the provided one, fail.";
          return false;
      }
    }
  } else { // has pair
    if (!MakeInitialPair3D(initial_pair_)) {
      OPENMVG_LOG_INFO << "Tried 2-view initialization from the provided one, fail.";
      if (!hasInitialTriplet())
        return false;
      OPENMVG_LOG_INFO << "Trying 3-view initialization from the provided one.";
      if (!MakeInitialTriplet3D(initial_triplet_)) {
        OPENMVG_LOG_INFO << "Tried 3-view initialization from the provided one, fail.";
        OPENMVG_LOG_INFO << "Tried 3-view initialization from the provided one, fail.";
        return false;
      }
    }
  }
  return true;
}


} // namespace sfm
} // namespace openMVG
