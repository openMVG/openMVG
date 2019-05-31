//
// Created by saurus on 18-7-26.
//

#include "sfm_data_PBA.hpp"
namespace openMVG {
  namespace sfm {
    bool sfm_data_PBA::DataToPBA(SfM_Data &sfm_data_,
                                 cameras::Intrinsic_Parameter_Type intrinsic_refinement_options_) {
      if (sfm_data_.intrinsics.size() == 0) {
        std::cerr << "SfMEnginePBA: Intrinsics number is 0" << std::endl;
        return false;
      }

      auto cam_type = sfm_data_.intrinsics.begin()->second->getType();
      using intrinsic_type = typename Hash_Map<IndexT, std::shared_ptr<cameras::IntrinsicBase>>::const_reference;
      bool single_camera_type = std::all_of(sfm_data_.intrinsics.begin(), sfm_data_.intrinsics.end(), [cam_type](intrinsic_type item) {
        return item.second->getType() == cam_type;
      });
      if (!single_camera_type)
      {
        std::cerr << "SfMEnginePBA: Only support one camera type" << std::endl;
        return false;
      }

      auto distortion_type = ParallelBA::DistortionT::PBA_NO_DISTORTION;
      if (cam_type == cameras::PINHOLE_CAMERA_RADIAL1) {
        distortion_type = ParallelBA::DistortionT::PBA_PROJECTION_DISTORTION;
      } else if (cam_type == cameras::PINHOLE_CAMERA_RADIAL1_PBA) {
        distortion_type = ParallelBA::DistortionT::PBA_MEASUREMENT_DISTORTION;
      } else {
        std::cerr << "SfMEnginePBA: Only support PINHOLE_CAMERA_RADIAL1 and PINHOLE_CAMERA_RADIAL1_PBA" << std::endl;
        return false;
      }

      if (intrinsic_refinement_options_ == cameras::Intrinsic_Parameter_Type::NONE) {
        pba.SetFixedIntrinsics(true);
      } else {
        if ((intrinsic_refinement_options_ & cameras::Intrinsic_Parameter_Type::ADJUST_DISTORTION)
            == (cameras::Intrinsic_Parameter_Type)0) {
          distortion_type = ParallelBA::DistortionT::PBA_NO_DISTORTION;
        }

        if ((intrinsic_refinement_options_ & cameras::Intrinsic_Parameter_Type::ADJUST_FOCAL_LENGTH)
            == (cameras::Intrinsic_Parameter_Type)0) {
          pba.SetFocalLengthFixed(true);
        }
      }
      pba.EnableRadialDistortion(distortion_type);

      int i = 0;
      //Data to pba begin
      // Calculate map from sfm_data_.poses to camera_data,
      // So that first n intrinsics in camera_data are non identical (n = different intrinsic num in sfm_data.poses)
      // (Required by pba)
      std::unordered_map<IndexT, IndexT> unique_intrinsic;        // from intrinsic_id to pose_id
      std::vector<std::pair<IndexT, IndexT>> duplicate_intrinsic; // pair: intrinsic_id, pose_id
      for (const auto &camera_openmvg : sfm_data_.poses) {
        auto intrinsic_id = sfm_data_.views[camera_openmvg.first]->id_intrinsic;
        auto pose_id = camera_openmvg.first;
        if (unique_intrinsic.count(intrinsic_id) == 0)
          unique_intrinsic[intrinsic_id] = pose_id;
        else
          duplicate_intrinsic.emplace_back(intrinsic_id, pose_id);
      }
      i = 0;
      pose_id2camera_id.clear();
      for (const auto &ele : unique_intrinsic)
        pose_id2camera_id[ele.second] = i++;
      for (const auto &pair : duplicate_intrinsic)
        pose_id2camera_id[pair.second] = i++;
      unique_intrinsic.clear();
      duplicate_intrinsic.clear();

      // cameras
      camera_data.resize(sfm_data_.poses.size());
      for (const auto &camera_openmvg : sfm_data_.poses) {
        auto idx = pose_id2camera_id[camera_openmvg.first];
        camera_data[idx] = CameraT();

        // extrinsics
        const Mat3 &camera_R = camera_openmvg.second.rotation();
        const Vec3 &camera_T = camera_openmvg.second.translation();
        for (int j = 0; j < 9; ++j) camera_data[idx].m[j / 3][j % 3] = (float) camera_R(j / 3, j % 3);
        for (int j = 0; j < 3; ++j) camera_data[idx].t[j] = (float) (camera_T(j));

        // intrinsics
        auto intrinsic_index = sfm_data_.views[camera_openmvg.first]->id_intrinsic;
        auto params = sfm_data_.intrinsics[intrinsic_index]->getParams();
        camera_data[idx].f = (float)params[0];
        if (distortion_type == ParallelBA::DistortionT::PBA_PROJECTION_DISTORTION) {
          camera_data[idx].SetProjectionDistortion(params[3]);
        } else if (distortion_type == ParallelBA::DistortionT::PBA_MEASUREMENT_DISTORTION) {
          camera_data[idx].SetMeasumentDistortion(params[3]);
        }
      }

      // points
      i = 0;
      point_data.resize(sfm_data_.structure.size());
      for (const auto &structure_openmvg : sfm_data_.structure) {
        double temp[3];
        for(int j = 0; j < 3; ++j) temp[j] = structure_openmvg.second.X[j];
        point_data[i].SetPoint(temp);
        ++i;
      }

      unsigned long sz = 0;
      for(const auto &structure_openmvg : sfm_data_.structure) sz += structure_openmvg.second.obs.size();
      measurements.resize(sz);
      camidx.resize(sz);
      ptidx.resize(sz);

      i = 0; sz = 0;
      int tot = 0;
      for (const auto &structure_openmvg : sfm_data_.structure) {
        for (auto &obs : structure_openmvg.second.obs) {
          camidx[tot] = pose_id2camera_id[obs.first];
          ptidx[tot] = i;
          auto intrinsic_index = sfm_data_.views[obs.first]->id_intrinsic;
          auto params = sfm_data_.intrinsics[intrinsic_index]->getParams();
          double principal_x = params[1];
          double principal_y = params[2];
          measurements[tot++].SetPoint2D(obs.second.x.x() - principal_x, obs.second.x.y() - principal_y);
        }
        ++i;
      }

      // focalmask
      // set mask to share model
      focalmask.resize(sfm_data_.poses.size());
      for (const auto &pose : sfm_data_.poses) {
        focalmask[pose_id2camera_id[pose.first]] = sfm_data_.views[pose.first]->id_intrinsic;
      }
      // transform focalmask from id_intrinsic --> 0 1 2 3
      i = 0;
      std::unordered_map<int, int> focalmask_map;
      for (auto &v : focalmask)
      {
        if (focalmask_map.count(v) == 0) {
          focalmask_map[v] = i;
          v = i++;
        } else {
          v = focalmask_map[v];
        }
      }
      focalmask_map.clear();

//    if you want to debug, you can export your data or load your own data use bellow
//    SaveBundlerModel("./test_output", camera_data, point_data, measurements, ptidx, camidx);
//    std::ifstream fin("./test_output");
//    LoadBundlerModel(fin, camera_data, point_data, measurements, ptidx, camidx);
      pba.SetCameraData(camera_data.size(), &camera_data[0]);
      pba.SetFocalMask(&focalmask[0], 1);

      pba.SetPointData(point_data.size(), &point_data[0]);                            //set 3D point data
      pba.SetProjection(measurements.size(), &measurements[0], &ptidx[0], &camidx[0]);//set the projections

      //Data to pba end
      return true;
    }


    bool sfm_data_PBA::Adjust(SfM_Data &sfm_data_){
      pba.RunBundleAdjustment();
      //Data to openmvg start
      for (auto &camera_openmvg : sfm_data_.poses) {
        // extrinsics
        Mat3 camera_R;
        Vec3 camera_T;
        int camera_id = pose_id2camera_id[camera_openmvg.first];
        for (int j = 0; j < 9; ++j) camera_R(j / 3, j % 3) = camera_data[camera_id].m[j / 3][j % 3];
        for (int j = 0; j < 3; ++j) camera_T(j) = camera_data[camera_id].t[j];
        camera_openmvg.second = Pose3(camera_R, -camera_R.transpose() * camera_T);

        // intrinsics
        auto intrinsic = sfm_data_.intrinsics[sfm_data_.views[camera_openmvg.first]->id_intrinsic];
        auto params = intrinsic->getParams();
        params[0] = camera_data[camera_id].f;
        params[3] = camera_data[camera_id].radial;
        intrinsic->updateFromParams(params);
      }

      // points
      int i = 0;
      for (auto &point_openmvg : sfm_data_.structure) {
        for (int j = 0; j < 3; ++j) point_openmvg.second.X[j] = point_data[i].xyz[j];
        ++i;
      }
      //Data to openmvg end
      return true;
    }
  }
}
