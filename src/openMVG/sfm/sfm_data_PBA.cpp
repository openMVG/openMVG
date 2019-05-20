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
      // cameras
      camera_data.resize(sfm_data_.poses.size());
      for (const auto &camera_openmvg : sfm_data_.poses) {
        camera_data[i] = CameraT();

        // extrinsics
        const Mat3 &camera_R = camera_openmvg.second.rotation();
        const Vec3 &camera_T = camera_openmvg.second.translation();
        for (int j = 0; j < 9; j++) camera_data[i].m[j / 3][j % 3] = (float) camera_R(j / 3, j % 3);
        for (int j = 0; j < 3; j++) camera_data[i].t[j] = (float) (camera_T(j));

        // intrinsics
        auto intrinsic_index = sfm_data_.views[camera_openmvg.first]->id_intrinsic;
        auto params = sfm_data_.intrinsics[intrinsic_index]->getParams();
        camera_data[i].f = (float)params[0];
        if (distortion_type == ParallelBA::DistortionT::PBA_PROJECTION_DISTORTION) {
          camera_data[i].SetProjectionDistortion(params[3]);
        } else if (distortion_type == ParallelBA::DistortionT::PBA_MEASUREMENT_DISTORTION) {
          camera_data[i].SetMeasumentDistortion(params[3]);
        }
        i++;
      }

      // points
      i = 0;
      point_data.resize(sfm_data_.structure.size());
      for (const auto &structure_openmvg : sfm_data_.structure) {
        double temp[3];
        for(int j = 0; j < 3; j++) temp[j] = structure_openmvg.second.X[j];
        point_data[i].SetPoint(temp);
        i++;
      }

      //transform camera and point id
      //example: 1 3 5 9 --> 0 1 2 3
      camera_map.clear();
      std::map<unsigned long, unsigned long> camera_id2pose_id;
      camera_id2pose_id.clear();
      vector<unsigned long> camera_id;
      camera_id.clear();

      unsigned long sz = 0;
      for(const auto &structure_openmvg : sfm_data_.structure) sz += structure_openmvg.second.obs.size();
      measurements.resize(sz);
      camidx.resize(sz);
      ptidx.resize(sz);

      i = 0; sz = 0;
      int tot = 0;
      for (const auto &structure_openmvg : sfm_data_.structure) {
        for (auto &obs : structure_openmvg.second.obs) {
          camidx[tot] = obs.first;
          if (camera_map.find(obs.first) == camera_map.end()) {
            camera_id.push_back(obs.first);
            camera_map[obs.first] = sz++;
          }
          ptidx[tot] = i;
          auto intrinsic_index = sfm_data_.views[obs.first]->id_intrinsic;
          auto params = sfm_data_.intrinsics[intrinsic_index]->getParams();
          double principal_x = params[1];
          double principal_y = params[2];
          measurements[tot++].SetPoint2D(obs.second.x.x() - principal_x, obs.second.x.y() - principal_y);
        }
        i++;
      }
      camera_map.clear(); sz = 0;
      sort(camera_id.begin(), camera_id.end());
      for (auto id : camera_id) if (camera_map.find(id) == camera_map.end())
      {
        camera_map[id] = sz;
        camera_id2pose_id[sz] = id;
        ++sz;
      }
      // init focalmask
      focalmask = std::vector<int>(sfm_data_.poses.size());                           //set mask to share model, mask must set after camera
      for (const auto &pose : sfm_data_.poses) {
        focalmask[camera_map[pose.first]] = sfm_data_.views[pose.first]->id_intrinsic;
      }
      // reorder focalmask, focalmask's first n+1 elements should be 0-n, n = max(focalmask)
      int max_ele = *std::max_element(focalmask.begin(), focalmask.end());
      for (i = 0; i <= max_ele; ++i) {
        if (focalmask[i] == i)
          continue;
        auto it = std::find(focalmask.begin() + i, focalmask.end(), i);
        if (it == focalmask.end())
          return false;

        int j = it - focalmask.begin();
        std::iter_swap(focalmask.begin() + i, it);
        std::iter_swap(camera_data.begin() + i, camera_data.begin() + j);
        std::swap(camera_map[camera_id2pose_id[i]], camera_map[camera_id2pose_id[j]]);
      }
      for (auto &id : camidx) id = static_cast<int>(camera_map[id]);
      camera_id.clear(); camera_id2pose_id.clear();

//    if you want to debug, you can export your data or load your own data use bellow
//    SaveBundlerModel("./test_output", camera_data, point_data, measurements, ptidx, camidx);
//    std::ifstream fin("./test_output");
//    LoadBundlerModel(fin, camera_data, point_data, measurements, ptidx, camidx);
      pba.SetCameraData(camera_data.size(), &camera_data[0]);
      pba.SetFocalMask(focalmask.data(), 1);

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
        int camera_id = camera_map[camera_openmvg.first];
        for (int j = 0; j < 9; j++) camera_R(j / 3, j % 3) = camera_data[camera_id].m[j / 3][j % 3];
        for (int j = 0; j < 3; j++) camera_T(j) = camera_data[camera_id].t[j];
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
        for (int j = 0; j < 3; j++) point_openmvg.second.X[j] = point_data[i].xyz[j];
        i++;
      }
      //Data to openmvg end
      return true;
    }
  }
}
