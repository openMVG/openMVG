//
// Created by saurus on 18-7-26.
//

#include "sfm_data_PBA.hpp"
namespace openMVG {
  namespace sfm {
    void sfm_data_PBA::DataToPBA(SfM_Data &sfm_data_,
                                 cameras::Intrinsic_Parameter_Type intrinsic_refinement_options_,
                                 cameras::EINTRINSIC cam_type_) {
      if(intrinsic_refinement_options_ == cameras::Intrinsic_Parameter_Type::NONE) pba.SetFixedIntrinsics(true);
      else {
        if(cam_type_ == cameras::PINHOLE_CAMERA_RADIAL1)
          pba.EnableRadialDistortion(ParallelBA::DistortionT::PBA_PROJECTION_DISTORTION);
        else if(cam_type_ == cameras::PINHOLE_CAMERA_RADIAL1_PBA)
          pba.EnableRadialDistortion(ParallelBA::DistortionT::PBA_MEASUREMENT_DISTORTION);
        if((intrinsic_refinement_options_ & cameras::Intrinsic_Parameter_Type::ADJUST_DISTORTION)
           == static_cast<cameras::Intrinsic_Parameter_Type>(0))
          pba.EnableRadialDistortion(ParallelBA::DistortionT::PBA_NO_DISTORTION);
        if((intrinsic_refinement_options_ & cameras::Intrinsic_Parameter_Type::ADJUST_FOCAL_LENGTH)
           == static_cast<cameras::Intrinsic_Parameter_Type>(0))
          pba.SetFocalLengthFixed(true);
      }


      int i = 0;
      //Data to pba begin
      camera_data.resize(sfm_data_.poses.size());
      for (const auto &camera_openmvg : sfm_data_.poses) {
        camera_data[i] = CameraT();
        const Mat3 &camera_R = camera_openmvg.second.rotation();
        const Vec3 &camera_T = camera_openmvg.second.translation();
        for (int j = 0; j < 9; j++) camera_data[i].m[j / 3][j % 3] = (float) camera_R(j / 3, j % 3);
        for (int j = 0; j < 3; j++) camera_data[i].t[j] = (float) (camera_T(j));
        camera_data[i].f = static_cast<float>(sfm_data_.intrinsics[0].get()->getParams().at(0));
        if(cam_type_ == cameras::PINHOLE_CAMERA_RADIAL1_PBA)
          camera_data[i].SetMeasumentDistortion(sfm_data_.intrinsics[0].get()->getParams().at(3));
        else if(cam_type_ == cameras::PINHOLE_CAMERA_RADIAL1)
          camera_data[i].SetProjectionDistortion((sfm_data_.intrinsics[0].get()->getParams().at(3)));
        i++;
      }
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
      std::map<unsigned long, unsigned long> camera_map;
      camera_map.clear();
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
          double principal_x = sfm_data_.intrinsics[0].get()->getParams().at(1);
          double principal_y = sfm_data_.intrinsics[0].get()->getParams().at(2);
          measurements[tot++].SetPoint2D(obs.second.x.x() - principal_x, obs.second.x.y() - principal_y);
        }
        i++;
      }
      camera_map.clear(); sz = 0;
      sort(camera_id.begin(), camera_id.end());
      for (auto id : camera_id) if (camera_map.find(id) == camera_map.end()) camera_map[id] = sz++;
      for (auto &id : camidx) id = static_cast<int>(camera_map[id]);
      camera_id.clear(); camera_map.clear();

//    if you want to debug, you can export your data or load your own data use bellow
//    SaveBundlerModel("./test_output", camera_data, point_data, measurements, ptidx, camidx);
//    std::ifstream fin("./test_output");
//    LoadBundlerModel(fin, camera_data, point_data, measurements, ptidx, camidx);
      pba.SetCameraData(camera_data.size(), &camera_data[0]);
      focalmask = new int[sfm_data_.poses.size()];                                        //set mask to share model, mask must set after camera
      for(i = 0; i < sfm_data_.poses.size(); i++) focalmask[i] = 0;
      pba.SetFocalMask(focalmask, 1);

      pba.SetPointData(point_data.size(), &point_data[0]);                            //set 3D point data
      pba.SetProjection(measurements.size(), &measurements[0], &ptidx[0], &camidx[0]);//set the projections

      //Data to pba end
    }


    bool sfm_data_PBA::Adjust(SfM_Data &sfm_data_){
      pba.RunBundleAdjustment();
      //Data to openmvg start
      int i = 0;
      for (auto &camera_openmvg : sfm_data_.poses) {
        Mat3 camera_R;
        Vec3 camera_T;
        for (int j = 0; j < 9; j++) camera_R(j / 3, j % 3) = camera_data[i].m[j / 3][j % 3];
        for (int j = 0; j < 3; j++) camera_T(j) = camera_data[i].t[j];
        camera_openmvg.second = Pose3(camera_R, -camera_R.transpose() * camera_T);
        i++;
      }
      vector<double> params(4);
      params[0] = camera_data[0].f; params[1] =  sfm_data_.intrinsics[0].get()->getParams().at(1);
      params[2] = sfm_data_.intrinsics[0].get()->getParams().at(2); params[3] = camera_data[0].radial;
      sfm_data_.intrinsics[0].get()->updateFromParams(params);
      i = 0;
      for (auto &point_openmvg : sfm_data_.structure) {
        for (int j = 0; j < 3; j++) point_openmvg.second.X[j] = point_data[i].xyz[j];
        i++;
      }
      //Data to openmvg end
      return true;
    }
  }
}
