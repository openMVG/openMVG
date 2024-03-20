// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2024 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_SFM_DATA_IO_RERUN_HPP
#define OPENMVG_SFM_SFM_DATA_IO_RERUN_HPP

#include "openMVG/cameras/Camera_Pinhole.hpp"
#include "openMVG/image/image_io.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"

#include <rerun.hpp>

namespace openMVG {
namespace sfm {

inline void Save_Rerun(
  const rerun::RecordingStream & rec_stream,
  const SfM_Data& sfm_data,
  const ESfM_Data flags_part)
{
  const bool b_structure = (flags_part & STRUCTURE) == STRUCTURE;
  const bool b_extrinsics = (flags_part & EXTRINSICS) == EXTRINSICS;
  const bool b_images = (flags_part & VIEWS) == VIEWS;
  const bool b_structure_per_image = (flags_part & STRUCTURE_PER_VIEW) == STRUCTURE_PER_VIEW;

  const std::string camera_entity{"world/Camera/"};
  const std::string points3D_entity{"world/Points3D"};

  if (b_extrinsics)
  {
    for (const auto & view : sfm_data.GetViews())
    {
      if (sfm_data.GetIntrinsics().find(view.second->id_intrinsic) == sfm_data.GetIntrinsics().end())
        continue;
      if (sfm_data.GetPoses().find(view.second->id_pose) == sfm_data.GetPoses().end())
        continue;

      const auto view_file_name = stlplus::filename_part(view.second->s_Img_path);
      const geometry::Pose3 pose = sfm_data.GetPoseOrDie(view.second.get());
      const auto intrinsic = sfm_data.GetIntrinsics().find(view.second->id_intrinsic)->second;

      const Vec3f camera_position = pose.translation().cast<float>();
      const auto camera_orientation = pose.rotation().cast<float>();

      const rerun::datatypes::Mat3x3 rr_rotation{
            {camera_orientation(0, 0),
             camera_orientation(1, 0),
             camera_orientation(2, 0),
             camera_orientation(0, 1),
             camera_orientation(1, 1),
             camera_orientation(2, 1),
             camera_orientation(0, 2),
             camera_orientation(1, 2),
             camera_orientation(2, 2)}};

      const rerun::datatypes::Vec3D rr_translation{
            camera_position.x(),
            camera_position.y(),
            camera_position.z()};

      rec_stream.log(
            camera_entity + view_file_name,
            rerun::archetypes::Transform3D(rr_translation, rr_rotation, true));

      std::shared_ptr<openMVG::cameras::Pinhole_Intrinsic> pinhole_intrinsic(
          dynamic_cast<openMVG::cameras::Pinhole_Intrinsic*>(intrinsic->clone()));
      if (pinhole_intrinsic)
      {
        const rerun::datatypes::Vec2D resolution{
            static_cast<float>(intrinsic->w()),
            static_cast<float>(intrinsic->h())};

        rec_stream.log(camera_entity + view_file_name,
                rerun::archetypes::Pinhole::from_focal_length_and_resolution(
                  pinhole_intrinsic->focal(),
                  resolution));
      }
      else
      {
        std::cerr << "Camera type not yet supported in the sfm_data to Rerun exporter" << std::endl;
      }

      if (b_images)
      {
        openMVG::image::Image<openMVG::image::RGBColor> img;
        const auto image_name = stlplus::create_filespec(sfm_data.s_root_path, view_file_name);
        const auto is_img_loaded =
            openMVG::image::ReadImage(image_name.c_str(), &img);
        if (is_img_loaded) {
          rec_stream.log(camera_entity + view_file_name,
                  rerun::Image({static_cast<uint64_t>(img.rows()),
                                static_cast<uint64_t>(img.cols()),
                                static_cast<uint64_t>(img.Depth())},
                               img.GetMat().data()->data()));
        }
      }
    }
  }

  if (b_structure || b_structure_per_image)
  {
    std::vector<rerun::components::Position3D> points3d;
    std::vector<rerun::components::KeypointId> track_ids;
    const auto& landmarks = sfm_data.GetLandmarks();
    points3d.reserve(landmarks.size());
    track_ids.reserve(landmarks.size());
    std::unordered_map<uint32_t, std::vector<rerun::components::Position2D>>
        points2d_per_img;
    for (const auto& landmark : landmarks) {
      points3d.emplace_back(landmark.second.X(0), landmark.second.X(1), landmark.second.X(2));
      track_ids.push_back(landmark.first);
      if (b_structure_per_image)
      {
        for (const auto& obs : landmark.second.obs) {
          points2d_per_img[obs.first].emplace_back(
              static_cast<float>(obs.second.x(0)),
              static_cast<float>(obs.second.x(1)));
        }
      }
    }
    if (b_structure)
    {
      rec_stream.log(points3D_entity,
                     rerun::archetypes::Points3D(points3d).with_keypoint_ids(track_ids));
    }
    if (b_structure_per_image)
    {
      for (const auto& view : sfm_data.views) {
        auto it = points2d_per_img.find(view.first);
        if (it != points2d_per_img.end()) {
          rec_stream.log(camera_entity + stlplus::filename_part(view.second->s_Img_path),
                  rerun::archetypes::Points2D(points2d_per_img.at(view.first)));
        }
      }
    }
  }
}

/// Save the structure and camera positions of a SfM_Data container as a rerun rrd file.
inline bool Save_Rerun
(
  const SfM_Data & sfm_data,
  const std::string & filename,
  ESfM_Data flags_part
)
{
  auto rec_stream = rerun::RecordingStream("openMVG_to_rerun");
  const auto error = rec_stream.save(filename.c_str());
  Save_Rerun(rec_stream, sfm_data, flags_part);
  return error.is_ok();
}

} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_SFM_DATA_IO_RERUN_HPP
