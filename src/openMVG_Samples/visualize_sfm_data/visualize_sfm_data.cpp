#include <rerun.hpp>

#include "openMVG/image/image_container.hpp"
#include "openMVG/image/image_io.hpp"
#include "openMVG/sfm/sfm.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"
#include "openMVG/sfm/sfm_data_utils.hpp"
#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

using namespace std::string_literals;

int main(int argc, char* argv[]) {
  CmdLine cmd;
  std::string input_sfm_data_path{};
  cmd.add(make_option('s', input_sfm_data_path, "sfm_data"));

  try {
    if (argc == 1) throw std::string("Invalid command line parameter.");
    cmd.process(argc, argv);
  } catch (const std::string& ex) {
    std::cerr << "Usage: " << argv[0] << '\n'
              << "[-s|--sfm_data] the path for the sfm_data in any "
              << "format openMVG::sfm::Save supports like .json or .bin\n";
    return EXIT_FAILURE;
  }

  // load sfm_data from a file.
  openMVG::sfm::SfM_Data sfm_data;
  std::string camera_entity{"world/Camera/"};
  std::string points3D_entity{"world/Points3D"};

  if (stlplus::is_file(input_sfm_data_path)) {
    openMVG::sfm::Load(sfm_data, input_sfm_data_path,
                       openMVG::sfm::ESfM_Data::ALL);
  } else {
    std::cerr << "The entered path is not valid.\n";
    return EXIT_FAILURE;
  }

  // rerun_sdk
  const auto rec = rerun::RecordingStream("openMVG_sfm_data_visualization");
  rec.spawn().exit_on_failure();

  const auto& views_poses = sfm_data.GetPoses();
  for (const auto& [view_id, view] : sfm_data.views) {
    /* be sure that the view->s_Img_path is the full path of the image when
    sfm_data got exported*/
    auto view_file_name = stlplus::filename_part(view->s_Img_path);
    const auto id_intrinsic = view->id_intrinsic;
    const auto id_pose = view->id_pose;

    if (auto it = views_poses.find(id_pose); it != views_poses.end()) {
      const auto& view_pose = views_poses.at(id_pose);
      Eigen::Vector3d view_translation = view_pose.translation();

      rerun::datatypes::Mat3x3 rr_rotation{
          {static_cast<float>(view_pose.rotation()(0, 0)),
           static_cast<float>(view_pose.rotation()(1, 0)),
           static_cast<float>(view_pose.rotation()(2, 0)),
           static_cast<float>(view_pose.rotation()(0, 1)),
           static_cast<float>(view_pose.rotation()(1, 1)),
           static_cast<float>(view_pose.rotation()(2, 1)),
           static_cast<float>(view_pose.rotation()(0, 2)),
           static_cast<float>(view_pose.rotation()(1, 2)),
           static_cast<float>(view_pose.rotation()(2, 2))}};

      rerun::datatypes::Vec3D rr_translation{
          static_cast<float>(view_translation(0)),
          static_cast<float>(view_translation(1)),
          static_cast<float>(view_translation(2))};

      rec.log(
          camera_entity + view_file_name,
          rerun::archetypes::Transform3D(rerun::datatypes::TranslationAndMat3x3(
              rr_translation, rr_rotation, true)));

      const rerun::datatypes::Vec2D resolution{
          static_cast<float>(view->ui_width),
          static_cast<float>(view->ui_height)};

      rec.log(camera_entity + view_file_name,
              rerun::archetypes::Pinhole::from_focal_length_and_resolution(
                  sfm_data.GetIntrinsics().at(id_intrinsic)->getParams()[0],
                  resolution));
      openMVG::image::Image<openMVG::image::RGBColor> img;
      const auto image_name = stlplus::create_filespec(sfm_data.s_root_path, view_file_name);
      auto is_img_loaded =
          openMVG::image::ReadImage(image_name.c_str(), &img);
      if (is_img_loaded) {
        rec.log(camera_entity + view_file_name,
                rerun::Image({static_cast<uint64_t>(img.rows()),
                              static_cast<uint64_t>(img.cols()),
                              static_cast<uint64_t>(img.Depth())},
                             img.GetMat().data()->data()));
      }
    }
  }

  const auto& landmarks = sfm_data.GetLandmarks();
  std::vector<rerun::components::Position3D> points3d;
  std::vector<rerun::components::KeypointId> track_ids;
  std::unordered_map<uint32_t, std::vector<rerun::components::Position2D>>
      points2d_per_img;
  for (const auto& landmark : landmarks) {
    points3d.emplace_back(landmark.second.X(0), landmark.second.X(1),
                          landmark.second.X(2));
    track_ids.push_back(landmark.first);
    for (const auto& obs : landmark.second.obs) {
      points2d_per_img[obs.first].push_back(
          {static_cast<float>(obs.second.x(0)),
           static_cast<float>(obs.second.x(1))});
    }
  }
  rec.log(points3D_entity,
          rerun::archetypes::Points3D(points3d).with_keypoint_ids(track_ids));

  for (const auto& view : sfm_data.views) {
    auto it = points2d_per_img.find(view.first);
    if (it != points2d_per_img.end()) {
      rec.log(camera_entity + stlplus::filename_part(view.second->s_Img_path),
              rerun::archetypes::Points2D(points2d_per_img.at(view.first)));
    }
  }
  return 0;
}