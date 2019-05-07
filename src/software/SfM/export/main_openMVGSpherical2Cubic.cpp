// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2019 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/cameras/Camera_Pinhole.hpp"
#include "openMVG/cameras/Camera_Spherical.hpp"
#include "openMVG/geometry/pose3.hpp"
#include "openMVG/image/image_io.hpp"
#include "openMVG/image/sample.hpp"
#include "openMVG/numeric/eigen_alias_definition.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"
#include "openMVG/sfm/sfm_landmark.hpp"
#include "openMVG/sfm/sfm_view.hpp"
#include "openMVG/types.hpp"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/progress/progress_display.hpp"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <array>
#include <cstdlib>
#include <cmath>
#include <iomanip>

using namespace openMVG;
using namespace openMVG::cameras;
using namespace openMVG::geometry;
using namespace openMVG::image;
using namespace openMVG::sfm;

/// Compute a rectilinear camera focal for a given angular desired FoV and image size
double FocalFromPinholeHeight
(
  int h,
  double thetaMax = openMVG::D2R(60) // Camera FoV
)
{
  float f = 1.f;
  while ( thetaMax < atan2( h / (2 * f) , 1))
  {
    ++f;
  }
  return f;
}

const std::array<openMVG::Mat3,6> GetCubicRotations()
{
  using namespace openMVG;
  return {
    RotationAroundY(D2R(0)),     // front
    RotationAroundY(D2R(-90)),   // right
    RotationAroundY(D2R(-180)),  // behind
    RotationAroundY(D2R(-270)),  // left
    RotationAroundX(D2R(-90)),   // up
    RotationAroundX(D2R(+90))    // down
  };
}

void ComputeCubicCameraIntrinsics(const int cubic_image_size,
                                  openMVG::cameras::Pinhole_Intrinsic & pinhole_camera)
{
  const double focal = FocalFromPinholeHeight(cubic_image_size, D2R(45));
  const double principal_point_xy = cubic_image_size / 2;
  pinhole_camera = Pinhole_Intrinsic(cubic_image_size,
                                     cubic_image_size,
                                     focal,
                                     principal_point_xy,
                                     principal_point_xy);
}

template <typename ImageT>
void SphericalToCubic
(
  const ImageT & equirectangular_image,
  const openMVG::cameras::Pinhole_Intrinsic & pinhole_camera,
  std::vector<ImageT> & cube_images
)
{
  using namespace openMVG;
  using namespace openMVG::cameras;

  const image::Sampler2d<image::SamplerLinear> sampler;

  //
  // Initialize a camera model for each image domain
  // - the equirectangular panorama
  const Intrinsic_Spherical sphere_camera(equirectangular_image.Width()  - 1,
                                          equirectangular_image.Height() - 1);
  // - the cube faces
  //
  // Perform backward/inverse rendering:
  // - For each cube face (rotation)
  // - Sample the panorama pixel by camera to camera bearing vector projection
  const int cubic_image_size = pinhole_camera.h();
  cube_images.resize(6, ImageT(cubic_image_size, cubic_image_size));

  // Initialize the rotation matrices corresponding to each cube face
  const std::array<Mat3, 6> rot_matrix = GetCubicRotations();

  for (const int i_rot : {0, 1, 2, 3, 4, 5})
  {
    auto & pinhole_image = cube_images[i_rot];
    const int image_width = pinhole_image.Width();
    const int image_height = pinhole_image.Height();

    // For every pinhole image pixels
    for (int x = 0; x < image_width; ++x)
    {
      for (int y = 0; y < image_height; ++y)
      { // Project the pinhole bearing vector to the spherical camera
        const Vec3 pinhole_bearing = rot_matrix[i_rot].transpose() * pinhole_camera(Vec2(x, y));
        const Vec2 sphere_proj = sphere_camera.project(pinhole_bearing);

        if (equirectangular_image.Contains(sphere_proj(1), sphere_proj(0)))
        {
          pinhole_image(y, x) = sampler(equirectangular_image, sphere_proj(1), sphere_proj(0));
        }
      }
    }
  }
}

int main(int argc, char *argv[]) {

  CmdLine cmd;
  std::string s_sfm_data_filename;
  std::string s_out_dir = "";
  int force_recompute_images = 1;

  cmd.add( make_option('i', s_sfm_data_filename, "sfmdata") );
  cmd.add( make_option('o', s_out_dir, "outdir") );
  cmd.add( make_option('f', force_recompute_images, "force_compute_cubic_images") );

  try {
      if (argc == 1) throw std::string("Invalid command line parameter.");
      cmd.process(argc, argv);
  } catch (const std::string& s) {
      std::cerr << "Usage: " << argv[0] << '\n'
      << "[-i|--sfmdata] filename, the SfM_Data file to convert\n"
      << "[-o|--outdir path]\n"
      << "[-f|--force_recompute_images] (default 1)\n"
      << std::endl;

      std::cerr << s << std::endl;
      return EXIT_FAILURE;
  }

  std::cout << "force_recompute_images = " << force_recompute_images << std::endl;

  // Create output dir
  if (!stlplus::folder_exists(s_out_dir))
      stlplus::folder_create(s_out_dir);

  SfM_Data sfm_data;
  if (!Load(sfm_data, s_sfm_data_filename, ESfM_Data(ALL))) {
      std::cerr << std::endl
      << "The input SfM_Data file \""<< s_sfm_data_filename << "\" cannot be read." << std::endl;
      return EXIT_FAILURE;
  }

  SfM_Data sfm_data_out; // the sfm_data that stores the cubical image list
  sfm_data_out.s_root_path = s_out_dir;

  // Convert every spherical view to cubic views
  {
    std::cout << "Generating cubic views:";
    C_Progress_display my_progress_bar(sfm_data.GetViews().size());
    const Views & views = sfm_data.GetViews();
    const Poses & poses = sfm_data.GetPoses();
    const Landmarks & structure = sfm_data.GetLandmarks();

    // generate views and camera poses for each new views
    int error_status = 0;
    #pragma omp parallel for shared(error_status) if(error_status < 1)
    for (int i = 0; i < static_cast<int>(views.size()); ++i)
    {
      ++my_progress_bar;
      auto view_it = views.begin();
      std::advance(view_it, i);
      const View * view = view_it->second.get();
      if (!sfm_data.IsPoseAndIntrinsicDefined(view))
          continue;

      Intrinsics::const_iterator iter_intrinsic = sfm_data.GetIntrinsics().find(view->id_intrinsic);

      const IntrinsicBase * cam = iter_intrinsic->second.get();
      if (cam && cam->getType() == CAMERA_SPHERICAL)
      {
        // We have a valid view with a corresponding camera & pose
        const std::string view_path = stlplus::create_filespec(sfm_data.s_root_path, view->s_Img_path);
        image::Image<image::RGBColor> spherical_image;
        if (!ReadImage(view_path.c_str(), &spherical_image))
        {
          std::cerr << "Cannot read the input panoramic image: " << view_path << std::endl;
          #pragma omp atomic
          ++error_status;
          continue;
        }

        openMVG::cameras::Pinhole_Intrinsic pinhole_camera;
        const int cubic_image_size = 1024;
        ComputeCubicCameraIntrinsics(cubic_image_size, pinhole_camera);

        // when cubical image computation is needed
        std::vector<image::Image<image::RGBColor>> cube_images;
        if (force_recompute_images)
          SphericalToCubic(spherical_image, pinhole_camera, cube_images);

        const std::array<Mat3,6> rot_matrix = GetCubicRotations();

        for (const int cubic_image_id : {0,1,2,3,4,5})
        {
          std::ostringstream os;
          os << std::setw(8) << std::setfill('0') << cubic_image_id;
          const std::string dst_cube_image = stlplus::create_filespec(
            stlplus::folder_append_separator(s_out_dir),
            stlplus::basename_part(view_path)
            + "_perspective_"
            + os.str(),
            "png");
          // when cubical image computation is needed
          if (force_recompute_images)
          {
            if (!WriteImage(dst_cube_image.c_str(), cube_images[cubic_image_id]))
            {
              std::cout << "Cannot export cubic images to: " << dst_cube_image << std::endl;
              #pragma omp atomic
              ++error_status;
              continue;
            }
          }

          const View v(
            stlplus::filename_part(dst_cube_image),
            view->id_view * 6 + cubic_image_id, // Id view
            0,  // Id intrinsic
            view->id_pose * 6 + cubic_image_id,  // Id pose
            pinhole_camera.w(),
            pinhole_camera.h());
          sfm_data_out.views[v.id_view] = std::make_shared<View>(v);

          // due to 360 cameras, rotation after BA might come with determinant -1
          // if so, negate the rotation for future use.
          Mat3 tmp_rotation = poses.at(view->id_pose).rotation();
          if (tmp_rotation.determinant() < 0)
          {
            std::cout << "Negative determinant" << std::endl;
            tmp_rotation = tmp_rotation*(-1.0f);
          }

          sfm_data_out.poses[v.id_pose] =
            Pose3(rot_matrix[cubic_image_id] * tmp_rotation,
                  poses.at(view->id_pose).center());

          if (sfm_data_out.GetIntrinsics().count(v.id_intrinsic) == 0)
            sfm_data_out.intrinsics[v.id_intrinsic] =
              std::make_shared<Pinhole_Intrinsic>(pinhole_camera);
      }
    }
    else
    {
      std::cout << "Loaded scene does not have spherical camera" << std::endl;
      #pragma omp atomic
      ++error_status;
      continue;
    }
  } // end of generate views and camera poses for each new views

  if (error_status > 0) // early exit
    return EXIT_FAILURE;

  // generate structure and associate it with new camera views
  {
    std::cout << "Creating cubic sfm_data structure:";

    C_Progress_display my_progress_bar(structure.size());
    for (const auto & it_structure : structure)
    {
      ++my_progress_bar;

      const Observations & obs = it_structure.second.obs;

      Landmark out_landmark;
      out_landmark.X = it_structure.second.X;

      // iterate across 360 views that can see the point
      for(const auto & it_obs : obs)
      {
        const IndexT pano_view_key = it_obs.first;
        const IndexT feature_key   = it_obs.second.id_feat;

        // get cubical camera ids and poses and reproject to see if the 3D point is inside the view
        bool is_reprojection_found = false;
        for (IndexT local_view_index = pano_view_key * 6; local_view_index < pano_view_key * 6 + 6; ++local_view_index)
        {
          const IndexT intrinsic_id = sfm_data_out.views[local_view_index]->id_intrinsic;
          const IndexT extrinsic_id = sfm_data_out.views[local_view_index]->id_pose;
          const Pose3 pose = sfm_data_out.poses[extrinsic_id];
          const auto cam = sfm_data_out.intrinsics[intrinsic_id];

          const int image_height = cam->h();
          const int image_width  = cam->w();

          const Vec2 projection = cam->project(pose(it_structure.second.X));

          if (projection.x() < 0 || projection.x() >= image_width ||
              projection.y() < 0 || projection.y() >= image_height)
              continue;

          const Vec3 point_to_cam_dir = (it_structure.second.X - pose.center()).normalized();
          const Vec3 cam_looking_dir = (pose.rotation().transpose() * Vec3(0, 0, 1)).normalized();

          const double angle = R2D(acos(point_to_cam_dir.dot(cam_looking_dir)));

          if (angle < 0 || angle > 90)
            continue;

          out_landmark.obs[local_view_index] = Observation(projection, feature_key);
          is_reprojection_found = true;
          break; // if one of the 6 views observe the 3D point, no other views from the 6 views should observe it
        } // end of looping 6 view of 1 pano image

        assert(is_reprojection_found); // make sure the observation is found
      } // end of observations from all panos
      sfm_data_out.structure.insert({it_structure.first, out_landmark});
    }
  }

  } // end of converting spherical view to cubical

  if (!Save(sfm_data_out,
            stlplus::create_filespec(stlplus::folder_append_separator(s_out_dir),
                                     "sfm_data_perspective.bin"),
            ESfM_Data(ALL))) {
    std::cerr << std::endl
    << "Cannot save the output sfm_data file" << std::endl;
    return EXIT_FAILURE;
  }

  std::cout
    << " #views: " << sfm_data_out.views.size() << "\n"
    << " #poses: " << sfm_data_out.poses.size() << "\n"
    << " #intrinsics: " << sfm_data_out.intrinsics.size() << "\n"
    << " #tracks: " << sfm_data_out.structure.size() << "\n" << std::endl;

  // Exit program
  return EXIT_SUCCESS;
}
