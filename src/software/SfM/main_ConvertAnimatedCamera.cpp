#include <openMVG/sfm/AlembicExporter.hpp>
#include <openMVG/sfm/sfm_data_io_gt.hpp>
#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"
#include <string>
#include <vector>

using namespace openMVG;
using namespace openMVG::sfm;

int main(int argc, char **argv)
{
  CmdLine cmd;

  std::string
    sSfM_Data_Filename_In,
    sSfM_Data_Filename_Out;

  cmd.add(make_option('i', sSfM_Data_Filename_In, "input_file"));
  cmd.add(make_option('o', sSfM_Data_Filename_Out, "output_file"));

  try {
      if (argc == 1) throw std::string("Invalid command line parameter.");
      cmd.process(argc, argv);
  } catch(const std::string& s) {
      std::cerr << "Usage: " << argv[0] << std::endl;
      std::cerr  << "[-i|--input_file] path to the input ground truth folder" << std::endl;
      std::cerr << "[-o|--output_file] path to the output Alembic file" << std::endl;

      std::cerr << s << std::endl;
      return EXIT_FAILURE;
  }

  if (sSfM_Data_Filename_In.empty() || sSfM_Data_Filename_Out.empty())
  {
    std::cerr << "Invalid input or output filename." << std::endl;
    return EXIT_FAILURE;
  }

  // Load input SfM_Data scene
  SfM_Data sfm_data;
  if (!readGt(sSfM_Data_Filename_In, sfm_data, false))
  {
    std::cerr << std::endl
      << "The input SfM_Data file \"" << sSfM_Data_Filename_In << "\" cannot be read." << std::endl;
    return EXIT_FAILURE;
  }

  // init alembic exporter
  dataio::AlembicExporter exporter( sSfM_Data_Filename_Out.c_str() );
  exporter.initAnimatedCamera("camera");

  for(const auto &iter : sfm_data.GetViews())
  {
    const auto &view = iter.second;
    const geometry::Pose3 pose_gt = sfm_data.GetPoses().at(view->id_pose);
    std::shared_ptr<cameras::IntrinsicBase> intrinsic_gt = std::make_shared<cameras::Pinhole_Intrinsic>();
    intrinsic_gt = sfm_data.GetIntrinsics().at(view->id_intrinsic);
    exporter.addCameraKeyframe(pose_gt, dynamic_cast<cameras::Pinhole_Intrinsic*>(intrinsic_gt.get()), view->s_Img_path, view->id_view, view->id_intrinsic);
  }
  return EXIT_SUCCESS;
}
