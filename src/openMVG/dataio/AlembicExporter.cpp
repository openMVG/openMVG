/* 
 * File:   Alembic_Exporter.cpp
 * Author: sgaspari
 * 
 * Created on September 24, 2015, 3:57 PM
 */

#if HAVE_ALEMBIC

#include "AlembicExporter.hpp"

namespace openMVG {
namespace dataio {


using namespace Alembic::Abc;
namespace AbcG = Alembic::AbcGeom;
using namespace AbcG;

AlembicExporter::~AlembicExporter()
{
}

void AlembicExporter::addPoints(const sfm::Landmarks &points)
{
  if(points.empty())
    return;
  // Fill vector with the values taken from OpenMVG 
  std::vector<V3f> positions;
  positions.reserve(points.size());

  // For all the 3d points in the hash_map
  for(const auto lm : points)
  {
    const openMVG::Vec3 &pt = lm.second.X;
    positions.emplace_back(pt[0], pt[1], pt[2]);
  }

  std::vector<Alembic::Util::uint64_t> ids(positions.size());
  iota(begin(ids), end(ids), 0);

  OPoints partsOut(topObj, "particleShape1");
  OPointsSchema &pSchema = partsOut.getSchema();

  OPointsSchema::Sample psamp(std::move(V3fArraySample(positions)), std::move(UInt64ArraySample(ids)));
  pSchema.set(psamp);
}

void AlembicExporter::appendCamera(const std::string &cameraName, 
                                   const geometry::Pose3 &pose,
                                   const cameras::Pinhole_Intrinsic *cam)
{
  const openMVG::Mat3 R = pose.rotation();
  const openMVG::Vec3 center = pose.center();
  // POSE
  // Compensate translation with rotation
  // Build transform matrix
  Abc::M44d xformMatrix;
  xformMatrix[0][0] = R(0, 0);
  xformMatrix[0][1] = R(0, 1);
  xformMatrix[0][2] = R(0, 2);
  xformMatrix[1][0] = R(1, 0);
  xformMatrix[1][1] = R(1, 1);
  xformMatrix[1][2] = R(1, 2);
  xformMatrix[2][0] = R(2, 0);
  xformMatrix[2][1] = R(2, 1);
  xformMatrix[2][2] = R(2, 2);
  xformMatrix[3][0] = center(0);
  xformMatrix[3][1] = center(1);
  xformMatrix[3][2] = center(2);
  xformMatrix[3][3] = 1.0;

  // Correct camera orientation for alembic
  M44d scale;
  scale[0][0] = 1;
  scale[1][1] = -1;
  scale[2][2] = -1;

  xformMatrix = scale*xformMatrix;

  XformSample xformsample;
  xformsample.setMatrix(xformMatrix);

  std::stringstream ss;
  ss << cameraName;
  Alembic::AbcGeom::OXform xform(topObj, "camxform_" + ss.str());
  xform.getSchema().set(xformsample);

  // Camera intrinsic parameters
  OCamera camObj(xform, "camera_" + ss.str());
  CameraSample camSample;

  // Take the max of the image size to handle the case where the image is in portrait mode 
  const float imgWidth = cam->w();
  const float imgHeight = cam->h();
  const float sensorWidthPix = std::max(imgWidth, imgHeight);
  const float sensorHeightPix = std::min(imgWidth, imgHeight);
  const float focalLengthPix = cam->focal();
  const float dx = cam->principal_point()(0);
  const float dy = cam->principal_point()(1);
  // Use a common sensor width as we don't have this information at this point
  // We chose a full frame 24x36 camera
  const float sensorWidth = 36.0; // 36mm per default
  const float sensorHeight = sensorWidth * sensorHeightPix / sensorWidthPix;
  const float focalLength = sensorWidth * focalLengthPix / sensorWidthPix;
  // Following values are in cm, hence the 0.1 multiplier
  const float hoffset = 0.1 * sensorWidth * (0.5 - dx / imgWidth);
  const float voffset = 0.1 * sensorHeight * (dy / imgHeight - 0.5) * sensorHeightPix / sensorWidthPix;
  const float haperture = 0.1 * sensorWidth * imgWidth / sensorWidthPix;
  const float vaperture = 0.1 * sensorWidth * imgHeight / sensorWidthPix;

  camSample.setFocalLength(focalLength);
  camSample.setHorizontalAperture(haperture);
  camSample.setVerticalAperture(vaperture);
  camSample.setHorizontalFilmOffset(hoffset);
  camSample.setVerticalFilmOffset(voffset);

  camObj.getSchema().set(camSample);
}

void AlembicExporter::addSfmData(const sfm::SfM_Data &sfm_data)
{
  for(const auto it : sfm_data.GetViews())
  {
    const sfm::View * view = it.second.get();
    if(!sfm_data.IsPoseAndIntrinsicDefined(view))
      continue;

    // OpenMVG Camera
    const openMVG::geometry::Pose3 pose = sfm_data.GetPoseOrDie(view);
    auto iterIntrinsic = sfm_data.GetIntrinsics().find(view->id_intrinsic);
    openMVG::cameras::Pinhole_Intrinsic *cam = static_cast<openMVG::cameras::Pinhole_Intrinsic*> (iterIntrinsic->second.get());
    
    const std::string cameraName = stlplus::basename_part(view->s_Img_path);
    appendCamera(cameraName, pose, cam);
  }
}

} //namespace dataio
} //namespace openMVG

#endif //HAVE_ALEMBIC