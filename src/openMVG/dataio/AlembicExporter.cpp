/* 
 * File:   Alembic_Exporter.cpp
 * Author: sgaspari
 * 
 * Created on September 24, 2015, 3:57 PM
 */

#if HAVE_ALEMBIC

#include "AlembicExporter.hpp"
#include "openMVG/sfm/sfm_view_metadata.hpp"

namespace openMVG {
namespace dataio {


using namespace Alembic::Abc;
namespace AbcG = Alembic::AbcGeom;
using namespace AbcG;

AlembicExporter::AlembicExporter(const std::string &filename)
: archive(Alembic::AbcCoreHDF5::WriteArchive(), filename)
, topObj(archive, Alembic::Abc::kTop)
, counter(0)
{
  // Create MVG hierarchy
  mvgRoot = Alembic::Abc::OObject(topObj, "mvgRoot");
  mvgCameras = Alembic::Abc::OObject(mvgRoot, "mvgCameras");
  mvgCloud = Alembic::Abc::OObject(mvgRoot, "mvgCloud");
  mvgPointCloud = Alembic::Abc::OObject(mvgCloud, "mvgPointCloud"); 

}

AlembicExporter::~AlembicExporter()
{
}

void AlembicExporter::addPoints(const sfm::Landmarks &landmarks, bool withVisibility)
{
  if(landmarks.empty())
    return;

  // Fill vector with the values taken from OpenMVG 
  std::vector<V3f> positions;
  positions.reserve(landmarks.size());

  // For all the 3d points in the hash_map
  for(const auto landmark : landmarks)
  {
    const openMVG::Vec3 &pt = landmark.second.X;
    positions.emplace_back(pt[0], pt[1], pt[2]);
  }

  std::vector<Alembic::Util::uint64_t> ids(positions.size());
  iota(begin(ids), end(ids), 0);

  OPoints partsOut(mvgPointCloud, "particleShape1");
  OPointsSchema &pSchema = partsOut.getSchema();

  OPointsSchema::Sample psamp(std::move(V3fArraySample(positions)), std::move(UInt64ArraySample(ids)));
  pSchema.set(psamp);

  if(withVisibility)
  {
    OCompoundProperty userProps = pSchema.getUserProperties();
    std::vector<uint32_t> visibilitySize;
    visibilitySize.reserve(positions.size());
    for(const auto landmark : landmarks)
    {
      visibilitySize.emplace_back(landmark.second.obs.size());
    }
    std::size_t nbObservations = std::accumulate(visibilitySize.begin(), visibilitySize.end(), 0);
    
    // Use std::vector<uint32_t> and std::vector<float> instead of std::vector<V2i> and std::vector<V2f>
    // Because Maya don't import them correctly
    std::vector<uint32_t> visibilityIds;
    visibilityIds.reserve(nbObservations*2);
    std::vector<float>featPos2d;
    featPos2d.reserve(nbObservations*2);

    for(sfm::Landmarks::const_iterator itLandmark = landmarks.cbegin(), itLandmarkEnd = landmarks.cend();
       itLandmark != itLandmarkEnd; ++itLandmark)
    {
      const sfm::Observations& observations = itLandmark->second.obs;
      for(const auto vObs: observations )
      {
        const sfm::Observation& obs = vObs.second;
        // (View ID, Feature ID)
        visibilityIds.emplace_back(vObs.first);
        visibilityIds.emplace_back(obs.id_feat);
        // Feature 2D position (x, y))
        featPos2d.emplace_back(obs.x[0]);
        featPos2d.emplace_back(obs.x[1]);
      }
    }

    OUInt32ArrayProperty propVisibilitySize( userProps, "mvg_visibilitySize" );
    propVisibilitySize.set(visibilitySize);

    // (viewID, featID)
    OUInt32ArrayProperty propVisibilityIds( userProps, "mvg_visibilityIds" );
    propVisibilityIds.set(visibilityIds);

    // Feature position (x,y)
    OFloatArrayProperty propFeatPos2d( userProps, "mvg_visibilityFeatPos" );
    propFeatPos2d.set(featPos2d);
  }
}

void AlembicExporter::appendCamera(const std::string &cameraName, 
                                   const geometry::Pose3 &pose,
                                   const cameras::Pinhole_Intrinsic *cam,
                                   const std::string &imagePath,
                                   const IndexT id_view,
                                   const IndexT id_intrinsic,
                                   const float sensorWidth_mm)
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
  M44d scale;   //  by default this is an identity matrix
  scale[0][0] = 1;
  scale[1][1] = -1;
  scale[2][2] = -1;

  xformMatrix = scale*xformMatrix;

  XformSample xformsample;
  xformsample.setMatrix(xformMatrix);

  std::stringstream ss;
  ss << cameraName;
  Alembic::AbcGeom::OXform xform(mvgCameras, "camxform_" + ss.str());
  xform.getSchema().set(xformsample);

  // Camera intrinsic parameters
  OCamera camObj(xform, "camera_" + ss.str());
  auto userProps = camObj.getSchema().getUserProperties();
  CameraSample camSample;

  // Take the max of the image size to handle the case where the image is in portrait mode 
  const float imgWidth = cam->w();
  const float imgHeight = cam->h();
  const float sensorWidth_pix = std::max(imgWidth, imgHeight);
  const float sensorHeight_pix = std::min(imgWidth, imgHeight);
  const float imgRatio = sensorHeight_pix / sensorWidth_pix;
  const float focalLength_pix = cam->focal();
  const float hoffset_pix = cam->principal_point()(0);
  const float voffset_pix = cam->principal_point()(1);
  

  const float sensorHeight_mm = sensorWidth_mm * imgRatio;
  const float focalLength_mm = sensorWidth_mm * focalLength_pix / sensorWidth_pix;
  const float pix2mm = sensorWidth_mm / sensorWidth_pix;

  // openMVG: origin is (top,left) corner and orientation is (bottom,right)
  // ABC: origin is centered and orientation is (up,right)
  // Following values are in cm, hence the 0.1 multiplier
  const float hoffset_cm = 0.1 * ((imgWidth*0.5) - hoffset_pix) * pix2mm;
  const float voffset_cm = -0.1 * ((imgHeight*0.5) - voffset_pix) * pix2mm; // vertical flip
  const float haperture_cm = 0.1 * imgWidth * pix2mm;
  const float vaperture_cm = 0.1 * imgHeight * pix2mm;

  camSample.setFocalLength(focalLength_mm);
  camSample.setHorizontalAperture(haperture_cm);
  camSample.setVerticalAperture(vaperture_cm);
  camSample.setHorizontalFilmOffset(hoffset_cm);
  camSample.setVerticalFilmOffset(voffset_cm);
  
  // Add sensor width (largest image side) in pixels as custom property
  OUInt32Property propSensorWidth_pix(userProps, "mvg_sensorWidth_pix");
  propSensorWidth_pix.set(sensorWidth_pix);

  // Add image path as custom property
  if(!imagePath.empty())
  {
    // Set camera image plane 
    OStringProperty imagePlane(userProps, "mvg_imagePath");
    imagePlane.set(imagePath.c_str());
  }

  OUInt32Property propViewId(userProps, "mvg_viewId");
  propViewId.set(id_view);

  OUInt32Property propIntrinsicId(userProps, "mvg_intrinsicId");
  propIntrinsicId.set(id_intrinsic);
  
  OStringProperty mvg_intrinsicType(userProps, "mvg_intrinsicType");
  mvg_intrinsicType.set(cam->getTypeStr());
  
  std::vector<double> intrinsicParams = cam->getParams();
  ODoubleArrayProperty mvg_intrinsicParams(userProps, "mvg_intrinsicParams");
  mvg_intrinsicParams.set(intrinsicParams);
  
  camObj.getSchema().set(camSample);
}

void AlembicExporter::initAnimatedCamera(const std::string& cameraName)
{
  // Sample the time in order to have one keyframe every frame
  // nb: it HAS TO be attached to EACH keyframed properties
  TimeSamplingPtr tsp( new TimeSampling(1.0 / 24.0, 1.0 / 24.0) );
  
  // Create the camera transform object
  std::stringstream ss;
  ss << cameraName;
  mxform = Alembic::AbcGeom::OXform(mvgCameras, "camxform_" + ss.str());
  mxform.getSchema().setTimeSampling(tsp);
  
  // Create the camera parameters object (intrinsics & custom properties)
  mcamObj = OCamera(mxform, "camera_" + ss.str());
  mcamObj.getSchema().setTimeSampling(tsp);
  
  // Add the custom properties
  auto userProps = mcamObj.getSchema().getUserProperties();
  // Sensor width
  mpropSensorWidth_pix = OUInt32Property(userProps, "mvg_sensorWidth_pix", tsp);
  // Image path
  mimagePlane = OStringProperty(userProps, "mvg_imagePath", tsp);
  // View id
  mpropViewId = OUInt32Property(userProps, "mvg_viewId", tsp);
  // Intrinsic id
  mpropIntrinsicId = OUInt32Property(userProps, "mvg_intrinsicId", tsp);
  // Intrinsic type (ex: PINHOLE_CAMERA_RADIAL3)
  mmvg_intrinsicType = OStringProperty(userProps, "mvg_intrinsicType", tsp);
  // Intrinsic parameters
  mmvg_intrinsicParams = ODoubleArrayProperty(userProps, "mvg_intrinsicParams", tsp);
}

void AlembicExporter::addCameraKeyframe(const geometry::Pose3 &pose,
                                          const cameras::Pinhole_Intrinsic *cam,
                                          const std::string &imagePath,
                                          const IndexT id_view,
                                          const IndexT id_intrinsic,
                                          const float sensorWidth_mm)
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

  // Create the XformSample
  XformSample xformsample;
  xformsample.setMatrix(xformMatrix);
  
  // Attach it to the schema of the OXform
  mxform.getSchema().set(xformsample);
  
  // Camera intrinsic parameters
  CameraSample camSample;

  // Take the max of the image size to handle the case where the image is in portrait mode 
  const float imgWidth = cam->w();
  const float imgHeight = cam->h();
  const float sensorWidth_pix = std::max(imgWidth, imgHeight);
  const float sensorHeight_pix = std::min(imgWidth, imgHeight);
  const float imgRatio = sensorHeight_pix / sensorWidth_pix;
  const float focalLength_pix = cam->focal();
  const float hoffset_pix = cam->principal_point()(0);
  const float voffset_pix = cam->principal_point()(1);
  

  const float sensorHeight_mm = sensorWidth_mm * imgRatio;
  const float focalLength_mm = sensorWidth_mm * focalLength_pix / sensorWidth_pix;
  const float pix2mm = sensorWidth_mm / sensorWidth_pix;

  // openMVG: origin is (top,left) corner and orientation is (bottom,right)
  // ABC: origin is centered and orientation is (up,right)
  // Following values are in cm, hence the 0.1 multiplier
  const float hoffset_cm = 0.1 * ((imgWidth*0.5) - hoffset_pix) * pix2mm;
  const float voffset_cm = -0.1 * ((imgHeight*0.5) - voffset_pix) * pix2mm; // vertical flip
  const float haperture_cm = 0.1 * imgWidth * pix2mm;
  const float vaperture_cm = 0.1 * imgHeight * pix2mm;

  camSample.setFocalLength(focalLength_mm);
  camSample.setHorizontalAperture(haperture_cm);
  camSample.setVerticalAperture(vaperture_cm);
  camSample.setHorizontalFilmOffset(hoffset_cm);
  camSample.setVerticalFilmOffset(voffset_cm);
  
  // Add sensor width (largest image side) in pixels as custom property
  mpropSensorWidth_pix.set(sensorWidth_pix);
  
  // Set custom attributes
  // Image path
  if(!imagePath.empty())
  {
    mimagePlane.set(imagePath.c_str());
  }
  // View id
  mpropViewId.set(id_view);
  // Intrinsic id
  mpropIntrinsicId.set(id_intrinsic);
  // Intrinsic type
  mmvg_intrinsicType.set(cam->getTypeStr());
  // Intrinsic parameters
  std::vector<double> intrinsicParams = cam->getParams();
  mmvg_intrinsicParams.set(intrinsicParams);
  
  // Attach intrinsic parameters to camera object
  mcamObj.getSchema().set(camSample);
}

void AlembicExporter::add(const sfm::SfM_Data &sfmdata, sfm::ESfM_Data flags_part)
{

  if(flags_part & sfm::ESfM_Data::VIEWS || flags_part & sfm::ESfM_Data::EXTRINSICS)
  {
    for(const auto it : sfmdata.GetViews())
    {
      const sfm::View * view = it.second.get();
      openMVG::geometry::Pose3 pose;
      std::shared_ptr<cameras::IntrinsicBase> cam = std::make_shared<cameras::Pinhole_Intrinsic>();
      if(sfmdata.IsPoseAndIntrinsicDefined(view))
      {
        // OpenMVG Camera
        pose = sfmdata.GetPoseOrDie(view);
        auto iterIntrinsic = sfmdata.GetIntrinsics().find(view->id_intrinsic);
        cam = iterIntrinsic->second;
      }
      else if(!(flags_part & sfm::ESfM_Data::VIEWS))
      {
        // If we don't export views, skip cameras without valid pose.
        continue;
      }
      const std::string cameraName = stlplus::basename_part(view->s_Img_path);
      const std::string sView_filename = stlplus::create_filespec(sfmdata.s_root_path, view->s_Img_path);
      
      // Use a common sensor width if we don't have this information.
      // We chose a full frame 24x36 camera
      float sensorWidth_mm = 36.0;
      const sfm::View_Metadata* viewMetadata = dynamic_cast<const sfm::View_Metadata*>(view);
      if(viewMetadata)
      {
        sensorWidth_mm = std::stof(viewMetadata->metadata.at(std::string("sensor_width")));
      }
      appendCamera(cameraName, pose, dynamic_cast<openMVG::cameras::Pinhole_Intrinsic*>(cam.get()), sView_filename, view->id_view, view->id_intrinsic, sensorWidth_mm);
    }
  }
  if(flags_part & sfm::ESfM_Data::STRUCTURE)
  {
    addPoints(sfmdata.GetLandmarks(), (flags_part & sfm::ESfM_Data::OBSERVATIONS));
  }
}


} //namespace dataio
} //namespace openMVG

#endif //HAVE_ALEMBIC
