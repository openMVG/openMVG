#if HAVE_ALEMBIC

#pragma once


#include <Alembic/AbcGeom/All.h>
#include <Alembic/AbcCoreHDF5/All.h>

#include <openMVG/sfm/sfm_data.hpp>
#include <openMVG/sfm/sfm_data_io.hpp>

namespace openMVG {
namespace dataio {

class AlembicExporter
{
public:

  AlembicExporter(const std::string &filename);

  /**
   * @brief Add a set of 3D points from a SFM scene
   * @param points The 3D points to add
   */
  void addPoints(const sfm::Landmarks &points,
                 bool withVisibility=true);

  /**
   * @brief Add a single camera
   * 
   * @param[in] name An identifier for the camera
   * @param[in] pose The camera pose
   * @param[in] cam The camera intrinsics
   * @param[in] imagePath Path to the image
   * @param[in] id_view View id
   * @param[in] id_intrinsic Intrinsic id
   * @param[in] sensorWidth_mm Width of the sensor in millimeters
   * @param[in] id_pose Pose id
   */
  void appendCamera(const std::string &name, 
                    const geometry::Pose3 &pose, 
                    const cameras::Pinhole_Intrinsic *cam,
                    const std::string &imagePath,
                    const IndexT id_view,
                    const IndexT id_intrinsic,
                    const float sensorWidth_mm=36.0,
                    const IndexT id_pose=UndefinedIndexT);
  
  /**
   * @brief Initiate an animated camera
   * 
   * @param cameraName An identifier for the camera
   */
  void initAnimatedCamera(const std::string &cameraName);
  
  /**
   * @brief Add a keyframe to the animated camera
   * 
   * @param pose The camera pose
   * @param cam The camera intrinsics parameters
   * @param imagePath The localized image path
   * @param id_view View id
   * @param id_intrinsic Intrinsic id
   * @param sensorWidth_mm Width of the sensor in millimeters
   */
  void addCameraKeyframe(const geometry::Pose3 &pose,
                           const cameras::Pinhole_Intrinsic *cam,
                           const std::string &imagePath,
                           const IndexT id_view,
                           const IndexT id_intrinsic,
                           const float sensorWidth_mm=36.0);
  
  /**
   * @brief Register keyframe on the previous values
   */
  void jumpKeyframe(const std::string &imagePath = std::string());
  
  /**
   * @brief Add SfM Data
   * 
   * @param sfmdata SfM_Data container
   * @param flags_part filter the elements to add
   */
  void add(const sfm::SfM_Data &sfmdata, sfm::ESfM_Data flags_part = sfm::ESfM_Data::ALL);

  virtual ~AlembicExporter();

private:
  Alembic::Abc::OArchive archive;
  Alembic::Abc::OObject topObj;
  Alembic::Abc::OObject mvgRoot;
  Alembic::Abc::OObject mvgCameras;
  Alembic::Abc::OObject mvgCloud;
  Alembic::Abc::OObject mvgPointCloud;
  unsigned int counter;
  Alembic::AbcGeom::OXform mxform;
  Alembic::AbcGeom::OCamera mcamObj;
  Alembic::AbcGeom::OUInt32ArrayProperty mpropSensorSize_pix;
  Alembic::AbcGeom::OStringProperty mimagePlane;
  Alembic::AbcGeom::OUInt32Property mpropViewId;
  Alembic::AbcGeom::OUInt32Property mpropIntrinsicId;
  Alembic::AbcGeom::OStringProperty mmvg_intrinsicType;
  Alembic::AbcGeom::ODoubleArrayProperty mmvg_intrinsicParams;
};

} // namespace data_io
} // namespace openMVG

#endif // HAVE_ALEMBIC

