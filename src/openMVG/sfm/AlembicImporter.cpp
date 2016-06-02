// Copyright (c) 2015 cpichard.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#if HAVE_ALEMBIC

#include "AlembicImporter.hpp"

#include <Alembic/AbcGeom/All.h>
#include <Alembic/AbcCoreFactory/All.h>


using namespace Alembic::Abc;
namespace AbcG = Alembic::AbcGeom;
using namespace AbcG;

namespace openMVG {
namespace dataio {

  
/**
 * @brief Retrieve an Abc property.
 *         Maya convert everything into arrays. So here we made a trick
 *         to retrieve the element directly or the first element
 *         if it's an array.
 * @param userProps
 * @param id
 * @param sampleTime
 * @return value
 */
template<class AbcProperty>
typename AbcProperty::traits_type::value_type getAbcProp(ICompoundProperty& userProps, const Alembic::Abc::PropertyHeader& propHeader, const std::string& id, chrono_t sampleTime)
{
  typedef typename AbcProperty::traits_type traits_type;
  typedef typename traits_type::value_type value_type;
  typedef typename Alembic::Abc::ITypedArrayProperty<traits_type> array_type;
  typedef typename array_type::sample_ptr_type array_sample_ptr_type;

  // Maya transforms everything into arrays
  if(propHeader.isArray())
  {
    Alembic::Abc::ITypedArrayProperty<traits_type> prop(userProps, id);
    array_sample_ptr_type sample;
    prop.get(sample, ISampleSelector(sampleTime));
    return (*sample)[0];
  }
  else
  {
    value_type v;
    AbcProperty prop(userProps, id);
    prop.get(v, ISampleSelector(sampleTime));
    return v;
  }
}


template<class ABCSCHEMA>
inline ICompoundProperty getAbcUserProperties(ABCSCHEMA& schema)
{
  ICompoundProperty userProps = schema.getUserProperties();
  if(userProps && userProps.getNumProperties() != 0)
    return userProps;

  // Maya always use ArbGeomParams instead of user properties.
  return schema.getArbGeomParams();
}


bool AlembicImporter::readPointCloud(IObject iObj, M44d mat, sfm::SfM_Data &sfmdata, sfm::ESfM_Data flags_part)
{

  using namespace openMVG::geometry;
  using namespace openMVG::sfm;

  IPoints points(iObj, kWrapExisting);
  IPointsSchema& ms = points.getSchema();
  P3fArraySamplePtr positions = ms.getValue().getPositions();

  ICompoundProperty userProps = getAbcUserProperties(ms);
  ICompoundProperty arbGeom = ms.getArbGeomParams();

  C3fArraySamplePtr sampleColors;
  if(arbGeom.getPropertyHeader("color"))
  {
    IC3fArrayProperty propColor(arbGeom, "color");
    propColor.get(sampleColors);
    if(sampleColors->size() != positions->size())
    {
      std::cerr << "[Alembic Importer] WARNING: colors will be ignored. Color vector size: " << sampleColors->size() << ", positions vector size: " << positions->size() << std::endl;
      sampleColors.reset();
    }
  }

  // Number of points before adding the Alembic data
  const std::size_t nbPointsInit = sfmdata.structure.size();
  for(std::size_t point3d_i = 0;
      point3d_i < positions->size();
      ++point3d_i)
  {
    const P3fArraySamplePtr::element_type::value_type & pos_i = positions->get()[point3d_i];
    Landmark& landmark = sfmdata.structure[nbPointsInit + point3d_i] = Landmark(Vec3(pos_i.x, pos_i.y, pos_i.z));
    if(sampleColors)
    {
      const P3fArraySamplePtr::element_type::value_type & color_i = sampleColors->get()[point3d_i];
      landmark.rgb = image::RGBColor(color_i[0], color_i[1], color_i[2]);
    }
  }

  if(userProps.getPropertyHeader("mvg_visibilitySize") &&
     userProps.getPropertyHeader("mvg_visibilityIds") &&
     userProps.getPropertyHeader("mvg_visibilityFeatPos")
     )
  {
    IUInt32ArrayProperty propVisibilitySize(userProps, "mvg_visibilitySize");
    UInt32ArraySamplePtr sampleVisibilitySize;
    propVisibilitySize.get(sampleVisibilitySize);

    IUInt32ArrayProperty propVisibilityIds(userProps, "mvg_visibilityIds");
    UInt32ArraySamplePtr sampleVisibilityIds;
    propVisibilityIds.get(sampleVisibilityIds);

    IFloatArrayProperty propFeatPos2d(userProps, "mvg_visibilityFeatPos");
    FloatArraySamplePtr sampleFeatPos2d;
    propFeatPos2d.get(sampleFeatPos2d);

    if( positions->size() != sampleVisibilitySize->size() )
    {
      std::cerr << "ABC Error: number of observations per 3D point should be identical to the number of 2D features." << std::endl;
      std::cerr << "Number of observations per 3D point size is " << sampleVisibilitySize->size() << std::endl;
      std::cerr << "Number of 3D points is " << positions->size() << std::endl;
      return false;
    }
    if( sampleVisibilityIds->size() != sampleFeatPos2d->size() )
    {
      std::cerr << "ABC Error: visibility Ids and features 2D pos should have the same size." << std::endl;
      std::cerr << "Visibility Ids size is " << sampleVisibilityIds->size() << std::endl;
      std::cerr << "Features 2d Pos size is " << sampleFeatPos2d->size() << std::endl;
      return false;
    }

    std::size_t obsGlobal_i = 0;
    for(std::size_t point3d_i = 0;
        point3d_i < positions->size();
        ++point3d_i)
    {
      Landmark& landmark = sfmdata.structure[nbPointsInit + point3d_i];
      // Number of observation for this 3d point
      const std::size_t visibilitySize = (*sampleVisibilitySize)[point3d_i];

      for(std::size_t obs_i = 0;
          obs_i < visibilitySize*2;
          obs_i+=2, obsGlobal_i+=2)
      {

        const int viewID = (*sampleVisibilityIds)[obsGlobal_i];
        const int featID = (*sampleVisibilityIds)[obsGlobal_i+1];
        Observation& obs = landmark.obs[viewID];
        obs.id_feat = featID;

        const float posX = (*sampleFeatPos2d)[obsGlobal_i];
        const float posY = (*sampleFeatPos2d)[obsGlobal_i+1];
        obs.x[0] = posX;
        obs.x[1] = posY;
      }
    }
  }
  return true;
}

bool AlembicImporter::readCamera(IObject iObj, M44d mat, sfm::SfM_Data &sfmdata, sfm::ESfM_Data flags_part, const chrono_t sampleTime)
{
  using namespace openMVG::geometry;
  using namespace openMVG::cameras;
  using namespace openMVG::sfm;

  ICamera camera(iObj, kWrapExisting);
  ICameraSchema cs = camera.getSchema();
  CameraSample camSample;
  if(sampleTime == 0)
    camSample = cs.getValue();
  else
    camSample = cs.getValue(ISampleSelector(sampleTime));

  // Check if we have an associated image plane
  ICompoundProperty userProps = getAbcUserProperties(cs);
  std::string imagePath;
  std::vector<int> sensorSize_pix = {2048, 2048};
  std::string mvg_intrinsicType = EINTRINSIC_enumToString(PINHOLE_CAMERA);
  std::vector<double> mvg_intrinsicParams;
  IndexT id_view = sfmdata.GetViews().size();
  IndexT id_intrinsic = sfmdata.GetIntrinsics().size();
  if(userProps)
  {
    if(const Alembic::Abc::PropertyHeader *propHeader = userProps.getPropertyHeader("mvg_imagePath"))
    {
      imagePath = getAbcProp<Alembic::Abc::IStringProperty>(userProps, *propHeader, "mvg_imagePath", sampleTime);
    }
    if(userProps.getPropertyHeader("mvg_sensorSizePix"))
    {
      Alembic::Abc::IUInt32ArrayProperty prop(userProps, "mvg_sensorSizePix");
      std::shared_ptr<UInt32ArraySample> sample;
      prop.get(sample, ISampleSelector(sampleTime));
      sensorSize_pix.assign(sample->get(), sample->get()+sample->size());
      assert(sensorSize_pix.size() == 2);
    }
    if(const Alembic::Abc::PropertyHeader *propHeader = userProps.getPropertyHeader("mvg_intrinsicType"))
    {
      mvg_intrinsicType = getAbcProp<Alembic::Abc::IStringProperty>(userProps, *propHeader, "mvg_intrinsicType", sampleTime);
    }
    if(userProps.getPropertyHeader("mvg_intrinsicParams"))
    {
      Alembic::Abc::IDoubleArrayProperty prop(userProps, "mvg_intrinsicParams");
      std::shared_ptr<DoubleArraySample> sample;
      prop.get(sample, ISampleSelector(sampleTime));
      mvg_intrinsicParams.assign(sample->get(), sample->get()+sample->size());
    }
    if(const Alembic::Abc::PropertyHeader *propHeader = userProps.getPropertyHeader("mvg_viewId"))
    {
      try {
        id_view = getAbcProp<Alembic::Abc::IUInt32Property>(userProps, *propHeader, "mvg_viewId", sampleTime);
      } catch(Alembic::Util::v7::Exception&)
      {
        id_view = getAbcProp<Alembic::Abc::IInt32Property>(userProps, *propHeader, "mvg_viewId", sampleTime);
      }
    }
    if(const Alembic::Abc::PropertyHeader *propHeader = userProps.getPropertyHeader("mvg_intrinsicId"))
    {
      try {
        id_intrinsic = getAbcProp<Alembic::Abc::IUInt32Property>(userProps, *propHeader, "mvg_intrinsicId", sampleTime);
      } catch(Alembic::Util::v7::Exception&)
      {
        id_intrinsic = getAbcProp<Alembic::Abc::IInt32Property>(userProps, *propHeader, "mvg_intrinsicId", sampleTime);
      }
    }
  }

  // OpenMVG Camera
  Mat3 cam_r;
  cam_r(0,0) = mat[0][0];
  cam_r(0,1) = mat[0][1];
  cam_r(0,2) = mat[0][2];
  cam_r(1,0) = mat[1][0];
  cam_r(1,1) = mat[1][1];
  cam_r(1,2) = mat[1][2];
  cam_r(2,0) = mat[2][0];
  cam_r(2,1) = mat[2][1];
  cam_r(2,2) = mat[2][2];
  Vec3 cam_t;
  cam_t(0) = mat[3][0];
  cam_t(1) = mat[3][1];
  cam_t(2) = mat[3][2];

  // Correct camera orientation from alembic
  const Mat3 scale = Vec3(1,-1,-1).asDiagonal();
  cam_r = scale*cam_r;

  Pose3 pose(cam_r, cam_t);

  // Get known values from alembic
  const float haperture_cm = camSample.getHorizontalAperture();
  const float vaperture_cm = camSample.getVerticalAperture();

  // Compute other needed values
  const float sensorWidth_mm = std::max(vaperture_cm, haperture_cm) * 10.0;
  const float mm2pix = sensorSize_pix.at(0) / sensorWidth_mm;
  const float imgWidth = haperture_cm * 10.0 * mm2pix;
  const float imgHeight = vaperture_cm * 10.0 * mm2pix;

  // Create intrinsic parameters object
  std::shared_ptr<Pinhole_Intrinsic> pinholeIntrinsic = createPinholeIntrinsic(EINTRINSIC_stringToEnum(mvg_intrinsicType));
  pinholeIntrinsic->setWidth(imgWidth);
  pinholeIntrinsic->setHeight(imgHeight);
  pinholeIntrinsic->updateFromParams(mvg_intrinsicParams);

  // Add imported data to the SfM_Data container TODO use UID
  sfmdata.views.emplace(id_view, std::make_shared<View>(imagePath, id_view, id_intrinsic, id_view, imgWidth, imgHeight));
  sfmdata.poses.emplace(id_view, pose);
  sfmdata.intrinsics.emplace(id_intrinsic, pinholeIntrinsic);

  return true;
}


// Top down read of 3d objects
void AlembicImporter::visitObject(IObject iObj, M44d mat, sfm::SfM_Data &sfmdata, sfm::ESfM_Data flags_part)
{
  // std::cout << "ABC visit: " << iObj.getFullName() << std::endl;
  
  const MetaData& md = iObj.getMetaData();
  if(IPoints::matches(md) && (flags_part & sfm::ESfM_Data::STRUCTURE))
  {
    readPointCloud(iObj, mat, sfmdata, flags_part);
  }
  else if(IXform::matches(md))
  {
    IXform xform(iObj, kWrapExisting);
    XformSample xs;
    if(xform.getSchema().getNumSamples() == 1)
    {
      xform.getSchema().get(xs);
      mat *= xs.getMatrix();
    }
    // If we have an animated camera we handle it with the xform here
    else
    {
      std::cout << xform.getSchema().getNumSamples() << " samples found in this animated xform." << std::endl;
      const float timebase = 1.0 / 24.0;
      const float timestep = 1.0 / 24.0;
      for(int frame = 0; frame < xform.getSchema().getNumSamples(); ++frame)
      {
        xform.getSchema().get(xs, ISampleSelector(frame * timestep + timebase));
        readCamera(iObj.getChild(0), mat * xs.getMatrix(), sfmdata, flags_part, frame * timestep + timebase);
      }
    }
  }
  else if(ICamera::matches(md) && (flags_part & sfm::ESfM_Data::EXTRINSICS))
  {
    ICamera check_cam(iObj, kWrapExisting);
    // If it's not an animated camera we add it here
    if(check_cam.getSchema().getNumSamples() == 1)
    {
      readCamera(iObj, mat, sfmdata, flags_part);
    }
  }

  // Recurse
  for(size_t i = 0; i < iObj.getNumChildren(); i++)
  {
    visitObject(iObj.getChild(i), mat, sfmdata, flags_part);
  }
}

AlembicImporter::AlembicImporter(const std::string &filename)
{
  Alembic::AbcCoreFactory::IFactory factory;
  Alembic::AbcCoreFactory::IFactory::CoreType coreType;
  Abc::IArchive archive = factory.getArchive(filename, coreType);

  // TODO : test if archive is correctly opened
  _rootEntity = archive.getTop();
}

void AlembicImporter::populate(sfm::SfM_Data &sfmdata, sfm::ESfM_Data flags_part)
{
  // TODO : handle the case where the archive wasn't correctly opened
  M44d xformMat;
  visitObject(_rootEntity, xformMat, sfmdata, flags_part);

  // TODO: fusion of common intrinsics
}

} // namespace data_io
} // namespace openMVG

#endif // WITH_ALEMBIC

