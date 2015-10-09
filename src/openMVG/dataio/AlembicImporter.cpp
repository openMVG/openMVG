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

// Top down insertion of 3d objects

void AlembicImporter::visitObject(IObject iObj, M44d mat, sfm::SfM_Data &sfmdata, sfm::ESfM_Data flags_part)
{
  using namespace openMVG::geometry;
  using namespace openMVG::cameras;
  using namespace openMVG::sfm;

  const MetaData& md = iObj.getMetaData();
  if(IPoints::matches(md) && (flags_part & sfm::ESfM_Data::STRUCTURE))
  {
    IPoints points(iObj, kWrapExisting);
    IPointsSchema ms = points.getSchema();
    P3fArraySamplePtr positions = ms.getValue().getPositions();

    //@todo use positions->get() and positions->size());
    //@todo put the points somewhere...
  }
  else if(IXform::matches(md))
  {
    IXform xform(iObj, kWrapExisting);
    XformSample xs;
    xform.getSchema().get(xs);
    mat *= xs.getMatrix();
  }
  else if(ICamera::matches(md) && (flags_part & sfm::ESfM_Data::EXTRINSICS))
  {
    ICamera camera(iObj, kWrapExisting);
    ICameraSchema cs = camera.getSchema();
    CameraSample camSample = cs.getValue();
    
    // Check if we have an associated image plane
    ICompoundProperty userProps = cs.getUserProperties();
    std::string imagePath;
    float sensorWidth_pix = 2048.0;
    std::string mvg_intrinsicType = "Pinhole_Intrinsic";
    std::vector<double> mvg_intrinsicParams;
    if(userProps)
    {
        if(userProps.getPropertyHeader("imagePath"))
        {
          Alembic::Abc::IStringProperty prop(userProps, "imagePath");
          prop.get(imagePath);
        }
        if(userProps.getPropertyHeader("sensorWidth_pix"))
        {
          Alembic::Abc::IFloatProperty prop(userProps, "sensorWidth_pix");
          prop.get(sensorWidth_pix);
        }
        if(userProps.getPropertyHeader("mvg_intrinsicType"))
        {
          Alembic::Abc::IStringProperty prop(userProps, "mvg_intrinsicType");
          prop.get(mvg_intrinsicType);
        }
        if(userProps.getPropertyHeader("mvg_intrinsicParams"))
        {
          Alembic::Abc::IDoubleArrayProperty prop(userProps, "mvg_intrinsicParams");
          std::shared_ptr<DoubleArraySample> sample;
          prop.get(sample);
          mvg_intrinsicParams.assign(sample->get(), sample->get()+sample->size());
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
    Mat3 scale;
    scale(0,0) = 1;
    scale(1,1) = -1;
    scale(2,2) = -1;
    cam_r = scale*cam_r;
    
    Pose3 pose(cam_r, cam_t);

    // Get known values from alembic
    const float haperture_cm = camSample.getHorizontalAperture();
    const float vaperture_cm = camSample.getVerticalAperture();
    const float hoffset_cm = camSample.getHorizontalFilmOffset();
    const float voffset_cm = camSample.getVerticalFilmOffset();
    const float focalLength_mm = camSample.getFocalLength();
    
    // Compute other needed values
    const float sensorWidth_mm = std::max(vaperture_cm, haperture_cm) * 10.0;
    const float mm2pix = sensorWidth_pix / sensorWidth_mm;
    const float imgWidth = haperture_cm * 10.0 * mm2pix;
    const float imgHeight = vaperture_cm * 10.0 * mm2pix;
    const float focalLength_pix = focalLength_mm * mm2pix;
    
    // Following values are in cm, hence the 10.0 multiplier
    const float hoffset_pix = (imgWidth*0.5) - (10.0 * hoffset_cm * mm2pix);
    const float voffset_pix = (imgHeight*0.5) + (10.0 * voffset_cm * mm2pix);


    IndexT id_view = sfmdata.GetViews().size(); // TODO real index
    IndexT id_pose = sfmdata.GetPoses().size();
    IndexT id_intrinsic = sfmdata.GetIntrinsics().size();
    // Create intrinsic parameters object
    std::shared_ptr<Pinhole_Intrinsic> pinholeIntrinsic = createPinholeIntrinsic(EINTRINSIC_stringToEnum(mvg_intrinsicType));
    pinholeIntrinsic->setWidth(imgWidth);
    pinholeIntrinsic->setHeight(imgHeight);
    pinholeIntrinsic->setK(focalLength_pix, hoffset_pix, voffset_pix);
    pinholeIntrinsic->updateFromParams(mvg_intrinsicParams);

    // Add imported data to the SfM_Data container TODO use UID
    sfmdata.views.emplace(id_view, std::make_shared<View>(imagePath, id_view, id_pose, id_intrinsic, imgWidth, imgHeight));
    sfmdata.poses.emplace(id_pose, pose);
    sfmdata.intrinsics.emplace(id_intrinsic, pinholeIntrinsic);
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

