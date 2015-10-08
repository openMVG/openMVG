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
    Pose3 pose(cam_r, cam_t);
    
    ICamera camera(iObj, kWrapExisting);
    ICameraSchema cs = camera.getSchema();
    CameraSample camSample = cs.getValue();

    double width = camSample.getHorizontalAperture(); // TODO
    double height = camSample.getVerticalAperture();

    Pinhole_Intrinsic intrinsic(
      width, height,
      camSample.getFocalLength(), // TODO: pix
      camSample.getHorizontalFilmOffset(), camSample.getVerticalFilmOffset());

    IndexT id_view = sfmdata.GetViews().size();
    IndexT id_pose = sfmdata.GetPoses().size();
    IndexT id_intrinsics = sfmdata.GetIntrinsics().size();
    
    const std::string objectName = iObj.getName();

    sfmdata.views.emplace(id_view, std::make_shared<View>(objectName, id_view, id_pose, id_intrinsics, width, height));
    sfmdata.poses.emplace(id_pose, pose);
    sfmdata.intrinsics.emplace(id_intrinsics, std::make_shared<Pinhole_Intrinsic>(intrinsic));
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

