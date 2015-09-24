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

void AlembicImport::visitObject(IObject iObj, M44d mat)
{
  const MetaData& md = iObj.getMetaData();
  if(IPoints::matches(md))
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
  else if(ICamera::matches(md))
  {
    ICamera camera(iObj, kWrapExisting);
    ICameraSchema cs = camera.getSchema();
    CameraSample matrix = cs.getValue();
    //@todo put the camera somewhere...
  }

  // Recurse
  for(size_t i = 0; i < iObj.getNumChildren(); i++)
  {
    visitObject(iObj.getChild(i), mat);
  }
}

AlembicImport::AlembicImport(const char* filename)
{
  Alembic::AbcCoreFactory::IFactory factory;
  Alembic::AbcCoreFactory::IFactory::CoreType coreType;
  Abc::IArchive archive = factory.getArchive(filename, coreType);

  // TODO : test if archive is correctly opened
  _rootEntity = archive.getTop();
}

void AlembicImport::populate()
{
  // TODO : handle the case where the archive wasn't correctly opened
  M44d xformMat;
  visitObject(_rootEntity, xformMat);
}

} // namespace data_io
} // namespace openMVG
#endif // WITH_ALEMBIC

