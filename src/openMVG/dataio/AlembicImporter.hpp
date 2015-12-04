// Copyright (c) 2015 cpichard.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#if HAVE_ALEMBIC
        
#pragma once

#include <Alembic/AbcGeom/All.h>
#include <Alembic/AbcCoreHDF5/All.h>

#include <openMVG/sfm/sfm_data.hpp>
#include <openMVG/sfm/sfm_data_io.hpp>

using namespace Alembic::Abc;
namespace AbcG = Alembic::AbcGeom;
using namespace AbcG;

namespace openMVG {
namespace dataio {

class AlembicImporter 
{
public:
  explicit AlembicImporter(const std::string &filename);
  ~AlembicImporter() = default;

  void populate(sfm::SfM_Data &sfmdata, sfm::ESfM_Data flags_part = sfm::ESfM_Data::ALL);

private:
  bool readPointCloud(IObject iObj, M44d mat, sfm::SfM_Data &sfmdata, sfm::ESfM_Data flags_part);
  bool readCamera(IObject iObj, M44d mat, sfm::SfM_Data &sfmdata, sfm::ESfM_Data flags_part, const chrono_t sampleTime = 0.0);
  
  //@todo complete the interface, also maybe parameters need to be passed by reference?
  void visitObject(IObject iObj, M44d mat, sfm::SfM_Data &sfmdata, sfm::ESfM_Data flags_part = sfm::ESfM_Data::ALL);
  IObject _rootEntity;
};

} // namespace mockup
} // namespace openMVG

#endif // WITH_ALEMBIC

