// Copyright (c) 2015 cpichard.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
        
#pragma once

#if HAVE_ALEMBIC

#include <Alembic/AbcGeom/All.h>
#include <Alembic/AbcCoreHDF5/All.h>

using namespace Alembic::Abc;
namespace AbcG = Alembic::AbcGeom;
using namespace AbcG;

namespace openMVG {
namespace dataio {

class AlembicImport 
{
public:
  explicit AlembicImport(const char* fileName);
  ~AlembicImport() = default;

  void populate(/*@todo put something here*/);

private:
  //@todo complete the interface, also maybe parameters need to be passed by reference?
  void visitObject(IObject iObj, M44d mat);
  IObject _rootEntity;
};

} // namespace mockup
} // namespace openMVG

#endif // WITH_ALEMBIC

