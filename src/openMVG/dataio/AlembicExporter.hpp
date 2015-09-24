#if HAVE_ALEMBIC

#pragma once


#include <Alembic/AbcGeom/All.h>
#include <Alembic/AbcCoreHDF5/All.h>

#include <openMVG/sfm/sfm_data.hpp>

namespace openMVG {
namespace dataio {

class AlembicExporter
{
public:

  AlembicExporter(const std::string &filename)
  : archive(Alembic::AbcCoreHDF5::WriteArchive(), filename)
  , topObj(archive, Alembic::Abc::kTop)
  , counter(0)
  {
  }

  void addPoints(const sfm::Landmarks &points);

  // Add a camera to the archive
  void appendCamera(const std::string &name, 
                    const geometry::Pose3 &pose, 
                    const cameras::Pinhole_Intrinsic *cam);
  
  void addSfmData(const sfm::SfM_Data &sfmdata);

  virtual ~AlembicExporter();

private:
  Alembic::Abc::OArchive archive;
  Alembic::Abc::OObject topObj;
  unsigned int counter;
};

} // namespace data_io
} // namespace openMVG

#endif // HAVE_ALEMBIC

