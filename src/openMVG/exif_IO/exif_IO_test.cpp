#ifdef USE_EXIV2
#include "exif_IO_Exiv2.hpp"
#endif
#include "exif_IO_openExif.hpp"

#include "testing/testing.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <iostream>
#include <memory>

using namespace std;
using namespace openMVG;

const std::string sImg = 
  stlplus::folder_part(
  stlplus::folder_part(
  stlplus::folder_up(string(THIS_SOURCE_DIR))))
    + "/openMVG_Samples/imageData/Exif_Test/100_7100.JPG";

#ifdef USE_EXIV2
TEST(Matching, Exif_IO_Exiv2_ReadDatas)
{
  std::auto_ptr<Exif_IO> exif_io (new Exif_IO_Exiv2( sImg ) );

  EXPECT_TRUE( exif_io->doesHaveExifInfo());

  EXPECT_EQ( "EASTMAN KODAK COMPANY", exif_io->getBrand());
  EXPECT_EQ( "KODAK Z612 ZOOM DIGITAL CAMERA", exif_io->getModel());

  EXPECT_EQ( 2832, exif_io->getWidth());
  EXPECT_EQ( 2128, exif_io->getHeight());
  EXPECT_NEAR( 5.85, exif_io->getFocal(), 1e-2);

  EXPECT_EQ( "", exif_io->getLensModel());
}
#endif 


TEST(Matching, Exif_IO_openExif_ReadDatas)
{
  std::auto_ptr<Exif_IO> exif_io ( new Exif_IO_OpenExif( sImg ) );

  EXPECT_TRUE( exif_io->doesHaveExifInfo());

  EXPECT_EQ( "EASTMAN KODAK COMPANY", exif_io->getBrand());
  EXPECT_EQ( "KODAK Z612 ZOOM DIGITAL CAMERA", exif_io->getModel());

  EXPECT_EQ( 2832, exif_io->getWidth());
  EXPECT_EQ( 2128, exif_io->getHeight());
  EXPECT_NEAR( 5.85, exif_io->getFocal(), 1e-2);

  EXPECT_EQ( "", exif_io->getLensModel());
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */

