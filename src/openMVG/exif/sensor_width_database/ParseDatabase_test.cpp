#include "testing/testing.h"
#include "ParseDatabase.hpp"

#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <string>

static const std::string sDatabase = "sensor_width_camera_database.txt";
TEST(Matching, InvalidDatabase)
{
  std::vector<Datasheet> vec_database;
  const std::string sfileDatabase = stlplus::create_filespec( std::string(THIS_SOURCE_DIR), "" );

  EXPECT_FALSE( parseDatabase( sfileDatabase, vec_database ) );
  EXPECT_TRUE( vec_database.empty() );
}

TEST(Matching, ValidDatabase)
{
  std::vector<Datasheet> vec_database;
  const std::string sfileDatabase = stlplus::create_filespec( std::string(THIS_SOURCE_DIR), sDatabase );

  EXPECT_TRUE( parseDatabase( sfileDatabase, vec_database ) );
  EXPECT_TRUE( !vec_database.empty() );

}

TEST(Matching, ParseDatabaseSD900)
{
  std::vector<Datasheet> vec_database;
  Datasheet datasheet;
  const std::string sfileDatabase = stlplus::create_filespec( std::string(THIS_SOURCE_DIR), sDatabase );
  const std::string sModel = "Canon PowerShot SD900";
  const std::string sBrand = "Canon";

  EXPECT_TRUE( parseDatabase( sfileDatabase, vec_database ) );
  EXPECT_TRUE( getInfo( sBrand, sModel, vec_database, datasheet ) );
  EXPECT_EQ( "Canon", datasheet._brand );
  EXPECT_EQ( "Canon PowerShot SD900", datasheet._model );
  EXPECT_EQ( 7.11, datasheet._sensorSize );
}

TEST(Matching, ParseDatabaseA710_IS)
{
  std::vector<Datasheet> vec_database;
  Datasheet datasheet;
  const std::string sfileDatabase = stlplus::create_filespec( std::string(THIS_SOURCE_DIR), sDatabase );
  const std::string sModel = "Canon PowerShot A710 IS";
  const std::string sBrand = "Canon";

  EXPECT_TRUE( parseDatabase( sfileDatabase, vec_database ) );
  EXPECT_TRUE( getInfo( sBrand, sModel, vec_database, datasheet ) );
  EXPECT_EQ( "Canon", datasheet._brand );
  EXPECT_EQ( "Canon PowerShot A710 IS", datasheet._model );
  EXPECT_EQ( 5.75, datasheet._sensorSize );
}

TEST(Matching, ParseDatabaseNotExist)
{
  std::vector<Datasheet> vec_database;
  Datasheet datasheet;
  const std::string sfileDatabase = stlplus::create_filespec( std::string(THIS_SOURCE_DIR), sDatabase );
  const std::string sModel = "NotExistModel";
  const std::string sBrand = "NotExistBrand";

  EXPECT_TRUE( parseDatabase( sfileDatabase, vec_database ) );
  EXPECT_FALSE( getInfo( sBrand, sModel, vec_database, datasheet ) );
}


TEST(Matching, ParseDatabaseCanon_EOS_550D)
{
  std::vector<Datasheet> vec_database;
  Datasheet datasheet;
  const std::string sfileDatabase = stlplus::create_filespec( std::string(THIS_SOURCE_DIR), sDatabase );
  const std::string sModel = "Canon EOS 550D";
  const std::string sBrand = "Canon";

  EXPECT_TRUE( parseDatabase( sfileDatabase, vec_database ) );
  EXPECT_TRUE( getInfo( sBrand, sModel, vec_database, datasheet ) );
  EXPECT_EQ( 22.3, datasheet._sensorSize );
}

TEST(Matching, ParseDatabaseCanon_EOS_5D_Mark_II)
{
  std::vector<Datasheet> vec_database;
  Datasheet datasheet;
  const std::string sfileDatabase = stlplus::create_filespec( std::string(THIS_SOURCE_DIR), sDatabase );
  const std::string sModel = "Canon EOS 5D Mark II";
  const std::string sBrand = "Canon";

  EXPECT_TRUE( parseDatabase( sfileDatabase, vec_database ) );
  EXPECT_TRUE( getInfo( sBrand, sModel, vec_database, datasheet ) );
  EXPECT_EQ( 36, datasheet._sensorSize );
}

TEST(Matching, ParseDatabaseCanon_EOS_1100D)
{
  std::vector<Datasheet> vec_database;
  Datasheet datasheet;
  const std::string sfileDatabase = stlplus::create_filespec( std::string(THIS_SOURCE_DIR), sDatabase );
  const std::string sModel = "Canon EOS 1100D";
  const std::string sBrand = "Canon";

  EXPECT_TRUE( parseDatabase( sfileDatabase, vec_database ) );
  EXPECT_TRUE( getInfo( sBrand, sModel, vec_database, datasheet ) );
  EXPECT_EQ( 22.2, datasheet._sensorSize );
}


/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
