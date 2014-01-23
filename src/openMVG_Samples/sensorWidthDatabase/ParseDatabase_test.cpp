#include "testing/testing.h"
#include "ParseDatabase.hpp"

#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <string>

TEST(Matching, ParseDatabaseSD900)
{
  std::vector<Datasheet> vec_database;
  Datasheet datasheet;
  std::string sfileDatabase = stlplus::create_filespec( std::string(THIS_SOURCE_DIR), "cameraGenerated.txt" );
  std::string sModel = "Canon PowerShot SD900";
  std::string sBrand = "Canon";

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
  std::string sfileDatabase = stlplus::create_filespec( std::string(THIS_SOURCE_DIR), "cameraGenerated.txt" );
  std::string sModel = "Canon PowerShot A710 IS";
  std::string sBrand = "Canon";

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
  std::string sfileDatabase = stlplus::create_filespec( std::string(THIS_SOURCE_DIR), "cameraGenerated.txt" );
  std::string sModel = "NotExistModel";
  std::string sBrand = "NotExistBrand";

  EXPECT_TRUE( parseDatabase( sfileDatabase, vec_database ) );
  EXPECT_FALSE( getInfo( sBrand, sModel, vec_database, datasheet ) );
}


TEST(Matching, ParseDatabaseCanon_EOS_550D)
{
  std::vector<Datasheet> vec_database;
  Datasheet datasheet;
  std::string sfileDatabase = stlplus::create_filespec( std::string(THIS_SOURCE_DIR), "cameraGenerated.txt" );
  std::string sModel = "Canon EOS 550D";
  std::string sBrand = "Canon";

  EXPECT_TRUE( parseDatabase( sfileDatabase, vec_database ) );
  EXPECT_TRUE( getInfo( sBrand, sModel, vec_database, datasheet ) );
  EXPECT_EQ( 22.3, datasheet._sensorSize );
}

TEST(Matching, ParseDatabaseCanon_EOS_5D_Mark_II)
{
  std::vector<Datasheet> vec_database;
  Datasheet datasheet;
  std::string sfileDatabase = stlplus::create_filespec( std::string(THIS_SOURCE_DIR), "cameraGenerated.txt" );
  std::string sModel = "Canon EOS 5D Mark II";
  std::string sBrand = "Canon";

  EXPECT_TRUE( parseDatabase( sfileDatabase, vec_database ) );
  EXPECT_TRUE( getInfo( sBrand, sModel, vec_database, datasheet ) );
  EXPECT_EQ( 36, datasheet._sensorSize );
}

TEST(Matching, ParseDatabaseCanon_EOS_1100D)
{
  std::vector<Datasheet> vec_database;
  Datasheet datasheet;
  std::string sfileDatabase = stlplus::create_filespec( std::string(THIS_SOURCE_DIR), "cameraGenerated.txt" );
  std::string sModel = "Canon EOS 1100D";
  std::string sBrand = "Canon";

  EXPECT_TRUE( parseDatabase( sfileDatabase, vec_database ) );
  EXPECT_TRUE( getInfo( sBrand, sModel, vec_database, datasheet ) );
  EXPECT_EQ( 22.2, datasheet._sensorSize );
}


/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
