// Copyright (c) 2013 Pierre Moulon, Bruno Duisit.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

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

  EXPECT_TRUE( parseDatabase( sfileDatabase, vec_database ) );
  EXPECT_TRUE( getInfo( sModel, vec_database, datasheet ) );
  EXPECT_EQ( "Canon PowerShot SD900", datasheet.model_ );
  EXPECT_EQ( 7.11, datasheet.sensorSize_ );
}

TEST(Matching, ParseDatabaseA710_IS)
{
  std::vector<Datasheet> vec_database;
  Datasheet datasheet;
  const std::string sfileDatabase = stlplus::create_filespec( std::string(THIS_SOURCE_DIR), sDatabase );
  const std::string sModel = "Canon PowerShot A710 IS";

  EXPECT_TRUE( parseDatabase( sfileDatabase, vec_database ) );
  EXPECT_TRUE( getInfo( sModel, vec_database, datasheet ) );
  EXPECT_EQ( "Canon PowerShot A710 IS", datasheet.model_ );
  EXPECT_EQ( 5.75, datasheet.sensorSize_ );
}

TEST(Matching, ParseDatabaseNotExist)
{
  std::vector<Datasheet> vec_database;
  Datasheet datasheet;
  const std::string sfileDatabase = stlplus::create_filespec( std::string(THIS_SOURCE_DIR), sDatabase );
  const std::string sModel = "NotExistModel";

  EXPECT_TRUE( parseDatabase( sfileDatabase, vec_database ) );
  EXPECT_FALSE( getInfo( sModel, vec_database, datasheet ) );
}

TEST(Matching, ParseDatabaseCanon_EOS_550D)
{
  std::vector<Datasheet> vec_database;
  Datasheet datasheet;
  const std::string sfileDatabase = stlplus::create_filespec( std::string(THIS_SOURCE_DIR), sDatabase );
  const std::string sModel = "Canon EOS 550D";

  EXPECT_TRUE( parseDatabase( sfileDatabase, vec_database ) );
  EXPECT_TRUE( getInfo( sModel, vec_database, datasheet ) );
  EXPECT_EQ( 22.3, datasheet.sensorSize_ );
}

TEST(Matching, ParseDatabaseCanon_EOS_5D_Mark_II)
{
  std::vector<Datasheet> vec_database;
  Datasheet datasheet;
  const std::string sfileDatabase = stlplus::create_filespec( std::string(THIS_SOURCE_DIR), sDatabase );
  const std::string sModel = "Canon EOS 5D Mark II";

  EXPECT_TRUE( parseDatabase( sfileDatabase, vec_database ) );
  EXPECT_TRUE( getInfo( sModel, vec_database, datasheet ) );
  EXPECT_EQ( 36, datasheet.sensorSize_ );
}

TEST(Matching, ParseDatabaseCanon_EOS_1100D)
{
  std::vector<Datasheet> vec_database;
  Datasheet datasheet;
  const std::string sfileDatabase = stlplus::create_filespec( std::string(THIS_SOURCE_DIR), sDatabase );
  const std::string sModel = "Canon EOS 1100D";

  EXPECT_TRUE( parseDatabase( sfileDatabase, vec_database ) );
  EXPECT_TRUE( getInfo( sModel, vec_database, datasheet ) );
  EXPECT_EQ( 22.2, datasheet.sensorSize_ );
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
