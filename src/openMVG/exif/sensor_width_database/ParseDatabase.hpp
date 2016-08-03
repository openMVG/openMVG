#pragma once

#include "openMVG/stl/split.hpp"
#include "datasheet.hpp"

#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <iterator>

#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

namespace openMVG {
namespace exif {
namespace sensordb {

// Parse the database
inline bool parseDatabase( const std::string& sfileDatabase, std::vector<Datasheet>& vec_database )
{
  std::ifstream iFilein( sfileDatabase.c_str() );
  if(!iFilein || !stlplus::is_file(sfileDatabase))
    return false;

  std::string line;
  while ( iFilein.good() )
  {
    getline( iFilein, line );
    if ( !line.empty() )
    {
      //std::stringstream sStream( line );
      if ( line[0] != '#' )
      {
        std::vector<std::string> values;
        stl::split(line, ";", values);
        if ( values.size() == 3 )
        {
          const std::string brand = values[0];
          const std::string model = values[1];
          const double sensorSize = atof( values[2].c_str() );
          vec_database.push_back( Datasheet( brand, model, sensorSize ) );
        }
      }
    }
  }
  return true;
}

// Get information for the given camera model
inline bool getInfo( const std::string& sBrand, const std::string& sModel, const std::vector<Datasheet>& vec_database, Datasheet& datasheetContent )
{
  bool existInDatabase = false;

  Datasheet refDatasheet( sBrand, sModel, -1. );
  std::vector<Datasheet>::const_iterator datasheet = std::find( vec_database.begin(), vec_database.end(), refDatasheet );
  if ( datasheet != vec_database.end() )
  {
    datasheetContent = *datasheet;
    existInDatabase = true;
  }

  return existInDatabase;
}

}
}
}
