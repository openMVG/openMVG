#ifndef PARSE_DATABASE_HPP
#define PARSE_DATABASE_HPP

#include "openMVG/split/split.hpp"
#include "datasheet.hpp"

#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <iterator>

// Parse the database
bool parseDatabase( const std::string& sfileDatabase, std::vector<Datasheet>& vec_database )
{
  bool createDatabase = false;
  std::ifstream iFilein( sfileDatabase.c_str() );
  if ( iFilein.is_open() )
  {
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
          split( line, ";", values );
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
    createDatabase = true;
  }
  else
  {
    std::cerr<< "Cannot open the database file: "
      << sfileDatabase << std::endl;
  }

  return createDatabase;
}

// Get information for the given camera model
bool getInfo( const std::string& sBrand, const std::string& sModel, const std::vector<Datasheet>& vec_database, Datasheet& datasheetContent )
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

#endif // PARSE_DATABASE_HPP
