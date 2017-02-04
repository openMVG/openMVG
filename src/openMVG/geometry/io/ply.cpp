#include "ply.hpp"
#include "ply_helper.hpp"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <functional>
#include <locale>

namespace openMVG
{
namespace geometry
{
  namespace io
  {
    /**
    * @brief Export a set of points to a PLY file 
    * @param pts A list of points 
    * @param nor A list of normal (one for each points) 
    * @param col A list of color (one for each points) 
    * @param path Path of the file to save 
    * @param binary Indicate if the file should be written to binary format 
    * @retval true if load is successful 
    * @retval false if load fails 
    */
    bool PLYWrite( const std::vector<Vec3> &pts,
                   const std::vector<Vec3> *nor,
                   const std::vector<Vec3uc> *col,
                   const std::string &path,
                   const bool binary )
    {
      std::ofstream file( path );
      if ( !file )
      {
        return false;
      }

      // Write header
      file << "ply" << std::endl;
      if ( binary )
      {
        file << "format binary_little_endian 1.0" << std::endl;
      }
      else
      {
        file << "format ascii 1.0" << std::endl;
      }

      const size_t nb_vert = pts.size();
      file << "element vertex " << nb_vert << std::endl;
      file << "property double x" << std::endl;
      file << "property double y" << std::endl;
      file << "property double z" << std::endl;
      if ( nor )
      {
        file << "property double nx" << std::endl;
        file << "property double ny" << std::endl;
        file << "property double nz" << std::endl;
      }
      if ( col )
      {
        file << "property uchar red" << std::endl;
        file << "property uchar green" << std::endl;
        file << "property uchar blue" << std::endl;
      }
      file << "end_header" << std::endl;

      // now write the points
      EndianAgnosticWriter<PLY_ASCII> asc_writer;
      EndianAgnosticWriter<PLY_LITTLE_ENDIAN> le_writer;

      for ( size_t id_vert = 0; id_vert < nb_vert; ++id_vert )
      {
        if ( binary )
        {
          le_writer.Write( file, pts[ id_vert ] );
          if ( nor )
          {
            le_writer.Write( file, ( *nor )[ id_vert ] );
          }
          if ( col )
          {
            le_writer.Write( file, ( *col )[ id_vert ] );
          }
        }
        else
        {
          asc_writer.Write( file, pts[ id_vert ] );
          if ( nor )
          {
            asc_writer.Write( file, ( *nor )[ id_vert ] );
          }
          if ( col )
          {
            asc_writer.Write( file, ( *col )[ id_vert ] );
          }
          file << std::endl;
        }
      }
      return true;
    }

    /**
    * @brief Export a set of points to a PLY file 
    * @param pts A list of points 
    * @param path Path of the file to save 
    * @param binary Indicate if the file should be written to binary format 
    * @retval true if load is successful 
    * @retval false if load fails 
    */
    bool PLYWrite( const std::vector<Vec3> &pts,
                   const std::string &path,
                   const bool binary )
    {
      return PLYWrite( pts, nullptr, nullptr, path, binary );
    }

    /**
    * @brief Load a ply file and gets is content points 
    * @param path Path of the file to load 
    * @param[out] pts The list of points in the given file 
    * @retval true if load is successful 
    * @retval false if load fails 
    */
    bool PLYRead( const std::string &path,
                  std::vector<Vec3> &pts )
    {
      return PLYRead( path, pts, nullptr, nullptr );
    }

    // trim from start
    static inline std::string &ltrim( std::string &s )
    {
      s.erase( s.begin(), std::find_if( s.begin(), s.end(),
                                        std::not1( std::ptr_fun<int, int>( std::isspace ) ) ) );
      return s;
    }

    // trim from end
    static inline std::string &rtrim( std::string &s )
    {
      s.erase( std::find_if( s.rbegin(), s.rend(),
                             std::not1( std::ptr_fun<int, int>( std::isspace ) ) )
                   .base(),
               s.end() );
      return s;
    }

    // trim from both ends
    static inline std::string &trim( std::string &s )
    {
      return ltrim( rtrim( s ) );
    }

    // from : http://stackoverflow.com/questions/53849/how-do-i-tokenize-a-string-in-c
    static inline std::vector<std::string> tokenize( const std::string &str )
    {
      // construct a stream from the string
      std::stringstream strstr( str );

      // use stream iterators to copy the stream to the vector as whitespace separated strings
      std::istream_iterator<std::string> it( strstr );
      std::istream_iterator<std::string> end;
      std::vector<std::string> results( it, end );

      return results;
    }

    /**
    * @brief Load a ply file and gets is content points 
    * @param path Path of the file to load 
    * @param[out] pts The list of points in the given file 
    * @param[out] nor The list of normals in the given file 
    * @param[out] col The list of colors in the given file 
    * @note If nor or col are equal to nullptr, nothing is retrived 
    * @note If file does not contains nor or col, corresponding vectors will be empty  
    * @retval true if load is successful 
    * @retval false if load fails 
    */
    bool PLYRead( const std::string &path,
                  std::vector<Vec3> &pts,
                  std::vector<Vec3> *nor,
                  std::vector<Vec3uc> *col )
    {
      std::ifstream file( path );
      if ( !file )
      {
        return false;
      }

      // Read header
      std::string line;
      std::getline( file, line );
      if ( trim( line ) != "ply" )
      {
        return false;
      }

      int nb_coord_pts = 0;
      int nb_coord_nor = 0;
      int nb_coord_col = 0;

      size_t property_pts_size = 0;
      size_t property_nor_size = 0;
      size_t property_col_size = 0;

      size_t nb_elt = 0;

      ply_endianness endianness = PLY_ASCII;

      while ( 1 )
      {
        std::getline( file, line );
        std::vector<std::string> tokens = tokenize( line );
        if ( tokens.size() == 0 )
        {
          continue;
        }
        if ( tokens[ 0 ] == "comment" )
        {
          continue;
        }
        else if ( tokens[ 0 ] == "format" )
        {
          if ( tokens.size() > 1 )
          {
            if ( tokens[ 1 ] == "ascii" )
            {
              endianness = PLY_ASCII;
            }
            else if ( tokens[ 1 ] == "binary_little_endian" )
            {
              endianness = PLY_LITTLE_ENDIAN;
            }
            else if ( tokens[ 1 ] == "binary_big_endian" )
            {
              endianness = PLY_BIG_ENDIAN;
            }
            else
            {
              return false;
            }
          }
          else
          {
            return false;
          }
        }
        else if ( tokens[ 0 ] == "element" )
        {
          if ( tokens.size() > 2 )
          {
            if ( tokens[ 1 ] == "vertex" )
            {
              Convert( tokens[ 2 ], nb_elt );
            }
          }
          else
          {
            return false;
          }
        }
        else if ( tokens[ 0 ] == "property" )
        {
          size_t current_size;
          if ( tokens.size() > 2 )
          {
            // read size
            const std::string type = tokens[ 1 ];
            if ( type == "uchar" || type == "char" )
            {
              current_size = 1;
            }
            else if ( type == "short" || type == "ushort" )
            {
              current_size = 2;
            }
            else if ( type == "int" || type == "uint" )
            {
              current_size = 4;
            }
            else if ( type == "float" )
            {
              current_size = 4;
            }
            else if ( type == "double" )
            {
              current_size = 8;
            }
            // read kind (x,y,z, nx,ny,nz, red, green, blue )
            const std::string item = tokens[ 2 ];
            if ( item == "x" || item == "y" || item == "z" )
            {
              property_pts_size = current_size;
              ++nb_coord_pts;
            }
            else if ( item == "nx" || item == "ny" || item == "nz" )
            {
              property_nor_size = current_size;
              ++nb_coord_nor;
            }
            else if ( item == "red" || item == "green" || item == "blue" )
            {
              property_col_size = current_size;
              ++nb_coord_col;
            }
          }
        }
        else if ( tokens[ 0 ] == "end_header" )
        {
          break;
        }
      }

      if ( nb_coord_pts != 3 )
      {
        std::cout << "nb_coord_pts != 3" << std::endl;
        return false;
      }
      if ( !( nb_coord_nor == 0 || nb_coord_nor == 3 ) )
      {
        std::cout << "!( nb_coord_nor == 0 || nb_coord_nor == 3 )" << __LINE__ << std::endl;
        return false;
      }
      if ( !( nb_coord_col == 0 || nb_coord_col == 3 ) )
      {
        std::cout << "!( nb_coord_col == 0 || nb_coord_col == 3 )" << __LINE__ << std::endl;
        return false;
      }
      if ( nb_elt == 0 )
      {
        std::cout << "nb_elt == 0" << __LINE__ << std::endl;
        return false;
      }
      if ( nb_coord_pts == 3 && property_pts_size == 0 )
      {
        std::cout << "nb_coord_pts == 3 && property_pts_size == 0" << __LINE__ << std::endl;
        return false;
      }
      if ( nb_coord_col == 3 && property_col_size == 0 )
      {
        std::cout << "nb_coord_col == 3 && property_col_size == 0" << __LINE__ << std::endl;
        return false;
      }
      if ( nb_coord_nor == 3 && property_nor_size == 0 )
      {
        std::cout << "nb_coord_nor == 3 && property_nor_size == 0" << __LINE__ << std::endl;
        return false;
      }

      // resize points
      pts.resize( nb_elt );
      if ( nor )
      {
        if ( nb_coord_nor == 3 )
        {
          nor->resize( nb_elt );
        }
        else
        {
          nor->clear();
        }
      }

      if ( col )
      {
        if ( nb_coord_col == 3 )
        {
          col->resize( nb_elt );
        }
        else
        {
          col->clear();
        }
      }

      EndianAgnosticReader<PLY_ASCII> ascii_reader;
      EndianAgnosticReader<PLY_LITTLE_ENDIAN> le_reader;
      EndianAgnosticReader<PLY_BIG_ENDIAN> be_reader;

      // now read the points
      for ( size_t id_point = 0; id_point < nb_elt; ++id_point )
      {
        if ( endianness == PLY_ASCII )
        {
          ascii_reader.Read( file, pts[ id_point ] );
          if ( nor && nb_coord_nor > 0 )
          {
            ascii_reader.Read( file, ( *nor )[ id_point ] );
          }
          if ( col && nb_coord_col > 0 )
          {
            ascii_reader.Read( file, ( *col )[ id_point ] );
          }
          // Skip everything until
          std::string tmp;
          std::getline( file, tmp );
        }
        else if ( endianness == PLY_LITTLE_ENDIAN )
        {
          le_reader.Read( file, pts[ id_point ] );
          if ( nor && nb_coord_pts > 0 )
          {
            le_reader.Read( file, ( *nor )[ id_point ] );
          }
          if ( col && nb_coord_col > 0 )
          {
            le_reader.Read( file, ( *col )[ id_point ] );
          }
        }
        else if ( endianness == PLY_BIG_ENDIAN )
        {
          be_reader.Read( file, pts[ id_point ] );
          if ( nor && nb_coord_pts > 0 )
          {
            be_reader.Read( file, ( *nor )[ id_point ] );
          }
          if ( col && nb_coord_col > 0 )
          {
            be_reader.Read( file, ( *col )[ id_point ] );
          }
        }
      }

      return true;
    }

  } // namespace io
} // namespace geometry
} // namespace openMVG