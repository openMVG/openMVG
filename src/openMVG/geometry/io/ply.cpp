#include "ply.hpp"
#include "ply_helper.hpp"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <functional>
#include <iterator>
#include <locale>
#include <string>

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

    // List of valid element in the file
    enum
    {
      // Points
      PLY_PT_X,
      PLY_PT_Y,
      PLY_PT_Z,
      // Normals
      PLY_NOR_X,
      PLY_NOR_Y,
      PLY_NOR_Z,
      // Color
      PLY_COL_R,
      PLY_COL_G,
      PLY_COL_B,
      // Undefined value
      PLY_UNDEFINED
    };

    bool CheckCanonicalOrder( const std::vector<int> &element_order )
    {
      // Check PT_X PT_Y PT_Z
      if ( element_order.size() == 3 )
      {
        return element_order[ 0 ] == PLY_PT_X &&
               element_order[ 1 ] == PLY_PT_Y &&
               element_order[ 2 ] == PLY_PT_Z;
      }

      // Check PT_X PT_Y PT_Z NOR_X NOR_Y NOR_Z or PT_X PT_Y PT_Z COL_R COL_G COL_B
      if ( element_order.size() == 6 )
      {
        return ( ( element_order[ 0 ] == PLY_PT_X &&
                   element_order[ 1 ] == PLY_PT_Y &&
                   element_order[ 2 ] == PLY_PT_Z &&
                   element_order[ 3 ] == PLY_NOR_X &&
                   element_order[ 4 ] == PLY_NOR_Y &&
                   element_order[ 5 ] == PLY_NOR_Z ) ||
                 ( element_order[ 0 ] == PLY_PT_X &&
                   element_order[ 1 ] == PLY_PT_Y &&
                   element_order[ 2 ] == PLY_PT_Z &&
                   element_order[ 3 ] == PLY_COL_R &&
                   element_order[ 4 ] == PLY_COL_G &&
                   element_order[ 5 ] == PLY_COL_B ) );
      }

      // All PT_X PT_Y PT_Z NOR_X NOR_Y NOR_Z COL_R COL_G COL_B
      if ( element_order.size() == 9 )
      {
        return element_order[ 0 ] == PLY_PT_X &&
               element_order[ 1 ] == PLY_PT_Y &&
               element_order[ 2 ] == PLY_PT_Z &&
               element_order[ 3 ] == PLY_NOR_X &&
               element_order[ 4 ] == PLY_NOR_Y &&
               element_order[ 5 ] == PLY_NOR_Z &&
               element_order[ 6 ] == PLY_COL_R &&
               element_order[ 7 ] == PLY_COL_G &&
               element_order[ 8 ] == PLY_COL_B;
      }

      // all other configurations are invalid (not really but we do not support it)
      return false;
    }

    /**
    * @brief Check if elements have consistent size (ie: component of a vector have same size) and if size is valid 
    * @param element_size Vector of size of each property 
    * @retval true id everything is ok 
    * @retval false if element are invalid 
    */
    bool CheckConsistentSize( const std::vector<int> &element_size )
    {
      if ( element_size.size() % 3 != 0 )
      {
        return false;
      }

      for ( int id_plate = 0; id_plate < element_size.size() / 3; ++id_plate )
      {
        const int s1 = element_size[ 3 * id_plate ];
        const int s2 = element_size[ 3 * id_plate + 1 ];
        const int s3 = element_size[ 3 * id_plate + 2 ];

        if ( !( s1 == s2 && s2 == s3 ) )
        {
          return false;
        }
        if ( s1 <= 0 || s1 > 8 )
        {
          return false;
        }
      }

      return true;
    }

    /**
    * @brief Indicate if a color is present in the list 
    * @retval true if a color is present 
    * @retval if no color element is present 
    */
    bool CheckHasColor( const std::vector<int> &element_order )
    {
      for ( size_t id_elt = 0; id_elt < element_order.size(); ++id_elt )
      {
        if ( element_order[ id_elt ] == PLY_COL_R ||
             element_order[ id_elt ] == PLY_COL_G ||
             element_order[ id_elt ] == PLY_COL_B )
        {
          return true;
        }
      }
      return false;
    }

    /**
    * @brief indicate if a sequence contain a normal information 
    * @retval true if a normal is present 
    * @retval if no normal element is present 
    */
    static inline bool CheckHasNormal( const std::vector<int> &element_order )
    {
      for ( size_t id_elt = 0; id_elt < element_order.size(); ++id_elt )
      {
        if ( element_order[ id_elt ] == PLY_NOR_X ||
             element_order[ id_elt ] == PLY_NOR_Y ||
             element_order[ id_elt ] == PLY_NOR_Z )
        {
          return true;
        }
      }
      return false;
    }

    /**
    * @brief Get byte size of each point component 
    * @param element_order Semantic element  
    * @param element_size Element size 
    * @return size of each component (in byte)
    * @note if nothing is present, return an empty vector 
    * @note Assume consistency (ie: vectors elements are consecutives)
    * @note Assume element_order and element_size at the same size
    */
    static inline std::vector<int> ExtractPointByteSize( const std::vector<int> &element_order,
                                                         const std::vector<int> &element_size )
    {
      std::vector<int> res;

      for ( int id_elt = 0; id_elt < element_order.size(); ++id_elt )
      {
        if ( element_order[ id_elt ] == PLY_PT_X )
        {
          res.push_back( element_size[ id_elt ] );
          res.push_back( element_size[ id_elt + 1 ] );
          res.push_back( element_size[ id_elt + 2 ] );
          return res;
        }
      }

      return res;
    }

    /**
    * @brief Get byte size of each normal component 
    * @param element_order Semantic element  
    * @param element_size Element size 
    * @return size of each component (in byte)
    * @note if nothing is present, return an empty vector 
    * @note Assume consistency (ie: vectors elements are consecutives)
    * @note Assume element_order and element_size at the same size
    */
    static inline std::vector<int> ExtractNormalByteSize( const std::vector<int> &element_order,
                                                          const std::vector<int> &element_size )
    {
      std::vector<int> res;

      for ( int id_elt = 0; id_elt < element_order.size(); ++id_elt )
      {
        if ( element_order[ id_elt ] == PLY_NOR_X )
        {
          res.push_back( element_size[ id_elt ] );
          res.push_back( element_size[ id_elt + 1 ] );
          res.push_back( element_size[ id_elt + 2 ] );
          return res;
        }
      }

      return res;
    }

    /**
    * @brief Get byte size of each color component 
    * @param element_order Semantic element  
    * @param element_size Element size 
    * @return size of each component (in byte)
    * @note if nothing is present, return an empty vector 
    * @note Assume consistency (ie: vectors elements are consecutives)
    * @note Assume element_order and element_size at the same size
    */
    static inline std::vector<int> ExtractColorByteSize( const std::vector<int> &element_order,
                                                         const std::vector<int> &element_size )
    {
      std::vector<int> res;

      for ( int id_elt = 0; id_elt < element_order.size(); ++id_elt )
      {
        if ( element_order[ id_elt ] == PLY_COL_R )
        {
          res.push_back( element_size[ id_elt ] );
          res.push_back( element_size[ id_elt + 1 ] );
          res.push_back( element_size[ id_elt + 2 ] );
          return res;
        }
      }

      return res;
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

      size_t nb_elt = 0;

      ply_endianness endianness = PLY_ASCII;

      std::vector<int> element_order;
      std::vector<int> element_size;

      while ( 1 )
      {
        if ( !std::getline( file, line ) )
        {
          break;
        }

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
              const bool ok = Convert( tokens[ 2 ], nb_elt );
              if ( !ok )
              {
                return false;
              }
            }
          }
          else
          {
            return false;
          }
        }
        else if ( tokens[ 0 ] == "property" )
        {
          if ( tokens.size() > 2 )
          {
            // read size
            int current_size       = -1;
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
            element_size.push_back( current_size );

            // read kind (x,y,z, nx,ny,nz, red, green, blue )
            const std::string item = tokens[ 2 ];
            if ( item == "x" )
            {
              element_order.push_back( PLY_PT_X );
            }
            else if ( item == "y" )
            {
              element_order.push_back( PLY_PT_Y );
            }
            else if ( item == "z" )
            {
              element_order.push_back( PLY_PT_Z );
            }
            else if ( item == "nx" )
            {
              element_order.push_back( PLY_NOR_X );
            }
            else if ( item == "ny" )
            {
              element_order.push_back( PLY_NOR_Y );
            }
            else if ( item == "nz" )
            {
              element_order.push_back( PLY_NOR_Z );
            }
            else if ( item == "red" )
            {
              element_order.push_back( PLY_COL_R );
            }
            else if ( item == "green" )
            {
              element_order.push_back( PLY_COL_G );
            }
            else if ( item == "blue" )
            {
              element_order.push_back( PLY_COL_B );
            }
          }
        }
        else if ( tokens[ 0 ] == "end_header" )
        {
          break;
        }
      }

      // Check if we have a canonical ordering
      if ( !CheckCanonicalOrder( element_order ) )
      {
        return false;
      }
      // Check if all elements in the same feature have same size
      if ( !CheckConsistentSize( element_size ) )
      {
        return false;
      }
      if ( nb_elt == 0 )
      {
        return false;
      }

      const bool has_color  = CheckHasColor( element_order );
      const bool has_normal = CheckHasNormal( element_order );

      const std::vector<int> pts_byte_size = ExtractPointByteSize( element_order, element_size );
      const std::vector<int> nor_byte_size = ExtractNormalByteSize( element_order, element_size );
      const std::vector<int> col_byte_size = ExtractColorByteSize( element_order, element_size );

      // resize points
      pts.resize( nb_elt );
      if ( nor )
      {
        if ( has_normal )
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
        if ( has_color )
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
      bool ok = true;
      for ( size_t id_point = 0; id_point < nb_elt; ++id_point )
      {
        if ( endianness == PLY_ASCII )
        {
          ok = ascii_reader.Read( file, pts[ id_point ], pts_byte_size );
          if ( !ok )
          {
            return false;
          }
          if ( nor && has_normal )
          {
            ok = ascii_reader.Read( file, ( *nor )[ id_point ], nor_byte_size );
            if ( !ok )
            {
              return false;
            }
          }
          if ( col && has_color )
          {
            ok = ascii_reader.Read( file, ( *col )[ id_point ], col_byte_size );
            if ( !ok )
            {
              return false;
            }
          }
          // Skip everything until
          std::string tmp;
          std::getline( file, tmp );
        }
        else if ( endianness == PLY_LITTLE_ENDIAN )
        {
          ok = le_reader.Read( file, pts[ id_point ], pts_byte_size );
          if ( !ok )
          {
            return false;
          }

          if ( nor && has_normal )
          {
            ok = le_reader.Read( file, ( *nor )[ id_point ], nor_byte_size );
            if ( !ok )
            {
              return false;
            }
          }
          if ( col && has_color )
          {
            ok = le_reader.Read( file, ( *col )[ id_point ], col_byte_size );
            if ( !ok )
            {
              return false;
            }
          }
        }
        else if ( endianness == PLY_BIG_ENDIAN )
        {
          ok = be_reader.Read( file, pts[ id_point ], pts_byte_size );
          if ( !ok )
          {
            return false;
          }
          if ( nor && has_normal )
          {
            ok = be_reader.Read( file, ( *nor )[ id_point ], nor_byte_size );
            if ( !ok )
            {
              return false;
            }
          }
          if ( col && has_color )
          {
            ok = be_reader.Read( file, ( *col )[ id_point ], col_byte_size );
            if ( !ok )
            {
              return false;
            }
          }
        }
      }

      return true;
    }

  } // namespace io
} // namespace geometry
} // namespace openMVG