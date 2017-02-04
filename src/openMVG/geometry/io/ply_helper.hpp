#ifndef OPENMVG_GEOMETRY_IO_PLY_HELPER_HPP
#define OPENMVG_GEOMETRY_IO_PLY_HELPER_HPP

#include "openMVG/numeric/numeric.h"

#include <fstream>

namespace openMVG
{
namespace geometry
{
  namespace io
  {
    /**
    * @brief Store different kind of endianess that can be found on a PLY file 
    */
    enum ply_endianness
    {
      PLY_ASCII,
      PLY_BIG_ENDIAN,
      PLY_LITTLE_ENDIAN
    };

    /**
    * @brief Compute runtime endianess 
    * @return Current system endianess 
    */
    static inline ply_endianness GetSystemEndianness()
    {
      short int tmp = 0x1;
      char *first   = (char *)&tmp;
      return ( first[ 0 ] == 1 ) ? PLY_LITTLE_ENDIAN : PLY_BIG_ENDIAN;
    }

    /**
    * @brief Fuction used to swap bytes (from LE to BE)
    * @param val Value to inverse 
    */
    template <typename T>
    T ByteSwap( T val )
    {
      T retVal;
      char *pVal    = (char *)&val;
      char *pRetVal = (char *)&retVal;
      int size      = sizeof( T );
      for ( int i = 0; i < size; i++ )
      {
        pRetVal[ size - 1 - i ] = pVal[ i ];
      }

      return retVal;
    }

    /**
    * @brief Class used to write a value using a specified PLY endianness 
    */
    template <int Endianness>
    struct EndianAgnosticWriter
    {
    public:
      /**
      * @brief Ctr 
      */
      EndianAgnosticWriter();

      /**
      * @brief Write a 3d vector of double in a file 
      * @param file File in which data will be written 
      * @param vec vector to write 
      */
      void Write( std::ofstream &file, const Vec3 &vec );

      /**
      * @brief Write a 3d vector of unsigned char in a file 
      * @param file File in which data will be written 
      * @param vec vector to write 
      */
      void Write( std::ofstream &file, const Vec3uc &vec );

    private:
      /// Current system endianess
      ply_endianness m_system_endianness;
    };

    // Specialization for ascii
    template <>
    struct EndianAgnosticWriter<PLY_ASCII>
    {
    public:
      /**
      * @brief Ctr 
      */
      EndianAgnosticWriter()
          : m_system_endianness( GetSystemEndianness() )
      {
      }

      /**
      * @brief Write a 3d vector of double in a file 
      * @param file File in which data will be written 
      * @param vec vector to write 
      */
      void Write( std::ofstream &file, const Vec3 &vec )
      {
        file << vec[ 0 ] << " " << vec[ 1 ] << " " << vec[ 2 ] << " ";
      }

      /**
      * @brief Write a 3d vector of unsigned char in a file 
      * @param file File in which data will be written 
      * @param vec vector to write 
      */
      void Write( std::ofstream &file, const Vec3uc &vec )
      {
        file << (int)vec[ 0 ] << " " << (int)vec[ 1 ] << " " << (int)vec[ 2 ] << " ";
      }

    private:
      /// Current system endianness
      ply_endianness m_system_endianness;
    };

    template <>
    struct EndianAgnosticWriter<PLY_LITTLE_ENDIAN>
    {
    public:
      /**
      * @brief Ctr 
      */
      EndianAgnosticWriter()
          : m_system_endianness( GetSystemEndianness() )
      {
      }

      /**
      * @brief Write a 3d vector of double in a file 
      * @param file File in which data will be written 
      * @param vec vector to write 
      */
      void Write( std::ofstream &file, const Vec3 &vec )
      {
        if ( m_system_endianness == PLY_LITTLE_ENDIAN )
        {
          // No conversion
          file.write( reinterpret_cast<const char *>( vec.data() ), sizeof( Vec3 ) );
        }
        else
        {
          // Conversion
          const Vec3 r( ByteSwap( vec[ 0 ] ), ByteSwap( vec[ 1 ] ), ByteSwap( vec[ 2 ] ) );
          file.write( reinterpret_cast<const char *>( r.data() ), sizeof( Vec3 ) );
        }
      }

      /**
      * @brief Write a 3d vector of unsigned char in a file 
      * @param file File in which data will be written 
      * @param vec vector to write 
      */
      void Write( std::ofstream &file, const Vec3uc &vec )
      {
        // No reverse because sizeof char == 1
        file.write( reinterpret_cast<const char *>( vec.data() ), sizeof( Vec3uc ) );
      }

    private:
      /// Current system endianness
      ply_endianness m_system_endianness;
    };

    template <>
    struct EndianAgnosticWriter<PLY_BIG_ENDIAN>
    {
    public:
      /**
      * @brief Ctr 
      */
      EndianAgnosticWriter()
          : m_system_endianness( GetSystemEndianness() )
      {
      }

      /**
      * @brief Write a 3d vector of double in a file 
      * @param file File in which data will be written 
      * @param vec vector to write 
      */
      void Write( std::ofstream &file, const Vec3 &vec )
      {
        if ( m_system_endianness == PLY_BIG_ENDIAN )
        {
          // No conversion
          file.write( reinterpret_cast<const char *>( vec.data() ), sizeof( Vec3 ) );
        }
        else
        {
          // Conversion
          const Vec3 r( ByteSwap( vec[ 0 ] ), ByteSwap( vec[ 1 ] ), ByteSwap( vec[ 2 ] ) );
          file.write( reinterpret_cast<const char *>( r.data() ), sizeof( Vec3 ) );
        }
      }

      /**
      * @brief Write a 3d vector of unsigned char in a file 
      * @param file File in which data will be written 
      * @param vec vector to write 
      */
      void Write( std::ofstream &file, const Vec3uc &vec )
      {
        // No reverse because sizeof char == 1
        file.write( reinterpret_cast<const char *>( vec.data() ), sizeof( Vec3uc ) );
      }

    private:
      /// Current system endianness
      ply_endianness m_system_endianness;
    };

    /**
    * @brief Convert a string value to a numeric value 
    * @param str String in which a numer is present 
    * @param[out] val Value to store 
    * @retval true if conversion is ok 
    * @retval false if conversion fails 
    */
    template <typename T>
    bool Convert( const std::string &str, T &val )
    {
      std::stringstream sstr( str );
      sstr >> val;
      return !sstr.fail();
    }

    // extract double value from raw binary array
    void ExtractDoubleFromRawArray( const char *data, const size_t elt_size, double &x, double &y, double &z )
    {
      if ( elt_size == 4 )
      {
        float *fdata = reinterpret_cast<float *>( const_cast<char *>( data ) );
        x            = fdata[ 0 ];
        y            = fdata[ 1 ];
        z            = fdata[ 2 ];
      }
      else if ( elt_size == 8 )
      {
        double *ddata = reinterpret_cast<double *>( const_cast<char *>( data ) );
        x             = ddata[ 0 ];
        y             = ddata[ 1 ];
        z             = ddata[ 2 ];
      }
    }

    /**
    * @brief Class used to read Vector data in agnostic way  
    */
    template <int Endianness>
    class EndianAgnosticReader
    {
    public:
      /**
        * @brief ctr 
        */
      EndianAgnosticReader();

      /**
        * @brief Read vector data in double format 
        * @param file Stream in which data is read 
        * @param[out] vec output vector 
        * @param elt_byte_size Byte size of each component to read 
        * @retval true if read is ok 
        * @retval false if read fails  
        */
      bool Read( std::ifstream &file, Vec3 &vec, const std::vector<int> &elt_byte_size );

      /**
        * @brief Read vector data in unsigned char format 
        * @param file Stream in which data is read 
        * @param elt_byte_size Byte size of each component to read 
        * @param[out] vec output vector 
        * @retval true if read is ok 
        * @retval false if read fails  
        */
      bool Read( std::ifstream &file, Vec3uc &vec, const std::vector<int> &elt_byte_size );

    private:
      /// Current system
      ply_endianness m_system_endianness;
    };

    /**
    * @brief Class used to read Vector data in agnostic way  
    * -> Specialization for ascii data 
    */
    template <>
    class EndianAgnosticReader<PLY_ASCII>
    {
    public:
      /**
        * @brief ctr 
        */
      EndianAgnosticReader()
          : m_system_endianness( GetSystemEndianness() )
      {
      }

      /**
        * @brief Read vector data in double format 
        * @param file Stream in which data is read 
        * @param[out] vec output vector 
        * @param elt_byte_size Byte size of each component to read 
        * @retval true if read is ok 
        * @retval false if read fails  
        */
      bool Read( std::ifstream &file, Vec3 &vec, const std::vector<int> &elt_byte_size )
      {
        double x, y, z;
        if ( file >> x >> y >> z )
        {
          vec = Vec3( x, y, z );
          return true;
        }
        else
        {
          return false;
        }
      }

      /**
        * @brief Read vector data in unsigned char format 
        * @param file Stream in which data is read 
        * @param elt_byte_size Byte size of each component to read 
        * @param[out] vec output vector 
        * @retval true if read is ok 
        * @retval false if read fails  
        */
      bool Read( std::ifstream &file, Vec3uc &vec, const std::vector<int> &elt_byte_size )
      {
        int a, b, c;
        if ( file >> a >> b >> c )
        {
          vec = Vec3uc( a, b, c );
          return true;
        }
        else
        {
          return false;
        }
      }

    private:
      /// Current system
      ply_endianness m_system_endianness;
    };

    /**
    * @brief Class used to read Vector data in agnostic way  
    * -> Specialization for ascii data 
    */
    template <>
    class EndianAgnosticReader<PLY_LITTLE_ENDIAN>
    {
    public:
      /**
        * @brief ctr 
        */
      EndianAgnosticReader()
          : m_system_endianness( GetSystemEndianness() )
      {
      }

      /**
        * @brief Read vector data in double format 
        * @param file Stream in which data is read 
        * @param[out] vec output vector 
        * @param elt_byte_size Byte size of each component to read 
        * @retval true if read is ok 
        * @retval false if read fails  
        */
      bool Read( std::ifstream &file, Vec3 &vec, const std::vector<int> &elt_byte_size )
      {
        char raw[ 3 * 8 ];
        file.read( reinterpret_cast<char *>( raw ), 3 * sizeof( double ) );
        if ( !file )
        {
          return false;
        }
        double data[3] ; 
        ExtractDoubleFromRawArray( raw, elt_byte_size[ 0 ], data[ 0 ], data[ 1 ], data[ 2 ] );

        if ( m_system_endianness != PLY_LITTLE_ENDIAN )
        {
          data[ 0 ] = ByteSwap( data[ 0 ] );
          data[ 1 ] = ByteSwap( data[ 1 ] );
          data[ 2 ] = ByteSwap( data[ 2 ] );
        }
        vec = Vec3( data[ 0 ], data[ 1 ], data[ 2 ] );
        return true;
      }

      /**
        * @brief Read vector data in unsigned char format 
        * @param file Stream in which data is read 
        * @param elt_byte_size Byte size of each component to read 
        * @param[out] vec output vector 
        * @retval true if read is ok 
        * @retval false if read fails  
        */
      bool Read( std::ifstream &file, Vec3uc &vec, const std::vector<int> &elt_byte_size )
      {
        unsigned char data[ 3 ];
        file.read( reinterpret_cast<char *>( data ), 3 * sizeof( unsigned char ) );
        if ( !file )
        {
          return false;
        }
        vec = Vec3uc( data[ 0 ], data[ 1 ], data[ 2 ] );
        return true;
      }

    private:
      /// Current system
      ply_endianness m_system_endianness;
    };

    /**
    * @brief Class used to read Vector data in agnostic way  
    * -> Specialization for ascii data 
    */
    template <>
    class EndianAgnosticReader<PLY_BIG_ENDIAN>
    {
    public:
      /**
        * @brief ctr 
        */
      EndianAgnosticReader()
          : m_system_endianness( GetSystemEndianness() )
      {
      }

      /**
        * @brief Read vector data in double format 
        * @param file Stream in which data is read 
        * @param[out] vec output vector 
        * @param elt_byte_size Byte size of each component to read 
        * @retval true if read is ok 
        * @retval false if read fails  
        */
      bool Read( std::ifstream &file, Vec3 &vec, const std::vector<int> &elt_byte_size )
      {
        char raw[ 3 * 8 ];
        file.read( reinterpret_cast<char *>( raw ), 3 * sizeof( double ) );
        if ( !file )
        {
          return false;
        }
        double data[3] ; 
        ExtractDoubleFromRawArray( raw, elt_byte_size[ 0 ], data[ 0 ], data[ 1 ], data[ 2 ] );

        if ( m_system_endianness != PLY_BIG_ENDIAN )
        {
          data[ 0 ] = ByteSwap( data[ 0 ] );
          data[ 1 ] = ByteSwap( data[ 1 ] );
          data[ 2 ] = ByteSwap( data[ 2 ] );
        }
        vec = Vec3( data[ 0 ], data[ 1 ], data[ 2 ] );

        return true;
      }

      /**
        * @brief Read vector data in unsigned char format 
        * @param file Stream in which data is read 
        * @param elt_byte_size Byte size of each component to read 
        * @param[out] vec output vector 
        * @retval true if read is ok 
        * @retval false if read fails  
        */
      bool Read( std::ifstream &file, Vec3uc &vec, const std::vector<int> &elt_byte_size )
      {
        unsigned char data[ 3 ];
        file.read( reinterpret_cast<char *>( data ), 3 * sizeof( unsigned char ) );
        if ( !file )
        {
          return false;
        }
        vec = Vec3uc( data[ 0 ], data[ 1 ], data[ 2 ] );
        return true;
      }

    private:
      /// Current system
      ply_endianness m_system_endianness;
    };

  } // namespace io
} // namespace geometry
} // namespace openMVG

#endif