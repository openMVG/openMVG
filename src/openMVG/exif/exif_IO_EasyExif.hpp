// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2013-2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_EXIF_EXIF_IO_EASYEXIF_HPP
#define OPENMVG_EXIF_EXIF_IO_EASYEXIF_HPP

#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include "openMVG/exif/exif_IO.hpp"
#include "third_party/easyexif/exif.h"

namespace openMVG
{
namespace exif
{

/**
* Remove all leading and trailing spaces from the input. The result is a trimmed copy of the input
* @todo move this class elsewhere since it's not exif related
*/
inline std::string trim_copy( const std::string& s )
{
  if (s.empty() )
  {
    return s;
  }

  std::string res( s );
  // remove leading and trailing spaces
  res.erase( 0, res.find_first_not_of( ' ' ) );
  res.erase( res.find_last_not_of( ' ' ) + 1 );
  // handle multiple trailing end character
  res = res.substr( 0, res.find( '\0' ) );
  return res;
}

/**
* @brief Class managing EXIF data using EasyEXIF library
*/
class Exif_IO_EasyExif : public Exif_IO
{
  public:

    /**
    * @brief Default constructor
    */
    Exif_IO_EasyExif(): bHaveExifInfo_( false )
    {
    }

    /**
    * @brief Constructor using a file
    * @param sFileName path of the image to analyze
    */
    explicit Exif_IO_EasyExif( const std::string & sFileName ): bHaveExifInfo_( false )
    {
      open( sFileName );
    }

    /**
    * @brief Open and populate EXIF data
    * @param sFileName path of the image to analyze
    * @retval true if image file could be parsed correctly
    * @retval false if image file could not be parsed and analysed
    */
    bool open( const std::string & sFileName ) override
    {
      // Read the file into a buffer
      FILE *fp = fopen( sFileName.c_str(), "rb" );
      if ( !fp )
      {
        return false;
      }
      fseek( fp, 0, SEEK_END );
      unsigned long fsize = ftell( fp );
      rewind( fp );
      std::vector<unsigned char> buf( fsize );
      if ( fread( &buf[0], 1, fsize, fp ) != fsize )
      {
        fclose( fp );
        return false;
      }
      fclose( fp );

      // Parse EXIF
      bHaveExifInfo_ = ( exifInfo_.parseFrom( &buf[0], fsize ) == PARSE_EXIF_SUCCESS );

      return bHaveExifInfo_;
    }

    /**
    * @brief Get image width
    * @return Width of the image (in pixel)
    */
    size_t getWidth() const override
    {
      return exifInfo_.ImageWidth;
    }

    /**
    * @brief Get image height
    * @return Height of the image (in pixel)
    */
    size_t getHeight() const override
    {
      return exifInfo_.ImageHeight;
    }

    /**
    * @brief Get Focal (in mm)
    * @return Focal of the lens when image was shot (in mm)
    */
    float getFocal() const override
    {
      return static_cast<float>( exifInfo_.FocalLength );
    }

    /**
    * @brief Get Brand of the camera
    * @return Brand name
    */
    std::string getBrand() const override
    {
      return trim_copy( exifInfo_.Make );
    }

    /**
    * @brief Get Model of the camera
    * @return Camera model name
    */
    std::string getModel() const override
    {
      return trim_copy( exifInfo_.Model );
    }

    /**
    * @brief Get Lens model
    * @return Lens model name
    */
    std::string getLensModel() const override
    {
      return trim_copy( exifInfo_.LensInfo.Model );
    }

    /**
    * @brief Get a unique indentifier for this image
    * @return Unique ID
    */
    std::string getImageUniqueID() const override
    {
      return exifInfo_.ImageUniqueID;
    }

    /**
      * @brief Verify if the file has metadata
      * @retval true If EXIF value are present
      * @retval false If EXIF value are not present
      */
    bool doesHaveExifInfo() const override
    {
      return bHaveExifInfo_;
    }

    /**
    * @brief Print all data
    * @return string containing all EXIF data
    */
    std::string allExifData() const override
    {
      std::ostringstream os;
      os
          << "Camera make       : " << exifInfo_.Make << "\n"
          << "Camera model      : " << exifInfo_.Model << "\n"
          << "Software          : " << exifInfo_.Software << "\n"
          << "Bits per sample   : " << exifInfo_.BitsPerSample << "\n"
          << "Image width       : " << exifInfo_.ImageWidth << "\n"
          << "Image height      : " << exifInfo_.ImageHeight << "\n"
          << "Image description : " << exifInfo_.ImageDescription << "\n"
          << "Image orientation : " << exifInfo_.Orientation << "\n"
          << "Image copyright   : " << exifInfo_.Copyright << "\n"
          << "Image date/time   : " << exifInfo_.DateTime << "\n"
          << "Original date/time: " << exifInfo_.DateTimeOriginal << "\n"
          << "Digitize date/time: " << exifInfo_.DateTimeDigitized << "\n"
          << "Subsecond time    : " << exifInfo_.SubSecTimeOriginal << "\n"
          << "Exposure time     : 1/" << ( unsigned ) ( 1.0 / exifInfo_.ExposureTime ) << "\n"
          << "F-stop            : " << exifInfo_.FNumber << "\n"
          << "ISO speed         : " << exifInfo_.ISOSpeedRatings << "\n"
          << "Subject distance  : " << exifInfo_.SubjectDistance << "\n"
          << "Exposure bias     : EV" << exifInfo_.ExposureBiasValue << "\n"
          << "Flash used?       : " << exifInfo_.Flash << "\n"
          << "Metering mode     : " << exifInfo_.MeteringMode << "\n"
          << "Lens focal length : mm\n" << exifInfo_.FocalLength << "\n"
          << "35mm focal length : mm\n" << exifInfo_.FocalLengthIn35mm << "\n"
          << "GPS Latitude      : deg ( deg, min, sec )\n" << "("
          <<  exifInfo_.GeoLocation.Latitude << ", "
          <<  exifInfo_.GeoLocation.LatComponents.degrees << ", "
          <<  exifInfo_.GeoLocation.LatComponents.minutes << ", "
          <<  exifInfo_.GeoLocation.LatComponents.seconds << ", "
          <<  exifInfo_.GeoLocation.LatComponents.direction << ")" << "\n"
          << "GPS Longitude      : deg ( deg, min, sec )\n" << "("
          <<  exifInfo_.GeoLocation.Longitude << ", "
          <<  exifInfo_.GeoLocation.LonComponents.degrees << ", "
          <<  exifInfo_.GeoLocation.LonComponents.minutes << ", "
          <<  exifInfo_.GeoLocation.LonComponents.seconds << ", "
          <<  exifInfo_.GeoLocation.LonComponents.direction << ")" << "\n"
          << "GPS Altitude       : m" << exifInfo_.GeoLocation.Altitude << "\n"
          << "Lens stop (min, max) : " << "("
          << exifInfo_.LensInfo.FStopMin << ", "
          << exifInfo_.LensInfo.FStopMax << ")"
          << "Lens focal (min, max) : " << "("
          << exifInfo_.LensInfo.FocalLengthMin << ", "
          << exifInfo_.LensInfo.FocalLengthMax << ")"
          << "Lens make : " <<  exifInfo_.LensInfo.Make << "\n"
          << "Lens model : " << exifInfo_.LensInfo.Model << "\n"
          << "Image Unique ID    : " << exifInfo_.ImageUniqueID << "\n";
      return os.str();
    }

        /**
    * @brief Try to read and save the EXIF GPS latitude
    * @return If GPS Latitude can be read & exported, return true
    */
    bool GPSLatitude(double * latitude) const override
    {
      if (exifInfo_.GeoLocation.Latitude != std::numeric_limits<double>::infinity())
      {
        (*latitude) = exifInfo_.GeoLocation.Latitude;
        return true;
      }
      return false;
    }

    /**
    * @brief Try to read and save the EXIF GPS longitude
    * @return If GPS Longitude can be read & exported, return true
    */
    bool GPSLongitude(double * longitude) const override
    {
      if (exifInfo_.GeoLocation.Longitude != std::numeric_limits<double>::infinity())
      {
        (*longitude) = exifInfo_.GeoLocation.Longitude;
        return true;
      }
      return false;
    }

   /**
    * @brief Try to read and save the EXIF GPS altitude
    * @return If GPS Altitude can be read & exported, return true
    */
    bool GPSAltitude(double * altitude) const  override
    {
      if (exifInfo_.GeoLocation.Altitude != std::numeric_limits<double>::infinity())
      {
        (*altitude) = exifInfo_.GeoLocation.Altitude;
        return true;
      }
      return false;
    }


  private:

    /// Internal data storing all exif data
    easyexif::EXIFInfo exifInfo_;

    /// Indicate if exifInfo_ is populated
    bool bHaveExifInfo_;
};

} // namespace exif
} // namespace openMVG

#endif // OPENMVG_EXIF_EXIF_IO_EASYEXIF_HPP
