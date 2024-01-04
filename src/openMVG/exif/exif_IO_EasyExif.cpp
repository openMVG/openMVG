// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2013-2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/exif/exif_IO_EasyExif.hpp"

#include <fstream>
#include <limits>
#include <sstream>
#include <vector>

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

class Exif_IO_EasyExif::EXIFInfoImpl {
  easyexif::EXIFInfo exif_info;
 public:
  easyexif::EXIFInfo& get() {return exif_info;}
};

Exif_IO_EasyExif::Exif_IO_EasyExif():
  bHaveExifInfo_( false ),
  pimpl_(new Exif_IO_EasyExif::EXIFInfoImpl())
{
}

Exif_IO_EasyExif::Exif_IO_EasyExif( const std::string & sFileName ):
  bHaveExifInfo_( false ),
  pimpl_(new Exif_IO_EasyExif::EXIFInfoImpl())
{
  open( sFileName );
}

bool Exif_IO_EasyExif::open( const std::string & sFileName )
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
  bHaveExifInfo_ = ( (*pimpl_).get().parseFrom( &buf[0], fsize ) == PARSE_EXIF_SUCCESS );

  return bHaveExifInfo_;
}

size_t Exif_IO_EasyExif::getWidth() const
{
  return (*pimpl_).get().ImageWidth;
}

size_t Exif_IO_EasyExif::getHeight() const
{
  return (*pimpl_).get().ImageHeight;
}

float Exif_IO_EasyExif::getFocal() const
{
  return static_cast<float>( (*pimpl_).get().FocalLength );
}

float Exif_IO_EasyExif::getFocalLengthIn35mm() const
{
  return static_cast<float>( (*pimpl_).get().FocalLengthIn35mm );
}

float Exif_IO_EasyExif::getFocalPlaneXResolution() const
{
  return static_cast<float>( (*pimpl_).get().LensInfo.FocalPlaneXResolution );
}

float Exif_IO_EasyExif::getFocalPlaneYResolution() const
{
  return static_cast<float>( (*pimpl_).get().LensInfo.FocalPlaneYResolution );
}

int Exif_IO_EasyExif::getFocalPlaneResolutionUnit() const
{
  return static_cast<int>( (*pimpl_).get().LensInfo.FocalPlaneResolutionUnit );
}

std::string Exif_IO_EasyExif::getBrand() const
{
  return trim_copy( (*pimpl_).get().Make );
}

std::string Exif_IO_EasyExif::getModel() const
{
  return trim_copy( (*pimpl_).get().Model );
}

std::string Exif_IO_EasyExif::getLensModel() const
{
  return trim_copy( (*pimpl_).get().LensInfo.Model );
}

std::string Exif_IO_EasyExif::getImageUniqueID() const
{
  return (*pimpl_).get().ImageUniqueID;
}

bool Exif_IO_EasyExif::doesHaveExifInfo() const
{
  return bHaveExifInfo_;
}

std::string Exif_IO_EasyExif::allExifData() const
{
  std::ostringstream os;
  os
      << "Camera make       : " << (*pimpl_).get().Make << "\n"
      << "Camera model      : " << (*pimpl_).get().Model << "\n"
      << "Software          : " << (*pimpl_).get().Software << "\n"
      << "Bits per sample   : " << (*pimpl_).get().BitsPerSample << "\n"
      << "Image width       : " << (*pimpl_).get().ImageWidth << "\n"
      << "Image height      : " << (*pimpl_).get().ImageHeight << "\n"
      << "Image description : " << (*pimpl_).get().ImageDescription << "\n"
      << "Image orientation : " << (*pimpl_).get().Orientation << "\n"
      << "Image copyright   : " << (*pimpl_).get().Copyright << "\n"
      << "Image date/time   : " << (*pimpl_).get().DateTime << "\n"
      << "Original date/time: " << (*pimpl_).get().DateTimeOriginal << "\n"
      << "Digitize date/time: " << (*pimpl_).get().DateTimeDigitized << "\n"
      << "Subsecond time    : " << (*pimpl_).get().SubSecTimeOriginal << "\n"
      << "Exposure time     : 1/" << ( unsigned ) ( 1.0 / (*pimpl_).get().ExposureTime ) << "\n"
      << "F-stop            : " << (*pimpl_).get().FNumber << "\n"
      << "ISO speed         : " << (*pimpl_).get().ISOSpeedRatings << "\n"
      << "Subject distance  : " << (*pimpl_).get().SubjectDistance << "\n"
      << "Exposure bias     : EV" << (*pimpl_).get().ExposureBiasValue << "\n"
      << "Flash used?       : " << (*pimpl_).get().Flash << "\n"
      << "Metering mode     : " << (*pimpl_).get().MeteringMode << "\n"
      << "Lens focal length : mm\n" << (*pimpl_).get().FocalLength << "\n"
      << "35mm focal length : mm\n" << (*pimpl_).get().FocalLengthIn35mm << "\n"
      << "GPS Latitude      : deg ( deg, min, sec )\n" << "("
      <<  (*pimpl_).get().GeoLocation.Latitude << ", "
      <<  (*pimpl_).get().GeoLocation.LatComponents.degrees << ", "
      <<  (*pimpl_).get().GeoLocation.LatComponents.minutes << ", "
      <<  (*pimpl_).get().GeoLocation.LatComponents.seconds << ", "
      <<  (*pimpl_).get().GeoLocation.LatComponents.direction << ")" << "\n"
      << "GPS Longitude      : deg ( deg, min, sec )\n" << "("
      <<  (*pimpl_).get().GeoLocation.Longitude << ", "
      <<  (*pimpl_).get().GeoLocation.LonComponents.degrees << ", "
      <<  (*pimpl_).get().GeoLocation.LonComponents.minutes << ", "
      <<  (*pimpl_).get().GeoLocation.LonComponents.seconds << ", "
      <<  (*pimpl_).get().GeoLocation.LonComponents.direction << ")" << "\n"
      << "GPS Altitude       : m" << (*pimpl_).get().GeoLocation.Altitude << "\n"
      << "Lens stop (min, max) : " << "("
      << (*pimpl_).get().LensInfo.FStopMin << ", "
      << (*pimpl_).get().LensInfo.FStopMax << ")"
      << "Lens focal (min, max) : " << "("
      << (*pimpl_).get().LensInfo.FocalLengthMin << ", "
      << (*pimpl_).get().LensInfo.FocalLengthMax << ")"
      << "Lens make : " <<  (*pimpl_).get().LensInfo.Make << "\n"
      << "Lens model : " << (*pimpl_).get().LensInfo.Model << "\n"
      << "Image Unique ID    : " << (*pimpl_).get().ImageUniqueID << "\n";
  return os.str();
}

bool Exif_IO_EasyExif::GPSLatitude(double * latitude) const
{
  if ((*pimpl_).get().GeoLocation.Latitude != std::numeric_limits<double>::infinity())
  {
    (*latitude) = (*pimpl_).get().GeoLocation.Latitude;
    return true;
  }
  return false;
}

bool Exif_IO_EasyExif::GPSLongitude(double * longitude) const
{
  if ((*pimpl_).get().GeoLocation.Longitude != std::numeric_limits<double>::infinity())
  {
    (*longitude) = (*pimpl_).get().GeoLocation.Longitude;
    return true;
  }
  return false;
}

bool Exif_IO_EasyExif::GPSAltitude(double * altitude) const
{
  if ((*pimpl_).get().GeoLocation.Altitude != std::numeric_limits<double>::infinity())
  {
    (*altitude) = (*pimpl_).get().GeoLocation.Altitude;
    return true;
  }
  return false;
}

} // namespace exif
} // namespace openMVG
