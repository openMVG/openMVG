// Copyright (c) 2013-2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef EXIF_IO_EASYEXIF_HPP
#define EXIF_IO_EASYEXIF_HPP

#include "openMVG/exif/exif_IO.hpp"
#include "third_party/easyexif/exif.h"

#include <fstream>
#include <sstream>
#include <vector>

namespace openMVG {
namespace exif  {

/// Remove all leading and trailing spaces from the input. The result is a trimmed copy of the input
inline std::string trim_copy(const std::string& s)
{
  if(s.empty())
    return s;

  std::string res(s);
  // remove leading and trailing spaces
  res.erase(0, res.find_first_not_of(' '));
  res.erase(res.find_last_not_of(' ')+1);
  // handle multiple trailing end character
  res = res.substr(0, res.find('\0'));
  return res;
}

class Exif_IO_EasyExif : public Exif_IO
{
  public:
    Exif_IO_EasyExif(): bHaveExifInfo_(false)
    {
    }

    Exif_IO_EasyExif( const std::string & sFileName ): bHaveExifInfo_(false)
    {
      open(sFileName);
    }

    bool open( const std::string & sFileName )
    {
      // Read the file into a buffer
      FILE *fp = fopen(sFileName.c_str(), "rb");
      if (!fp) {
        return false;
      }
      fseek(fp, 0, SEEK_END);
      unsigned long fsize = ftell(fp);
      rewind(fp);
      std::vector<unsigned char> buf(fsize);
      if (fread(&buf[0], 1, fsize, fp) != fsize) {
        fclose(fp);
        return false;
      }
      fclose(fp);

      // Parse EXIF
      bHaveExifInfo_ = (exifInfo_.parseFrom(&buf[0], fsize) == PARSE_EXIF_SUCCESS);

      return bHaveExifInfo_;
    }

    size_t getWidth() const
    {
      return exifInfo_.ImageWidth;
    }

    size_t getHeight() const
    {
      return exifInfo_.ImageHeight;
    }

    float getFocal() const
    {
      return static_cast<float>(exifInfo_.FocalLength);
    }

    std::string getBrand() const
    {
      return trim_copy(exifInfo_.Make);
    }

    std::string getModel() const
    {
      return trim_copy(exifInfo_.Model);
    }

    std::string getLensModel() const
    {
      return trim_copy(exifInfo_.LensInfo.Model);
    }

     std::string getImageUniqueID() const
    {
      return exifInfo_.ImageUniqueID;
    }

    /**Verify if the file has metadata*/
    bool doesHaveExifInfo() const
    {
      return bHaveExifInfo_;
    }

    /** Print all data*/
    std::string allExifData() const
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
        << "Exposure time     : 1/" << (unsigned) (1.0/exifInfo_.ExposureTime) << "\n"
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

  private:
    easyexif::EXIFInfo exifInfo_;
    bool bHaveExifInfo_;
};

} // namespace exif
} // namespace openMVG

#endif //EXIF_IO_EASYEXIF_HPP
