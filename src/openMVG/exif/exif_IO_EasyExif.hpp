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
        return false;
      }
      fclose(fp);

      // Parse EXIF
      int code = exifInfo_.parseFrom(&buf[0], fsize);
      if (code)
        bHaveExifInfo_ = false;
      else
        bHaveExifInfo_ = true;

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

    std::string getImageUniqueID() const
    {
      return trim_copy(exifInfo_.ImageUniqueID);
    }

    std::string getSerialNumber() const
    {
      return trim_copy(exifInfo_.SerialNumber);
    }

    std::string getLensModel() const
    {
      return trim_copy(exifInfo_.LensModel);
    }

    std::string getLensSerialNumber() const
    {
      return trim_copy(exifInfo_.LensSerialNumber);
    }

    // File change date and time
    std::string getDateTime() const
    {
      return trim_copy(exifInfo_.DateTime);
    }

    // Original file date and time (may not exist)
    std::string getDateTimeOriginal() const
    {
      return trim_copy(exifInfo_.DateTimeOriginal);
    }

    // Digitization date and time (may not exist)
    std::string getDateTimeDigitized() const
    {
      return trim_copy(exifInfo_.DateTimeDigitized);
    }

    // Sub-second time that original picture was taken
    std::string getSubSecTimeOriginal() const
    {
      return trim_copy(exifInfo_.SubSecTimeOriginal);
    }

    /**Verify if the file has metadata*/
    bool doesHaveExifInfo() const
    {
      return bHaveExifInfo_;
    }

    /** Print all data*/
    std::string getExifDataString() const
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
        << "Lens make : " <<  exifInfo_.LensInfo.Make
        << "Lens model : " << exifInfo_.LensInfo.Model
        << "Image Unique ID    : " << exifInfo_.ImageUniqueID << "\n"
        << "Serial Number      : " << exifInfo_.SerialNumber << "\n"
        << "Lens Serial Number : " << exifInfo_.LensSerialNumber << "\n";
      return os.str();
    }

    std::map<std::string, std::string> getExifData() const
    {
      std::map<std::string, std::string> allExifData;
      if(!exifInfo_.Make.empty())
        allExifData.emplace("camera_make", exifInfo_.Make);
      if(!exifInfo_.Model.empty())
        allExifData.emplace("camera_model", exifInfo_.Model);
      if(!exifInfo_.Software.empty())
        allExifData.emplace("software", exifInfo_.Software);
      if(exifInfo_.BitsPerSample > 0)
        allExifData.emplace("bits_per_sample", std::to_string(exifInfo_.BitsPerSample));
      if(exifInfo_.ImageWidth > 0)
        allExifData.emplace("image_width", std::to_string(exifInfo_.ImageWidth));
      if(exifInfo_.ImageHeight > 0)
        allExifData.emplace("image_height", std::to_string(exifInfo_.ImageHeight));
      if(!exifInfo_.ImageDescription.empty())
        allExifData.emplace("image_description", exifInfo_.ImageDescription);
      if(exifInfo_.Orientation > 0)
        allExifData.emplace("image_orientation", std::to_string(exifInfo_.Orientation));
      if(!exifInfo_.Copyright.empty())
        allExifData.emplace("image_copyright", exifInfo_.Copyright);
      if(!exifInfo_.DateTime.empty())
        allExifData.emplace("image_date_time", exifInfo_.DateTime);
      if(!exifInfo_.DateTimeOriginal.empty())
        allExifData.emplace("original_date_time", exifInfo_.DateTimeOriginal);
      if(!exifInfo_.DateTimeDigitized.empty())
        allExifData.emplace("digitize_date_time", exifInfo_.DateTimeDigitized);
      if(!exifInfo_.SubSecTimeOriginal.empty())
        allExifData.emplace("subsecond_time", exifInfo_.SubSecTimeOriginal);
      if(exifInfo_.ExposureTime > 0)
        allExifData.emplace("exposure_time", std::to_string(exifInfo_.ExposureTime));
      if(exifInfo_.FNumber > 0)
        allExifData.emplace("f_stop", std::to_string(exifInfo_.FNumber));
      if(exifInfo_.ISOSpeedRatings > 0)
        allExifData.emplace("iso_speed", std::to_string(exifInfo_.ISOSpeedRatings));
      if(exifInfo_.SubjectDistance > 0)
        allExifData.emplace("subject_distance", std::to_string(exifInfo_.SubjectDistance));
      if(exifInfo_.ExposureBiasValue > 0)
        allExifData.emplace("exposure_bias", std::to_string(exifInfo_.ExposureBiasValue));
      if(exifInfo_.Flash > 0)
        allExifData.emplace("flash_used", std::to_string(exifInfo_.Flash));
      if(exifInfo_.MeteringMode > 0)
        allExifData.emplace("metering_mode", std::to_string(exifInfo_.MeteringMode));
      if(exifInfo_.FocalLength > 0)
        allExifData.emplace("lens_focal_length", std::to_string(exifInfo_.FocalLength));
      if(exifInfo_.FocalLengthIn35mm > 0)
        allExifData.emplace("35mm_focal_length", std::to_string(exifInfo_.FocalLengthIn35mm));
      if(exifInfo_.GeoLocation.Latitude > 0
        || exifInfo_.GeoLocation.LatComponents.degrees > 0
        || exifInfo_.GeoLocation.LatComponents.minutes > 0
        || exifInfo_.GeoLocation.LatComponents.seconds > 0
        || exifInfo_.GeoLocation.LatComponents.direction > 0)
      {
        allExifData.emplace("gps_latitude", "("
          + std::to_string(exifInfo_.GeoLocation.Latitude) + ", "
          + std::to_string(exifInfo_.GeoLocation.LatComponents.degrees) + ", "
          + std::to_string(exifInfo_.GeoLocation.LatComponents.minutes) + ", "
          + std::to_string(exifInfo_.GeoLocation.LatComponents.seconds) + ", "
          + std::to_string(exifInfo_.GeoLocation.LatComponents.direction)  + ")");
      }
      if(exifInfo_.GeoLocation.Longitude > 0
        || exifInfo_.GeoLocation.LonComponents.degrees > 0
        || exifInfo_.GeoLocation.LonComponents.minutes > 0
        || exifInfo_.GeoLocation.LonComponents.seconds > 0
        || exifInfo_.GeoLocation.LonComponents.direction > 0)
      {
        allExifData.emplace("gps_longitude", "("
          + std::to_string(exifInfo_.GeoLocation.Longitude) + ", "
          + std::to_string(exifInfo_.GeoLocation.LonComponents.degrees) + ", "
          + std::to_string(exifInfo_.GeoLocation.LonComponents.minutes) + ", "
          + std::to_string(exifInfo_.GeoLocation.LonComponents.seconds) + ", "
          + std::to_string(exifInfo_.GeoLocation.LonComponents.direction)  + ")");
      }
      if(exifInfo_.GeoLocation.Altitude > 0)
        allExifData.emplace("gps_altitude", std::to_string(exifInfo_.GeoLocation.Altitude));
      if(exifInfo_.LensInfo.FStopMin > 0)
        allExifData.emplace("lens_stop_min", std::to_string(exifInfo_.LensInfo.FStopMin));
      if(exifInfo_.LensInfo.FStopMax > 0)
        allExifData.emplace("lens_stop_max", std::to_string(exifInfo_.LensInfo.FStopMax));
      if(exifInfo_.LensInfo.FocalLengthMin > 0)
        allExifData.emplace("lens_focal_min", std::to_string(exifInfo_.LensInfo.FocalLengthMin));
      if(exifInfo_.LensInfo.FocalLengthMax > 0)
        allExifData.emplace("lens focal_max", std::to_string(exifInfo_.LensInfo.FocalLengthMax));
      if(!exifInfo_.LensInfo.Make.empty())
        allExifData.emplace("lens_make", exifInfo_.LensInfo.Make);
      if(!exifInfo_.LensInfo.Model.empty())
        allExifData.emplace("lens_model", exifInfo_.LensInfo.Model);

      return allExifData;
    }

  private:
    EXIFInfo exifInfo_;
    bool bHaveExifInfo_;
};

} // namespace exif
} // namespace openMVG

#endif //EXIF_IO_EASYEXIF_HPP
