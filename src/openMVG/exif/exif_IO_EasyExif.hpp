#ifndef EXIF_IO_EASYEXIF_HPP
#define EXIF_IO_EASYEXIF_HPP

#include "openMVG/exif/exif_IO.hpp"
#include "third_party/easyexif/exif.h"

#include <fstream>
#include <sstream>
#include <vector>

namespace openMVG {
namespace exif  {

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
      std::string sbrand = exifInfo_.Make;
      // remove leading and trailing spaces
      sbrand.erase(0, sbrand.find_first_not_of(' '));
      sbrand.erase(sbrand.find_last_not_of(' '));
      return sbrand;
    }

    std::string getModel() const
    {
      std::string smodel = exifInfo_.Model;
      // remove leading and trailing spaces
      smodel.erase(0, smodel.find_first_not_of(' '));
      smodel.erase(smodel.find_last_not_of(' '));
      return smodel;
    }

    std::string getLensModel() const
    {
      return "";
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
        << "Camera make       : " << exifInfo_.Make
        << "Camera model      : " << exifInfo_.Model
        << "Software          : " << exifInfo_.Software
        << "Bits per sample   : " << exifInfo_.BitsPerSample
        << "Image width       : " << exifInfo_.ImageWidth
        << "Image height      : " << exifInfo_.ImageHeight
        << "Image description : " << exifInfo_.ImageDescription
        << "Image orientation : " << exifInfo_.Orientation
        << "Image copyright   : " << exifInfo_.Copyright
        << "Image date/time   : " << exifInfo_.DateTime
        << "Original date/time: " << exifInfo_.DateTimeOriginal
        << "Digitize date/time: " << exifInfo_.DateTimeDigitized
        << "Subsecond time    : " << exifInfo_.SubSecTimeOriginal
        << "Exposure time     : 1/time " << (unsigned) (1.0/exifInfo_.ExposureTime)
        << "F-stop            : " << exifInfo_.FNumber
        << "ISO speed         : " << exifInfo_.ISOSpeedRatings
        << "Subject distance  : " << exifInfo_.SubjectDistance
        << "Exposure bias     : EV" << exifInfo_.ExposureBiasValue
        << "Flash used?       : " << exifInfo_.Flash
        << "Metering mode     : " << exifInfo_.MeteringMode
        << "Lens focal length : mm\n" << exifInfo_.FocalLength
        << "35mm focal length : mm\n" << exifInfo_.FocalLengthIn35mm
        << "GPS Latitude      : deg ( deg, min, sec )\n" << "("
        <<  exifInfo_.GeoLocation.Latitude << ", " 
        <<  exifInfo_.GeoLocation.LatComponents.degrees << ", " 
        <<  exifInfo_.GeoLocation.LatComponents.minutes << ", " 
        <<  exifInfo_.GeoLocation.LatComponents.seconds << ", " 
        <<  exifInfo_.GeoLocation.LatComponents.direction << ")"
        << "GPS Longitude     : deg ( deg, min, sec )\n" << "("
        <<  exifInfo_.GeoLocation.Longitude << ", " 
        <<  exifInfo_.GeoLocation.LonComponents.degrees << ", " 
        <<  exifInfo_.GeoLocation.LonComponents.minutes << ", " 
        <<  exifInfo_.GeoLocation.LonComponents.seconds << ", " 
        <<  exifInfo_.GeoLocation.LonComponents.direction << ")"
        << "GPS Altitude      : m" << exifInfo_.GeoLocation.Altitude;
      return os.str();
    }

  private:
    EXIFInfo exifInfo_;
    bool bHaveExifInfo_;
};

} // namespace exif
} // namespace openMVG

#endif //EXIF_IO_EASYEXIF_HPP
