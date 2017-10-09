/**************************************************************************
  exif.h  -- A simple ISO C++ library to parse basic EXIF
             information from a JPEG file.

  Based on the description of the EXIF file format at:
  -- http://park2.wakwak.com/~tsuruzoh/Computer/Digicams/exif-e.html
  -- http://www.media.mit.edu/pia/Research/deepview/exif.html
  -- http://www.exif.org/Exif2-2.PDF

  Copyright (c) 2010-2015 Mayank Lahiri
  mlahiri@gmail.com
  All rights reserved.

  VERSION HISTORY:
  ================

  2.2: Release December 2014
       --

  2.1: Released July 2013
       -- fixed a bug where JPEGs without an EXIF SubIFD would not be parsed
       -- fixed a bug in parsing GPS coordinate seconds
       -- fixed makefile bug
       -- added two pathological test images from Matt Galloway
          http://www.galloway.me.uk/2012/01/uiimageorientation-exif-orientation-sample-images/
       -- split main parsing routine for easier integration into Firefox

  2.0: Released February 2013
       -- complete rewrite
       -- no new/delete
       -- added GPS support

  1.0: Released 2010

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are met:

  -- Redistributions of source code must retain the above copyright notice,
     this list of conditions and the following disclaimer.
  -- Redistributions in binary form must reproduce the above copyright notice,
     this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ``AS IS'' AND ANY EXPRESS
   OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
   OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN
   NO EVENT SHALL THE FREEBSD PROJECT OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
   INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
   BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
   OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
   EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
#ifndef __EXIF_H
#define __EXIF_H

#include <string>

namespace easyexif {

//
// Class responsible for storing and parsing EXIF information from a JPEG blob
//
class EXIFInfo {
 public:
  // Parsing function for an entire JPEG image buffer.
  //
  // PARAM 'data': A pointer to a JPEG image.
  // PARAM 'length': The length of the JPEG image.
  // RETURN:  PARSE_EXIF_SUCCESS (0) on success with 'result' filled out
  //          error code otherwise, as defined by the PARSE_EXIF_ERROR_* macros
  int parseFrom(const unsigned char *data, unsigned length);
  int parseFrom(const std::string &data);

  // Parsing function for an EXIF segment. This is used internally by parseFrom()
  // but can be called for special cases where only the EXIF section is
  // available (i.e., a blob starting with the bytes "Exif\0\0").
  int parseFromEXIFSegment(const unsigned char *buf, unsigned len);

  // Set all data members to default values.
  void clear();

  // Data fields filled out by parseFrom()
  char ByteAlign;                   // 0 = Motorola byte alignment, 1 = Intel
  std::string ImageDescription;     // Image description
  std::string Make;                 // Camera manufacturer's name
  std::string Model;                // Camera model
  unsigned short Orientation;       // Image orientation, start of data corresponds to
                                    // 0: unspecified in EXIF data
                                    // 1: upper left of image
                                    // 3: lower right of image
                                    // 6: upper right of image
                                    // 8: lower left of image
                                    // 9: undefined
  unsigned short BitsPerSample;     // Number of bits per component
  std::string Software;             // Software used
  std::string DateTime;             // File change date and time
  std::string DateTimeOriginal;     // Original file date and time (may not exist)
  std::string DateTimeDigitized;    // Digitization date and time (may not exist)
  std::string SubSecTimeOriginal;   // Sub-second time that original picture was taken
  std::string Copyright;            // File copyright information
  std::string ImageUniqueID;        // Unique identifier assigned to the picture
  double ExposureTime;              // Exposure time in seconds
  double FNumber;                   // F/stop
  unsigned short ExposureProgram;   // Exposure program
                                    // 0: Not defined
                                    // 1: Manual
                                    // 2: Normal program
                                    // 3: Aperture priority
                                    // 4: Shutter priority
                                    // 5: Creative program
                                    // 6: Action program
                                    // 7: Portrait mode
                                    // 8: Landscape mode
  unsigned short ISOSpeedRatings;   // ISO speed
  double ShutterSpeedValue;         // Shutter speed (reciprocal of exposure time)
  double ExposureBiasValue;         // Exposure bias value in EV
  double SubjectDistance;           // Distance to focus point in meters
  double FocalLength;               // Focal length of lens in millimeters
  unsigned short FocalLengthIn35mm; // Focal length in 35mm film
  char Flash;                       // 0 = no flash, 1 = flash used
  unsigned short FlashReturnedLight;// Flash returned light status
                                    // 0: No strobe return detection function
                                    // 1: Reserved
                                    // 2: Strobe return light not detected
                                    // 3: Strobe return light detected
  unsigned short FlashMode;         // Flash mode
                                    // 0: Unknown
                                    // 1: Compulsory flash firing
                                    // 2: Compulsory flash suppression
                                    // 3: Automatic mode
  unsigned short MeteringMode;      // Metering mode
                                    // 1: average
                                    // 2: center weighted average
                                    // 3: spot
                                    // 4: multi-spot
                                    // 5: multi-segment
  unsigned ImageWidth;              // Image width reported in EXIF data
  unsigned ImageHeight;             // Image height reported in EXIF data
  struct Geolocation_t {            // GPS information embedded in file
    double Latitude;                  // Image latitude expressed as decimal
    double Longitude;                 // Image longitude expressed as decimal
    double Altitude;                  // Altitude in meters, relative to sea level
    char AltitudeRef;                 // 0 = above sea level, -1 = below sea level
    double DOP;                       // GPS degree of precision (DOP)
    struct Coord_t {
      double degrees;
      double minutes;
      double seconds;
      char direction;
    } LatComponents, LonComponents;   // Latitude, Longitude expressed in deg/min/sec
  } GeoLocation;
  struct LensInfo_t {               // Lens information
    double FStopMin;                // Min aperture (f-stop)
    double FStopMax;                // Max aperture (f-stop)
    double FocalLengthMin;          // Min focal length (mm)
    double FocalLengthMax;          // Max focal length (mm)
    double FocalPlaneXResolution;   // Focal plane X-resolution
    double FocalPlaneYResolution;   // Focal plane Y-resolution
    unsigned short FocalPlaneResolutionUnit; // Focal plane resolution unit
                                             // 1: No absolute unit of measurement.
                                             // 2: Inch.
                                             // 3: Centimeter.
                                             // 4: Millimeter.
                                             // 5: Micrometer.
    std::string Make;               // Lens manufacturer
    std::string Model;              // Lens model
  } LensInfo;


  EXIFInfo() {
    clear();
  }
};

}

// Parse was successful
#define PARSE_EXIF_SUCCESS                    0
// No JPEG markers found in buffer, possibly invalid JPEG file
#define PARSE_EXIF_ERROR_NO_JPEG              1982
// No EXIF header found in JPEG file.
#define PARSE_EXIF_ERROR_NO_EXIF              1983
// Byte alignment specified in EXIF file was unknown (not Motorola or Intel).
#define PARSE_EXIF_ERROR_UNKNOWN_BYTEALIGN    1984
// EXIF header was found, but data was corrupted.
#define PARSE_EXIF_ERROR_CORRUPT              1985

#endif
