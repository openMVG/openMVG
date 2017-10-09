/**************************************************************************
  exif.cpp  -- A simple ISO C++ library to parse basic EXIF
               information from a JPEG file.

  Copyright (c) 2010-2016 Mayank Lahiri
  mlahiri@gmail.com
  All rights reserved (BSD License).

  See exif.h for version history.

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
#include "exif.h"

#include <algorithm>
#include <cstdint>
#include <stdio.h>
#include <limits>
#include <memory>
#include <vector>

using std::string;

namespace {

struct Rational {
  uint32_t numerator, denominator;
  operator double() const {
    if (denominator < 1e-20) {
      return 0;
    }
    return static_cast<double>(numerator) / static_cast<double>(denominator);
  }
};

// IF Entry
class IFEntry {
 public:
  using byte_vector = std::vector<uint8_t>;
  using ascii_vector = std::string;
  using short_vector = std::vector<uint16_t>;
  using long_vector = std::vector<uint32_t>;
  using rational_vector = std::vector<Rational>;

  IFEntry()
      : tag_(0xFF), format_(0xFF), data_(0), length_(0), val_byte_(nullptr) {}
  IFEntry(const IFEntry &) = delete;
  IFEntry &operator=(const IFEntry &) = delete;
  IFEntry(IFEntry &&other)
      : tag_(other.tag_),
        format_(other.format_),
        data_(other.data_),
        length_(other.length_),
        val_byte_(other.val_byte_) {
    other.tag_ = 0xFF;
    other.format_ = 0xFF;
    other.data_ = 0;
    other.length_ = 0;
    other.val_byte_ = nullptr;
  }
  ~IFEntry() { delete_union(); }
  unsigned short tag() const { return tag_; }
  void tag(unsigned short tag) { tag_ = tag; }
  unsigned short format() const { return format_; }
  bool format(unsigned short format) {
    switch (format) {
      case 0x01:
      case 0x02:
      case 0x03:
      case 0x04:
      case 0x05:
      case 0x07:
      case 0x09:
      case 0x0a:
      case 0xff:
        break;
      default:
        return false;
    }
    delete_union();
    format_ = format;
    new_union();
    return true;
  }
  unsigned data() const { return data_; }
  void data(unsigned data) { data_ = data; }
  unsigned length() const { return length_; }
  void length(unsigned length) { length_ = length; }

  // functions to access the data
  //
  // !! it's CALLER responsibility to check that format !!
  // !! is correct before accessing it's field          !!
  //
  // - getters are use here to allow future addition
  //   of checks if format is correct
  byte_vector &val_byte() { return *val_byte_; }
  ascii_vector &val_string() { return *val_string_; }
  short_vector &val_short() { return *val_short_; }
  long_vector &val_long() { return *val_long_; }
  rational_vector &val_rational() { return *val_rational_; }

 private:
  // Raw fields
  unsigned short tag_;
  unsigned short format_;
  unsigned data_;
  unsigned length_;

  // Parsed fields
  union {
    byte_vector *val_byte_;
    ascii_vector *val_string_;
    short_vector *val_short_;
    long_vector *val_long_;
    rational_vector *val_rational_;
  };

  void delete_union() {
    switch (format_) {
      case 0x1:
        delete val_byte_;
        val_byte_ = nullptr;
        break;
      case 0x2:
        delete val_string_;
        val_string_ = nullptr;
        break;
      case 0x3:
        delete val_short_;
        val_short_ = nullptr;
        break;
      case 0x4:
        delete val_long_;
        val_long_ = nullptr;
        break;
      case 0x5:
        delete val_rational_;
        val_rational_ = nullptr;
        break;
      case 0xff:
        break;
      default:
        // should not get here
        // should I throw an exception or ...?
        break;
    }
  }
  void new_union() {
    switch (format_) {
      case 0x1:
        val_byte_ = new byte_vector();
        break;
      case 0x2:
        val_string_ = new ascii_vector();
        break;
      case 0x3:
        val_short_ = new short_vector();
        break;
      case 0x4:
        val_long_ = new long_vector();
        break;
      case 0x5:
        val_rational_ = new rational_vector();
        break;
      case 0xff:
        break;
      default:
        // should not get here
        // should I throw an exception or ...?
        break;
    }
  }
};

// Helper functions
template <typename T, bool alignIntel>
T parse(const unsigned char *buf);

template <>
uint8_t parse<uint8_t, false>(const unsigned char *buf) {
  return *buf;
}

template <>
uint8_t parse<uint8_t, true>(const unsigned char *buf) {
  return *buf;
}

template <>
uint16_t parse<uint16_t, false>(const unsigned char *buf) {
  return (static_cast<uint16_t>(buf[0]) << 8) | buf[1];
}

template <>
uint16_t parse<uint16_t, true>(const unsigned char *buf) {
  return (static_cast<uint16_t>(buf[1]) << 8) | buf[0];
}

template <>
uint32_t parse<uint32_t, false>(const unsigned char *buf) {
  return (static_cast<uint32_t>(buf[0]) << 24) |
         (static_cast<uint32_t>(buf[1]) << 16) |
         (static_cast<uint32_t>(buf[2]) << 8) | buf[3];
}

template <>
uint32_t parse<uint32_t, true>(const unsigned char *buf) {
  return (static_cast<uint32_t>(buf[3]) << 24) |
         (static_cast<uint32_t>(buf[2]) << 16) |
         (static_cast<uint32_t>(buf[1]) << 8) | buf[0];
}

template <>
Rational parse<Rational, true>(const unsigned char *buf) {
  Rational r;
  r.numerator = parse<uint32_t, true>(buf);
  r.denominator = parse<uint32_t, true>(buf + 4);
  return r;
}

template <>
Rational parse<Rational, false>(const unsigned char *buf) {
  Rational r;
  r.numerator = parse<uint32_t, false>(buf);
  r.denominator = parse<uint32_t, false>(buf + 4);
  return r;
}

/**
 * Try to read entry.length() values for this entry.
 *
 * Returns:
 *  true  - entry.length() values were read
 *  false - something went wrong, vec's content was not touched
 */
template <typename T, bool alignIntel, typename C>
bool extract_values(C &container, const unsigned char *buf, const unsigned base,
                    const unsigned len, const IFEntry &entry) {
  const unsigned char *data;
  uint32_t reversed_data;
  // if data fits into 4 bytes, they are stored directly in
  // the data field in IFEntry
  if (sizeof(T) * entry.length() <= 4) {
    if (alignIntel) {
      reversed_data = entry.data();
    } else {
      reversed_data = entry.data();
      // this reversing works, but is ugly
      unsigned char *data = reinterpret_cast<unsigned char *>(&reversed_data);
      unsigned char tmp;
      tmp = data[0];
      data[0] = data[3];
      data[3] = tmp;
      tmp = data[1];
      data[1] = data[2];
      data[2] = tmp;
    }
    data = reinterpret_cast<const unsigned char *>(&(reversed_data));
  } else {
    data = buf + base + entry.data();
    if (data + sizeof(T) * entry.length() > buf + len) {
      return false;
    }
  }
  container.resize(entry.length());
  for (size_t i = 0; i < entry.length(); ++i) {
    container[i] = parse<T, alignIntel>(data + sizeof(T) * i);
  }
  return true;
}

template <bool alignIntel>
void parseIFEntryHeader(const unsigned char *buf, unsigned short &tag,
                        unsigned short &format, unsigned &length,
                        unsigned &data) {
  // Each directory entry is composed of:
  // 2 bytes: tag number (data field)
  // 2 bytes: data format
  // 4 bytes: number of components
  // 4 bytes: data value or offset to data value
  tag = parse<uint16_t, alignIntel>(buf);
  format = parse<uint16_t, alignIntel>(buf + 2);
  length = parse<uint32_t, alignIntel>(buf + 4);
  data = parse<uint32_t, alignIntel>(buf + 8);
}

template <bool alignIntel>
void parseIFEntryHeader(const unsigned char *buf, IFEntry &result) {
  unsigned short tag;
  unsigned short format;
  unsigned length;
  unsigned data;

  parseIFEntryHeader<alignIntel>(buf, tag, format, length, data);

  result.tag(tag);
  result.format(format);
  result.length(length);
  result.data(data);
}

template <bool alignIntel>
IFEntry parseIFEntry_temp(const unsigned char *buf, const unsigned offs,
                          const unsigned base, const unsigned len) {
  IFEntry result;

  // check if there even is enough data for IFEntry in the buffer
  if (buf + offs + 12 > buf + len) {
    result.tag(0xFF);
    return result;
  }

  parseIFEntryHeader<alignIntel>(buf + offs, result);

  // Parse value in specified format
  switch (result.format()) {
    case 1:
      if (!extract_values<uint8_t, alignIntel>(result.val_byte(), buf, base,
                                               len, result)) {
        result.tag(0xFF);
      }
      break;
    case 2:
      // string is basically sequence of uint8_t (well, according to EXIF even
      // uint7_t, but
      // we don't have that), so just read it as bytes
      if (!extract_values<uint8_t, alignIntel>(result.val_string(), buf, base,
                                               len, result)) {
        result.tag(0xFF);
      }
      // and cut zero byte at the end, since we don't want that in the
      // std::string
      if (!result.val_string().empty())  {
        result.val_string().resize(result.val_string().length() - 1);
      }
      break;
    case 3:
      if (!extract_values<uint16_t, alignIntel>(result.val_short(), buf, base,
                                                len, result)) {
        result.tag(0xFF);
      }
      break;
    case 4:
      if (!extract_values<uint32_t, alignIntel>(result.val_long(), buf, base,
                                                len, result)) {
        result.tag(0xFF);
      }
      break;
    case 5:
      if (!extract_values<Rational, alignIntel>(result.val_rational(), buf,
                                                base, len, result)) {
        result.tag(0xFF);
      }
      break;
    case 7:
    case 9:
    case 10:
      break;
    default:
      result.tag(0xFF);
  }
  return result;
}

// helper functions for convinience
template <typename T>
T parse_value(const unsigned char *buf, bool alignIntel) {
  if (alignIntel) {
    return parse<T, true>(buf);
  } else {
    return parse<T, false>(buf);
  }
}

void parseIFEntryHeader(const unsigned char *buf, bool alignIntel,
                        unsigned short &tag, unsigned short &format,
                        unsigned &length, unsigned &data) {
  if (alignIntel) {
    parseIFEntryHeader<true>(buf, tag, format, length, data);
  } else {
    parseIFEntryHeader<false>(buf, tag, format, length, data);
  }
}

IFEntry parseIFEntry(const unsigned char *buf, const unsigned offs,
                     const bool alignIntel, const unsigned base,
                     const unsigned len) {
  if (alignIntel) {
    return parseIFEntry_temp<true>(buf, offs, base, len);
  } else {
    return parseIFEntry_temp<false>(buf, offs, base, len);
  }
}
}

//
// Locates the EXIF segment and parses it using parseFromEXIFSegment
//
int easyexif::EXIFInfo::parseFrom(const unsigned char *buf, unsigned len) {
  // Sanity check: all JPEG files start with 0xFFD8.
  if (!buf || len < 4) return PARSE_EXIF_ERROR_NO_JPEG;
  if (buf[0] != 0xFF || buf[1] != 0xD8) return PARSE_EXIF_ERROR_NO_JPEG;

  // Sanity check: some cameras pad the JPEG image with some bytes at the end.
  // Normally, we should be able to find the JPEG end marker 0xFFD9 at the end
  // of the image buffer, but not always. As long as there are some bytes
  // except 0xD9 at the end of the image buffer, keep decrementing len until
  // an 0xFFD9 is found. If JPEG end marker 0xFFD9 is not found,
  // then we can be reasonably sure that the buffer is not a JPEG.
  while (len > 2) {
    if (buf[len - 1] == 0xD9 && buf[len - 2] == 0xFF)
      break;
    len--;
  }
  if (len <= 2)
    return PARSE_EXIF_ERROR_NO_JPEG;

  clear();

  // Scan for EXIF header (bytes 0xFF 0xE1) and do a sanity check by
  // looking for bytes "Exif\0\0". The marker length data is in Motorola
  // byte order, which results in the 'false' parameter to parse16().
  // The marker has to contain at least the TIFF header, otherwise the
  // EXIF data is corrupt. So the minimum length specified here has to be:
  //   2 bytes: section size
  //   6 bytes: "Exif\0\0" string
  //   2 bytes: TIFF header (either "II" or "MM" string)
  //   2 bytes: TIFF magic (short 0x2a00 in Motorola byte order)
  //   4 bytes: Offset to first IFD
  // =========
  //  16 bytes
  unsigned offs = 0;  // current offset into buffer
  for (offs = 0; offs < len - 1; offs++)
    if (buf[offs] == 0xFF && buf[offs + 1] == 0xE1) break;
  if (offs + 4 > len) return PARSE_EXIF_ERROR_NO_EXIF;
  offs += 2;
  unsigned short section_length = parse_value<uint16_t>(buf + offs, false);
  if (offs + section_length > len || section_length < 16)
    return PARSE_EXIF_ERROR_CORRUPT;
  offs += 2;

  return parseFromEXIFSegment(buf + offs, len - offs);
}

int easyexif::EXIFInfo::parseFrom(const string &data) {
  return parseFrom(
      reinterpret_cast<const unsigned char *>(data.data()), static_cast<unsigned>(data.length()));
}

//
// Main parsing function for an EXIF segment.
//
// PARAM: 'buf' start of the EXIF TIFF, which must be the bytes "Exif\0\0".
// PARAM: 'len' length of buffer
//
int easyexif::EXIFInfo::parseFromEXIFSegment(const unsigned char *buf,
                                             unsigned len) {
  bool alignIntel = true;  // byte alignment (defined in EXIF header)
  unsigned offs = 0;       // current offset into buffer
  if (!buf || len < 6) return PARSE_EXIF_ERROR_NO_EXIF;

  if (!std::equal(buf, buf + 6, "Exif\0\0")) return PARSE_EXIF_ERROR_NO_EXIF;
  offs += 6;

  // Now parsing the TIFF header. The first two bytes are either "II" or
  // "MM" for Intel or Motorola byte alignment. Sanity check by parsing
  // the unsigned short that follows, making sure it equals 0x2a. The
  // last 4 bytes are an offset into the first IFD, which are added to
  // the global offset counter. For this block, we expect the following
  // minimum size:
  //  2 bytes: 'II' or 'MM'
  //  2 bytes: 0x002a
  //  4 bytes: offset to first IDF
  // -----------------------------
  //  8 bytes
  if (offs + 8 > len) return PARSE_EXIF_ERROR_CORRUPT;
  unsigned tiff_header_start = offs;
  if (buf[offs] == 'I' && buf[offs + 1] == 'I')
    alignIntel = true;
  else {
    if (buf[offs] == 'M' && buf[offs + 1] == 'M')
      alignIntel = false;
    else
      return PARSE_EXIF_ERROR_UNKNOWN_BYTEALIGN;
  }
  this->ByteAlign = alignIntel;
  offs += 2;
  if (0x2a != parse_value<uint16_t>(buf + offs, alignIntel))
    return PARSE_EXIF_ERROR_CORRUPT;
  offs += 2;
  unsigned first_ifd_offset = parse_value<uint32_t>(buf + offs, alignIntel);
  offs += first_ifd_offset - 4;
  if (offs >= len) return PARSE_EXIF_ERROR_CORRUPT;

  // Now parsing the first Image File Directory (IFD0, for the main image).
  // An IFD consists of a variable number of 12-byte directory entries. The
  // first two bytes of the IFD section contain the number of directory
  // entries in the section. The last 4 bytes of the IFD contain an offset
  // to the next IFD, which means this IFD must contain exactly 6 + 12 * num
  // bytes of data.
  if (offs + 2 > len) return PARSE_EXIF_ERROR_CORRUPT;
  int num_entries = parse_value<uint16_t>(buf + offs, alignIntel);
  if (offs + 6 + 12 * num_entries > len) return PARSE_EXIF_ERROR_CORRUPT;
  offs += 2;
  unsigned exif_sub_ifd_offset = len;
  unsigned gps_sub_ifd_offset = len;
  while (--num_entries >= 0) {
    IFEntry result =
        parseIFEntry(buf, offs, alignIntel, tiff_header_start, len);
    offs += 12;
    switch (result.tag()) {
      case 0x102:
        // Bits per sample
        if (result.format() == 3 && result.val_short().size())
          this->BitsPerSample = result.val_short().front();
        break;

      case 0x10E:
        // Image description
        if (result.format() == 2) this->ImageDescription = result.val_string();
        break;

      case 0x10F:
        // Digicam make
        if (result.format() == 2) this->Make = result.val_string();
        break;

      case 0x110:
        // Digicam model
        if (result.format() == 2) this->Model = result.val_string();
        break;

      case 0x112:
        // Orientation of image
        if (result.format() == 3 && result.val_short().size())
          this->Orientation = result.val_short().front();
        break;

      case 0x131:
        // Software used for image
        if (result.format() == 2) this->Software = result.val_string();
        break;

      case 0x132:
        // EXIF/TIFF date/time of image modification
        if (result.format() == 2) this->DateTime = result.val_string();
        break;

      case 0x8298:
        // Copyright information
        if (result.format() == 2) this->Copyright = result.val_string();
        break;

      case 0x8825:
        // GPS IFS offset
        gps_sub_ifd_offset = tiff_header_start + result.data();
        break;

      case 0x8769:
        // EXIF SubIFD offset
        exif_sub_ifd_offset = tiff_header_start + result.data();
        break;
    }
  }

  // Jump to the EXIF SubIFD if it exists and parse all the information
  // there. Note that it's possible that the EXIF SubIFD doesn't exist.
  // The EXIF SubIFD contains most of the interesting information that a
  // typical user might want.
  if (exif_sub_ifd_offset + 4 <= len) {
    offs = exif_sub_ifd_offset;
    int num_entries = parse_value<uint16_t>(buf + offs, alignIntel);
    if (offs + 6 + 12 * num_entries > len) return PARSE_EXIF_ERROR_CORRUPT;
    offs += 2;
    while (--num_entries >= 0) {
      IFEntry result =
          parseIFEntry(buf, offs, alignIntel, tiff_header_start, len);
      switch (result.tag()) {
        case 0x829a:
          // Exposure time in seconds
          if (result.format() == 5 && result.val_rational().size())
            this->ExposureTime = result.val_rational().front();
          break;

        case 0x829d:
          // FNumber
          if (result.format() == 5 && result.val_rational().size())
            this->FNumber = result.val_rational().front();
          break;

      case 0x8822:
        // Exposure Program
        if (result.format() == 3 && result.val_short().size())
          this->ExposureProgram = result.val_short().front();
        break;

        case 0x8827:
          // ISO Speed Rating
          if (result.format() == 3 && result.val_short().size())
            this->ISOSpeedRatings = result.val_short().front();
          break;

        case 0x9003:
          // Original date and time
          if (result.format() == 2)
            this->DateTimeOriginal = result.val_string();
          break;

        case 0x9004:
          // Digitization date and time
          if (result.format() == 2)
            this->DateTimeDigitized = result.val_string();
          break;

        case 0x9201:
          // Shutter speed value
          if (result.format() == 5 && result.val_rational().size())
            this->ShutterSpeedValue = result.val_rational().front();
          break;

        case 0x9204:
          // Exposure bias value
          if (result.format() == 5 && result.val_rational().size())
            this->ExposureBiasValue = result.val_rational().front();
          break;

        case 0x9206:
          // Subject distance
          if (result.format() == 5 && result.val_rational().size())
            this->SubjectDistance = result.val_rational().front();
          break;

        case 0x9209:
          // Flash used
          if (result.format() == 3 && result.val_short().size()) {
            uint16_t data = result.val_short().front();
            
            this->Flash = data & 1;
            this->FlashReturnedLight = (data & 6) >> 1;
            this->FlashMode = (data & 24) >> 3;
          }
          break;

        case 0x920a:
          // Focal length
          if (result.format() == 5 && result.val_rational().size())
            this->FocalLength = result.val_rational().front();
          break;

        case 0x9207:
          // Metering mode
          if (result.format() == 3 && result.val_short().size())
            this->MeteringMode = result.val_short().front();
          break;

        case 0x9291:
          // Subsecond original time
          if (result.format() == 2)
            this->SubSecTimeOriginal = result.val_string();
          break;

        case 0xa002:
          // EXIF Image width
          if (result.format() == 4 && result.val_long().size())
            this->ImageWidth = result.val_long().front();
          if (result.format() == 3 && result.val_short().size())
            this->ImageWidth = result.val_short().front();
          break;

        case 0xa003:
          // EXIF Image height
          if (result.format() == 4 && result.val_long().size())
            this->ImageHeight = result.val_long().front();
          if (result.format() == 3 && result.val_short().size())
            this->ImageHeight = result.val_short().front();
          break;

        case 0xa20e:
          // EXIF Focal plane X-resolution
          if (result.format() == 5) {
            this->LensInfo.FocalPlaneXResolution = result.val_rational()[0];
          }
          break;

        case 0xa20f:
          // EXIF Focal plane Y-resolution
          if (result.format() == 5) {
            this->LensInfo.FocalPlaneYResolution = result.val_rational()[0];
          }
          break;

          case 0xa210:
            // EXIF Focal plane resolution unit
            if (result.format() == 3 && result.val_short().size()) {
               this->LensInfo.FocalPlaneResolutionUnit = result.val_short().front();
            }
          break;

        case 0xa405:
          // Focal length in 35mm film
          if (result.format() == 3 && result.val_short().size())
            this->FocalLengthIn35mm = result.val_short().front();
          break;

        case 0xa432:
          // Focal length and FStop.
          if (result.format() == 5) {
            int sz = static_cast<unsigned>(result.val_rational().size());
            if (sz)
              this->LensInfo.FocalLengthMin = result.val_rational()[0];
            if (sz > 1)
              this->LensInfo.FocalLengthMax = result.val_rational()[1];
            if (sz > 2)
              this->LensInfo.FStopMin = result.val_rational()[2];
            if (sz > 3)
              this->LensInfo.FStopMax = result.val_rational()[3];
          }
          break;

        case 0xa433:
          // Lens make.
          if (result.format() == 2) {
            this->LensInfo.Make = result.val_string();
          }
          break;

        case 0xa434:
          // Lens model.
          if (result.format() == 2) {
            this->LensInfo.Model = result.val_string();
          }
          break;
      }
      offs += 12;
    }
  }

  // Jump to the GPS SubIFD if it exists and parse all the information
  // there. Note that it's possible that the GPS SubIFD doesn't exist.
  if (gps_sub_ifd_offset + 4 <= len) {
    offs = gps_sub_ifd_offset;
    int num_entries = parse_value<uint16_t>(buf + offs, alignIntel);
    if (offs + 6 + 12 * num_entries > len) return PARSE_EXIF_ERROR_CORRUPT;
    offs += 2;
    while (--num_entries >= 0) {
      unsigned short tag, format;
      unsigned length, data;
      parseIFEntryHeader(buf + offs, alignIntel, tag, format, length, data);
      switch (tag) {
        case 1:
          // GPS north or south
          this->GeoLocation.LatComponents.direction = *(buf + offs + 8);
          if (this->GeoLocation.LatComponents.direction == 0) {
            this->GeoLocation.LatComponents.direction = '?';
          }
          if ('S' == this->GeoLocation.LatComponents.direction) {
            this->GeoLocation.Latitude = -this->GeoLocation.Latitude;
          }
          break;

        case 2:
          // GPS latitude
          if ((format == 5 || format == 10) && length == 3) {
            this->GeoLocation.LatComponents.degrees = parse_value<Rational>(
                buf + data + tiff_header_start, alignIntel);
            this->GeoLocation.LatComponents.minutes = parse_value<Rational>(
                buf + data + tiff_header_start + 8, alignIntel);
            this->GeoLocation.LatComponents.seconds = parse_value<Rational>(
                buf + data + tiff_header_start + 16, alignIntel);
            this->GeoLocation.Latitude =
                this->GeoLocation.LatComponents.degrees +
                this->GeoLocation.LatComponents.minutes / 60 +
                this->GeoLocation.LatComponents.seconds / 3600;
            if ('S' == this->GeoLocation.LatComponents.direction) {
              this->GeoLocation.Latitude = -this->GeoLocation.Latitude;
            }
          }
          break;

        case 3:
          // GPS east or west
          this->GeoLocation.LonComponents.direction = *(buf + offs + 8);
          if (this->GeoLocation.LonComponents.direction == 0) {
            this->GeoLocation.LonComponents.direction = '?';
          }
          if ('W' == this->GeoLocation.LonComponents.direction) {
            this->GeoLocation.Longitude = -this->GeoLocation.Longitude;
          }
          break;

        case 4:
          // GPS longitude
          if ((format == 5 || format == 10) && length == 3) {
            this->GeoLocation.LonComponents.degrees = parse_value<Rational>(
                buf + data + tiff_header_start, alignIntel);
            this->GeoLocation.LonComponents.minutes = parse_value<Rational>(
                buf + data + tiff_header_start + 8, alignIntel);
            this->GeoLocation.LonComponents.seconds = parse_value<Rational>(
                buf + data + tiff_header_start + 16, alignIntel);
            this->GeoLocation.Longitude =
                this->GeoLocation.LonComponents.degrees +
                this->GeoLocation.LonComponents.minutes / 60 +
                this->GeoLocation.LonComponents.seconds / 3600;
            if ('W' == this->GeoLocation.LonComponents.direction)
              this->GeoLocation.Longitude = -this->GeoLocation.Longitude;
          }
          break;

        case 5:
          // GPS altitude reference (below or above sea level)
          this->GeoLocation.AltitudeRef = *(buf + offs + 8);
          if (1 == this->GeoLocation.AltitudeRef) {
            this->GeoLocation.Altitude = -this->GeoLocation.Altitude;
          }
          break;

        case 6:
          // GPS altitude
          if ((format == 5 || format == 10)) {
            this->GeoLocation.Altitude = parse_value<Rational>(
                buf + data + tiff_header_start, alignIntel);
            if (1 == this->GeoLocation.AltitudeRef) {
              this->GeoLocation.Altitude = -this->GeoLocation.Altitude;
            }
          }
          break;

        case 11:
          // GPS degree of precision (DOP)
          if ((format == 5 || format == 10)) {
            this->GeoLocation.DOP = parse_value<Rational>(
                buf + data + tiff_header_start, alignIntel);
          }
          break;
      }
      offs += 12;
    }
  }

  return PARSE_EXIF_SUCCESS;
}

void easyexif::EXIFInfo::clear() {
  // Strings
  ImageDescription = "";
  Make = "";
  Model = "";
  Software = "";
  DateTime = "";
  DateTimeOriginal = "";
  DateTimeDigitized = "";
  SubSecTimeOriginal = "";
  Copyright = "";

  // Shorts / unsigned / double
  ByteAlign = 0;
  Orientation = 0;

  BitsPerSample = 0;
  ExposureTime = 0;
  FNumber = 0;
  ExposureProgram = 0;
  ISOSpeedRatings = 0;
  ShutterSpeedValue = 0;
  ExposureBiasValue = 0;
  SubjectDistance = 0;
  FocalLength = 0;
  FocalLengthIn35mm = 0;
  Flash = 0;
  FlashReturnedLight = 0;
  FlashMode = 0;
  MeteringMode = 0;
  ImageWidth = 0;
  ImageHeight = 0;

  // Geolocation
  GeoLocation.Latitude = std::numeric_limits<double>::infinity();
  GeoLocation.Longitude = std::numeric_limits<double>::infinity();
  GeoLocation.Altitude = std::numeric_limits<double>::infinity();
  GeoLocation.AltitudeRef = 0;
  GeoLocation.DOP = std::numeric_limits<double>::infinity();
  GeoLocation.LatComponents.degrees = std::numeric_limits<double>::infinity();
  GeoLocation.LatComponents.minutes = std::numeric_limits<double>::infinity();
  GeoLocation.LatComponents.seconds = std::numeric_limits<double>::infinity();
  GeoLocation.LatComponents.direction = '?';
  GeoLocation.LonComponents.degrees = std::numeric_limits<double>::infinity();
  GeoLocation.LonComponents.minutes = std::numeric_limits<double>::infinity();
  GeoLocation.LonComponents.seconds = std::numeric_limits<double>::infinity();
  GeoLocation.LonComponents.direction = '?';

  // LensInfo
  LensInfo.FocalLengthMax = 0;
  LensInfo.FocalLengthMin = 0;
  LensInfo.FStopMax = 0;
  LensInfo.FStopMin = 0;
  LensInfo.FocalPlaneYResolution = 0;
  LensInfo.FocalPlaneXResolution = 0;
  LensInfo.FocalPlaneResolutionUnit= 0;
  LensInfo.Make = "";
  LensInfo.Model = "";
}
