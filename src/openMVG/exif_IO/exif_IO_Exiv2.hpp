#ifndef EXIF_IO_EXIV2_HPP
#define EXIF_IO_EXIV2_HPP

#include "exif_IO.hpp"

#include <exiv2/exiv2.hpp>

#include <fstream>
#include <string>
#include <sstream>

class Exif_IO_Exiv2 : public Exif_IO
{
  public:
    Exif_IO_Exiv2()
    {
    }

    Exif_IO_Exiv2( const std::string & sFileName  )
    {
      open( sFileName );
    }

    size_t getWidth() const
    {
      std::ostringstream oss;
      oss << "Exif.Photo.PixelXDimension";
      Exiv2::ExifData &exifData = image->exifData();
      try
      {
        std::stringstream ss;
        ss << exifData[oss.str()].value().toString();
        size_t width;
        ss >> width;
        return width;
      }
      catch (Exiv2::AnyError& e)
      {
          std::cout << "Caught Exiv2 exception '" << e << "'\n";
          return -1;
      }
    }

    size_t getHeight() const
    {
      std::ostringstream oss;
      oss << "Exif.Photo.PixelYDimension";
      Exiv2::ExifData &exifData = image->exifData();
      try
      {
        std::stringstream ss;
        ss << exifData[oss.str()].value().toString();
        size_t height;
        ss >> height;
        return height;
      }
      catch (Exiv2::AnyError& e)
      {
          std::cout << "Caught Exiv2 exception '" << e << "'\n";
          return -1;
      }
    }

    float getFocal() const
    {
      std::ostringstream oss;
      oss << "Exif.Photo.FocalLength";
      Exiv2::ExifData &exifData = image->exifData();
      try
      {
        return exifData[oss.str()].value().toFloat();
      }
      catch (Exiv2::AnyError& e)
      {
          std::cout << "Caught Exiv2 exception '" << e << "'\n";
          return -1;
      }
    }

    std::string getBrand() const
    {
      std::ostringstream oss;
      oss << "Exif.Image.Make";
      Exiv2::ExifData &exifData = image->exifData();
      try
      {
        return exifData[oss.str()].value().toString();
      }
      catch (Exiv2::AnyError& e)
      {
          std::cout << "Caught Exiv2 exception '" << e << "'\n";
          return "";
      }
    }

    std::string getModel() const
    {
      std::ostringstream oss;
      oss << "Exif.Image.Model";
      Exiv2::ExifData &exifData = image->exifData();
      try
      {
        return exifData[oss.str()].value().toString();
      }
      catch (Exiv2::AnyError& e)
      {
          std::cout << "Caught Exiv2 exception '" << e << "'\n";
          return "";
      }
    }

    std::string getLensModel() const
    {
      std::ostringstream oss;
      oss << "Exif." << getBrand() << ".LensModel";
      Exiv2::ExifData &exifData = image->exifData();
      try
      {
        return exifData[oss.str()].value().toString();
      }
      catch (Exiv2::AnyError& e)
      {
          std::cout << "Caught Exiv2 exception '" << e << "'\n";
          return "";
      }
    }

    /** Open the file for checking and parsing */
    bool open( const std::string & sFileName )
    {
      bool isOpen = false;
      image = Exiv2::ImageFactory::open( sFileName.c_str() );
      Exiv2::IptcData &iptcData = image->iptcData();
      if (!iptcData.empty())
	isOpen = true;

      image->readMetadata();
      return isOpen;
    }

    /**Verify if the file has metadata*/
    bool doesHaveExifInfo() const
    {
      Exiv2::ExifData &exifData = image->exifData();
      return !exifData.empty();
    }
    /** Print all data*/
    std::string allExifData() const
    {
      try
      {
        Exiv2::ExifData &exifData = image->exifData();
        Exiv2::ExifData::const_iterator end = exifData.end();

        std::ostringstream oss;
        for (Exiv2::ExifData::const_iterator i = exifData.begin(); i != end; ++i)
        {
          const char* tn = i->typeName();
          oss << std::setw(44) << std::setfill(' ') << std::left
              << i->key() << " "
              << "0x" << std::setw(4) << std::setfill('0') << std::right
              << std::hex << i->tag() << " "
              << std::setw(9) << std::setfill(' ') << std::left
              << (tn ? tn : "Unknown") << " "
              << std::dec << std::setw(3)
              << std::setfill(' ') << std::right
              << i->count() << "  "
              << std::dec << i->value()
              << "\n";
        }
        return oss.str();
      }
      catch (Exiv2::AnyError& e)
      {
          std::cout << "Caught Exiv2 exception '" << e << "'\n";
          return "";
      }
    }

  private :
     Exiv2::Image::AutoPtr image;
};
#endif //EXIF_IO_EXIV2_HPP
