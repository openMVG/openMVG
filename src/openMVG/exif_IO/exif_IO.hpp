#ifndef EXIF_IO_HPP
#define EXIF_IO_HPP

#include <string>

class Exif_IO
{
  public:
    virtual size_t getWidth() const = 0;

    virtual size_t getHeight() const = 0;

    virtual float getFocal() const = 0;

    virtual std::string getBrand() const = 0;

    virtual std::string getModel() const = 0;

    virtual std::string getLensModel() const = 0;

    /** Open the file for checking and parsing */
    virtual bool open( const std::string & sFileName ) = 0;

    /**Verify if the file has metadata*/
    virtual bool doesHaveExifInfo() const = 0;

    /** Print all data*/
    virtual std::string allExifData() const = 0;

};
#endif //EXIF_IO_HPP

