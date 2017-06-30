// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2013-2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_EXIF_EXIF_IO_HPP
#define OPENMVG_EXIF_EXIF_IO_HPP

#include <string>

namespace openMVG
{
/**
* @brief namespace containing various classes and functions used to manage EXIF data
*/
namespace exif
{

/**
* @brief Standard class for EXchangeable Image file Format (EXIF) manipulation
*/
class Exif_IO
{
  public:

    /**
    * @brief Get image width
    * @return Width of the image (in pixel)
    */
    virtual size_t getWidth() const = 0;

    /**
    * @brief Get image height
    * @return Height of the image (in pixel)
    */
    virtual size_t getHeight() const = 0;

    /**
    * @brief Get Focal (in mm)
    * @return Focal of the lens when image was shot (in mm)
    */
    virtual float getFocal() const = 0;

    /**
    * @brief Get Brand of the camera
    * @return Brand name
    */
    virtual std::string getBrand() const = 0;

    /**
    * @brief Get Model of the camera
    * @return Camera model name
    */
    virtual std::string getModel() const = 0;

    /**
    * @brief Get Lens model
    * @return Lens model name
    */
    virtual std::string getLensModel() const = 0;

    /**
    * @brief Get a unique identifier for this image
    * @return Unique ID
    */
    virtual std::string getImageUniqueID() const = 0;

    /**
    * @brief Open the file for checking and parsing
    * @param sFileName path of the file to open
    * @retval true If file could be opened without error
    * @retval false If there was an error during opening
    */
    virtual bool open( const std::string & sFileName ) = 0;

    /**
    * @brief Verify if the file has metadata
    * @retval true If EXIF value are present
    * @retval false If EXIF value are not present
    */
    virtual bool doesHaveExifInfo() const = 0;

    /**
    * @brief Print all data
    * @return string containing all EXIF data
    */
    virtual std::string allExifData() const = 0;

    /**
    * @brief Try to read and save the EXIF GPS latitude
    * @return If GPS Latitude can be read & exported, return true
    */
    virtual bool GPSLatitude(double * latitude) const = 0;

    /**
    * @brief Try to read and save the EXIF GPS longitude
    * @return If GPS Longitude can be read & exported, return true
    */
    virtual bool GPSLongitude(double * longitude) const = 0;

   /**
    * @brief Try to read and save the EXIF GPS altitude
    * @return If GPS Altitude can be read & exported, return true
    */
    virtual bool GPSAltitude(double * altitude) const = 0;

};

} // namespace exif
} // namespace openMVG

#endif // OPENMVG_EXIF_EXIF_IO_HPP

