// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2013-2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_EXIF_EXIF_IO_EASYEXIF_HPP
#define OPENMVG_EXIF_EXIF_IO_EASYEXIF_HPP

#include <memory>
#include <string>

#include "openMVG/exif/exif_IO.hpp"

namespace openMVG
{
namespace exif
{

/**
* @brief Class managing EXIF data using EasyEXIF library
*/
class Exif_IO_EasyExif : public Exif_IO
{
  public:

    /**
    * @brief Default constructor
    */
    Exif_IO_EasyExif();

    /**
    * @brief Constructor using a file
    * @param sFileName path of the image to analyze
    */
    explicit Exif_IO_EasyExif( const std::string & sFileName );

    /**
    * @brief Open and populate EXIF data
    * @param sFileName path of the image to analyze
    * @retval true if image file could be parsed correctly
    * @retval false if image file could not be parsed and analysed
    */
    bool open( const std::string & sFileName ) override;

    /**
    * @brief Get image width
    * @return Width of the image (in pixel)
    */
    size_t getWidth() const override;

    /**
    * @brief Get image height
    * @return Height of the image (in pixel)
    */
    size_t getHeight() const override;

    /**
    * @brief Get Focal (in mm)
    * @return Focal of the lens when image was shot (in mm)
    */
    float getFocal() const override;

    /**
    * @brief Get FocalLengthIn35mm (in mm)
    * @return The equivalent focal length assuming a 35mm film camera, in mm.
    */
    float getFocalLengthIn35mm() const override;

    /**
    * @brief Get FocalPlaneXResolution
    * @return Number of pixels in the image width (X) direction per
    *           FocalPlaneResolutionUnit on the camera focal plane.
    */
    float getFocalPlaneXResolution() const override;

    /**
    * @brief Get FocalPlaneYResolution
    * @return Number of pixels in the image height (Y) direction per
    *           FocalPlaneResolutionUnit on the camera focal plane.
    */
    float getFocalPlaneYResolution() const override;

    /**
    * @brief Get FocalPlaneResolutionUnit
    *        Unit -> 2: inch, 3: centimeter, 4: millimeter, 5: micrometer.
    * @return Indicates the unit for measuring FocalPlaneXResolution and
    *          FocalPlaneYResolution.
    */
    int getFocalPlaneResolutionUnit() const override;

    /**
    * @brief Get Brand of the camera
    * @return Brand name
    */
    std::string getBrand() const override;

    /**
    * @brief Get Model of the camera
    * @return Camera model name
    */
    std::string getModel() const override;

    /**
    * @brief Get Lens model
    * @return Lens model name
    */
    std::string getLensModel() const override;

    /**
    * @brief Get a unique indentifier for this image
    * @return Unique ID
    */
    std::string getImageUniqueID() const override;

    /**
      * @brief Verify if the file has metadata
      * @retval true If EXIF value are present
      * @retval false If EXIF value are not present
      */
    bool doesHaveExifInfo() const override;

    /**
    * @brief Print all data
    * @return string containing all EXIF data
    */
    std::string allExifData() const override;

    /**
    * @brief Try to read and save the EXIF GPS latitude
    * @return If GPS Latitude can be read & exported, return true
    */
    bool GPSLatitude(double * latitude) const override;

    /**
    * @brief Try to read and save the EXIF GPS longitude
    * @return If GPS Longitude can be read & exported, return true
    */
    bool GPSLongitude(double * longitude) const override;

   /**
    * @brief Try to read and save the EXIF GPS altitude
    * @return If GPS Altitude can be read & exported, return true
    */
    bool GPSAltitude(double * altitude) const  override;


  private:

    /// Hide the easyexif::EXIFInfo to a Pimp "Pointer to implementation".
    class EXIFInfoImpl;
    std::unique_ptr<EXIFInfoImpl> pimpl_;

    /// Indicate if exifInfo_ is populated
    bool bHaveExifInfo_;
};

} // namespace exif
} // namespace openMVG

#endif // OPENMVG_EXIF_EXIF_IO_EASYEXIF_HPP
