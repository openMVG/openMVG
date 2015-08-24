// Copyright (c) 2013-2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "exif_IO.hpp"

#include "openMVG/stl/hash.hpp"


namespace openMVG {
namespace exif  {

std::size_t computeUID(const Exif_IO& exifReader, const std::string& imageFilename)
{
  std::size_t uid = 0;

  if( !exifReader.getImageUniqueID().empty() ||
      !exifReader.getSerialNumber().empty() ||
      !exifReader.getLensSerialNumber().empty()
    )
  {
    stl::hash_combine(uid, exifReader.getImageUniqueID());
    stl::hash_combine(uid, exifReader.getSerialNumber());
    stl::hash_combine(uid, exifReader.getLensSerialNumber());
  }
  else
  {
    // No metadata to identify the image, fallback to the filename
    stl::hash_combine(uid, imageFilename);
  }

  if( !exifReader.getSubSecTimeOriginal().empty() )
  {
    stl::hash_combine(uid, exifReader.getSubSecTimeOriginal());
  }
  else
  {
    // If no original date/time, fallback to the file date/time
    stl::hash_combine(uid, exifReader.getDateTime());
  }

  stl::hash_combine(uid, exifReader.getWidth());
  stl::hash_combine(uid, exifReader.getHeight());

  return uid;
}

} // namespace exif
} // namespace openMVG
