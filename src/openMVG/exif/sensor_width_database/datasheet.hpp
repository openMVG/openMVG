// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2013 Pierre Moulon, Bruno Duisit.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_EXIF_SENSOR_WIDTH_DATABASE_DATASHEET_HPP
#define OPENMVG_EXIF_SENSOR_WIDTH_DATABASE_DATASHEET_HPP

#include <algorithm>
#include <string>
#include <vector>

#include "openMVG/stl/split.hpp"

// Database structure to store camera model and sensor size
struct Datasheet
{
  Datasheet() = default;

  Datasheet
  (
    const std::string& model,
    const double& sensorSize
  ):
    model_(model),
    sensorSize_(sensorSize)
  {}

  bool operator==(const Datasheet& rhs) const
  {
    // Check if brand & model are equals or not
    // Camera model is construct as the following:
    // MAKER MODELNAME
    // There is some point to care about for the comparison:
    // since camera model can be formatted differently:
    //  -> "KODAK Z612 ZOOM DIGITAL CAMERA" to "Kodak EasyShare Z612"
    // - we use lower case comparison
    // - we ensure that DIGIT based substring have a match

    std::string this_model = model_;
    std::transform(this_model.begin(), this_model.end(), this_model.begin(), ::tolower);
    std::string rhs_model = rhs.model_;
    std::transform(rhs_model.begin(), rhs_model.end(), rhs_model.begin(), ::tolower);

    // retrieve camera maker sustring [0->first space]
    std::string this_maker( this_model );
    this_maker = this_maker.substr( 0, this_maker.find( ' ' ) );
    std::string rhs_maker( rhs_model );
    rhs_maker = rhs_maker.substr( 0, rhs_maker.find( ' ' ) );

    if (this_maker.compare(rhs_maker) == 0)
    {
      // Extract digit substring from database entry
      std::vector<std::string> vec_db_model;
      stl::split(this_model, ' ', vec_db_model);
      std::string db_model_digit_substring = "";
      for (const std::string & db_sub_model : vec_db_model)
      {
        if ( std::find_if(db_sub_model.begin(), db_sub_model.end(), isdigit) != db_sub_model.end() )
        {
          db_model_digit_substring = db_sub_model;
          break;
        }
      }

      // Search for digit substring in image camera model and if there is a match

      const bool has_digit =
        std::find_if(rhs_model.begin(), rhs_model.end(), isdigit) != rhs_model.end();
      if (!has_digit || db_model_digit_substring.empty() )
      {
        // If there is no digit in either the db entry or the image camera then compare the full model string
        return (this_model.compare(rhs_model) == 0);
      }
      else
      {
        // Split the model string and ensure that the digit part is found
        std::vector<std::string> vec_model_rhs;
        stl::split(rhs_model, ' ', vec_model_rhs);
        for (const std::string & rhs_sub_model : vec_model_rhs)
        {
          // Check if we have a digit based part
          const bool sub_has_digit =
            std::find_if(rhs_sub_model.begin(), rhs_sub_model.end(), isdigit) != rhs_sub_model.end();
          // Check that substring containing the camera digit model is corresponding
          if (sub_has_digit && ( db_model_digit_substring == rhs_sub_model ) )
          {
            return true; // substring that contains digit got a match
          }
        }
      }
    }
    return false;
  }

  std::string model_;
  double sensorSize_;
};


#endif // OPENMVG_EXIF_SENSOR_WIDTH_DATABASE_DATASHEET_HPP
