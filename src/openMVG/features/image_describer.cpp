
// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <openMVG/features/image_describer.hpp>

#include <exception>

namespace openMVG{
namespace features{

EDESCRIBER_PRESET describerPreset_stringToEnum(const std::string& sPreset)
{
  if(sPreset == "LOW")
    return LOW_PRESET;
  if (sPreset == "MEDIUM")
    return MEDIUM_PRESET;
  if(sPreset == "NORMAL")
    return NORMAL_PRESET;
  if (sPreset == "HIGH")
    return HIGH_PRESET;
  if (sPreset == "ULTRA")
    return ULTRA_PRESET;
  throw std::invalid_argument("Invalid descriptor preset: " + sPreset);
}

std::string describerPreset_enumToString(const EDESCRIBER_PRESET preset)
{
  if(preset == LOW_PRESET)
    return "LOW";
  if (preset == MEDIUM_PRESET)
    return "MEDIUM";
  if(preset == NORMAL_PRESET)
    return "NORMAL";
  if (preset == HIGH_PRESET)
    return "HIGH";
  if (preset == ULTRA_PRESET)
    return "ULTRA";
  throw std::invalid_argument("Unrecognized EDESCRIBER_PRESET "+std::to_string(preset));
}

std::ostream& operator<<(std::ostream& os, EDESCRIBER_PRESET p)
{
    return os << describerPreset_enumToString(p);
}

std::istream& operator>>(std::istream& in, EDESCRIBER_PRESET& p)
{
    std::string token;
    in >> token;
    p = describerPreset_stringToEnum(token);
    return in;
}

}//namespace features
}//namespace openMVG