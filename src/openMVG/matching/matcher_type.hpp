// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_MATCHING_MATCHER_TYPE_HPP
#define OPENMVG_MATCHING_MATCHER_TYPE_HPP

#include <stdexcept>

namespace openMVG{
namespace matching{

enum EMatcherType : unsigned char
{
  ANN_L2,
  BRUTE_FORCE_L2,
  BRUTE_FORCE_HAMMING,
  CASCADE_HASHING_L2
};

/**
 * @brief convert a string type to its corresponding enum EMatcherType
 * @param type a string corresponding to a matching method
 * @return the enum EMatcherType corresponding to the string
 */
inline EMatcherType StringToEnum(const std::string & type)
{
  if(type == "ANN_L2")                   return EMatcherType::ANN_L2;
  if(type == "BRUTE_FORCE_L2")           return EMatcherType::BRUTE_FORCE_L2;
  if(type == "BRUTE_FORCE_HAMMING")      return EMatcherType::BRUTE_FORCE_HAMMING;
  if(type == "CASCADE_HASHING_L2")       return EMatcherType::CASCADE_HASHING_L2;
  throw std::out_of_range("Invalid EMatcherType : " + type);
}

/**
 * @brief convert an enum EMatcherType to its corresponding string
 * @param matcherType the matcher type
 * @return the corresponding enum string equivalent
 */
inline std::string EnumToString(const EMatcherType & type)
{
  switch(type)
  {
    case EMatcherType::ANN_L2:                  return "ANN_L2";
    case EMatcherType::BRUTE_FORCE_L2:          return "BRUTE_FORCE_L2";
    case EMatcherType::BRUTE_FORCE_HAMMING:     return "BRUTE_FORCE_HAMMING";
    case EMatcherType::CASCADE_HASHING_L2:      return "CASCADE_HASHING_L2";
  }
  throw std::out_of_range("Invalid EMatcherType");
}

} // namespace matching
} // namespace openMVG

#endif // OPENMVG_MATCHING_MATCHER_TYPE_HPP
