
// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_FEATURES_FEATURE_HPP
#define OPENMVG_FEATURES_FEATURE_HPP

#include "openMVG/numeric/numeric.h"
#include <iostream>
#include <iterator>
#include <fstream>
#include <string>
#include <vector>

namespace openMVG {

/**
 * Abstract base class for features.
 */
class FeatureBase {
public:
  virtual inline ~FeatureBase() {};

  //this pure virtual method will be called to print the derived class' object.
  virtual std::ostream& print(std::ostream& output) const = 0;

  //this pure virtual method will be called to read the derived class' object.
  virtual std::istream& read(std::istream& input) = 0;
};

//with overloaded operators:
inline std::ostream& operator<<(std::ostream& out, const FeatureBase& obj)
{
  return obj.print(out); //simply call the print method.
}

inline std::istream& operator>>(std::istream& in, FeatureBase& obj)
{
  return obj.read(in); //simply call the read method.
}

/**
 * Base class for Point features.
 * Store position of the feature point.
 */
class PointFeature : public FeatureBase {
public:
  virtual inline ~PointFeature() {};

  inline PointFeature(float x=0.0f, float y=0.0f)
   : _coords(x, y) {}

  inline float x() const { return _coords(0); }
  inline float y() const { return _coords(1); }
  inline Vec2f coords() const { return _coords;}

  inline float& x() { return _coords(0); }
  inline float& y() { return _coords(1); }
  inline Vec2f& coords() { return _coords;}

  inline std::ostream& print(std::ostream& os) const
  { return os << _coords(0) << " " << _coords(1); }

  inline std::istream& read(std::istream& in)
  { return in >> _coords(0) >> _coords(1); }

protected:
  Vec2f _coords;  // (x, y).
};

/**
 * Base class for ScaleInvariant Oriented Point features.
 * Add scale and orientation description to basis PointFeature.
 */
class SIOPointFeature : public PointFeature {
public:
  virtual ~SIOPointFeature() {};

  SIOPointFeature(float x=0.0f, float y=0.0f,
                  float scale=0.0f, float orient=0.0f)
    : PointFeature(x,y)
    , _scale(scale)
    , _orientation(orient) {}

  inline float scale() const { return _scale; }
  inline float& scale() { return _scale; }
  inline float orientation() const { return _orientation; }
  inline float& orientation() { return _orientation; }

  bool operator ==(const SIOPointFeature& b) const {
    return (_scale == b.scale()) &&
           (_orientation == b.orientation()) &&
           (x() == b.x()) && (y() == b.y()) ;
  };

  std::ostream& print(std::ostream& os) const
  {
    return PointFeature::print(os) << " " << _scale << " " << _orientation;
  }

  std::istream& read(std::istream& in)
  {
    return PointFeature::read(in) >> _scale >> _orientation;
  }

protected:
  float _scale;        // In pixels.
  float _orientation;  // In radians.
};

/// Read feats from file
template<typename FeaturesT >
static bool loadFeatsFromFile(
  const std::string & sfileNameFeats,
  FeaturesT & vec_feat)
{
  vec_feat.clear();
  bool bOk = false;

  std::ifstream fileIn(sfileNameFeats.c_str());
  std::copy(
    std::istream_iterator<typename FeaturesT::value_type >(fileIn),
    std::istream_iterator<typename FeaturesT::value_type >(),
    std::back_inserter(vec_feat));
  bOk = !fileIn.bad();
  fileIn.close();
  return bOk;
}

/// Write feats to file
template<typename FeaturesT >
static bool saveFeatsToFile(
  const std::string & sfileNameFeats,
  FeaturesT & vec_feat)
{
  std::ofstream file(sfileNameFeats.c_str());
  std::copy(vec_feat.begin(), vec_feat.end(),
            std::ostream_iterator<typename FeaturesT::value_type >(file,"\n"));
  bool bOk = file.good();
  file.close();
  return bOk;
}

} // namespace openMVG

#endif // OPENMVG_FEATURES_FEATURE_HPP
