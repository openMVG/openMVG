
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
namespace features {

/**
 * Base class for Point features.
 * Store position of a feature point.
 */
class PointFeature {

  friend std::ostream& operator<<(std::ostream& out, const PointFeature& obj);
  friend std::istream& operator>>(std::istream& in, PointFeature& obj);

public:
  PointFeature(float x=0.0f, float y=0.0f)
   : coords_(x, y) {}

  float x() const { return coords_(0); }
  float y() const { return coords_(1); }
  const Vec2f & coords() const { return coords_;}

  float& x() { return coords_(0); }
  float& y() { return coords_(1); }
  Vec2f& coords() { return coords_;}

  template<class Archive>
  void serialize(Archive & ar)
  {
    ar (coords_(0), coords_(1));
  }

protected:
  Vec2f coords_;  // (x, y).
};

typedef std::vector<PointFeature> PointFeatures;

//with overloaded operators:
inline std::ostream& operator<<(std::ostream& out, const PointFeature& obj)
{
  return out << obj.coords_(0) << " " << obj.coords_(1);
}

inline std::istream& operator>>(std::istream& in, PointFeature& obj)
{
  return in >> obj.coords_(0) >> obj.coords_(1);
}

/**
 * Base class for ScaleInvariant Oriented Point features.
 * Add scale and orientation description to basis PointFeature.
 */
class SIOPointFeature : public PointFeature {

  friend std::ostream& operator<<(std::ostream& out, const SIOPointFeature& obj);
  friend std::istream& operator>>(std::istream& in, SIOPointFeature& obj);

public:
  SIOPointFeature(float x=0.0f, float y=0.0f,
                  float scale=0.0f, float orient=0.0f)
    : PointFeature(x,y)
    , scale_(scale)
    , orientation_(orient) {}

  float scale() const { return scale_; }
  float& scale() { return scale_; }
  float orientation() const { return orientation_; }
  float& orientation() { return orientation_; }

  bool operator ==(const SIOPointFeature& b) const {
    return (scale_ == b.scale()) &&
           (orientation_ == b.orientation()) &&
           (x() == b.x()) && (y() == b.y()) ;
  };

  bool operator !=(const SIOPointFeature& b) const {
    return !((*this)==b);
  };

  template<class Archive>
  void serialize(Archive & ar)
  {
    ar (
      coords_(0), coords_(1),
      scale_,
      orientation_);
  }

protected:
  float scale_;        // In pixels.
  float orientation_;  // In radians.
};

//
inline std::ostream& operator<<(std::ostream& out, const SIOPointFeature& obj)
{
  const PointFeature *pf = static_cast<const PointFeature*>(&obj);
  return out << *pf << " " << obj.scale_ << " " << obj.orientation_;
}

inline std::istream& operator>>(std::istream& in, SIOPointFeature& obj)
{
  PointFeature *pf = static_cast<PointFeature*>(&obj);
  return in >> *pf >> obj.scale_ >> obj.orientation_;
}

/// Read feats from file
template<typename FeaturesT >
static bool loadFeatsFromFile(
  const std::string & sfileNameFeats,
  FeaturesT & vec_feat)
{
  vec_feat.clear();

  std::ifstream fileIn(sfileNameFeats.c_str());
  if (!fileIn.is_open())
    return false;

  std::copy(
    std::istream_iterator<typename FeaturesT::value_type >(fileIn),
    std::istream_iterator<typename FeaturesT::value_type >(),
    std::back_inserter(vec_feat));
  const bool bOk = !fileIn.bad();
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
  if (!file.is_open())
    return false;
  std::copy(vec_feat.begin(), vec_feat.end(),
            std::ostream_iterator<typename FeaturesT::value_type >(file,"\n"));
  const bool bOk = file.good();
  file.close();
  return bOk;
}

/// Export point feature based vector to a matrix [(x,y)'T, (x,y)'T]
template< typename FeaturesT, typename MatT >
void PointsToMat(
  const FeaturesT & vec_feats,
  MatT & m)
{
  m.resize(2, vec_feats.size());
  typedef typename FeaturesT::value_type ValueT; // Container type

  size_t i = 0;
  for( typename FeaturesT::const_iterator iter = vec_feats.begin();
    iter != vec_feats.end(); ++iter, ++i)
  {
    const ValueT & feat = *iter;
    m.col(i) << feat.x(), feat.y();
  }
}

} // namespace features
} // namespace openMVG

#endif // OPENMVG_FEATURES_FEATURE_HPP
