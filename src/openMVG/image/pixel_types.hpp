// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_IMAGE_PIXELTYPES_HPP
#define OPENMVG_IMAGE_PIXELTYPES_HPP

#include "openMVG/numeric/numeric.h"

namespace openMVG {
namespace image {

/// RGB template pixel type
template <typename T>
class Rgb : public Eigen::Matrix<T, 3, 1, 0, 3, 1>
{
  typedef Eigen::Matrix<T, 3, 1, 0, 3, 1> Base;
  typedef T TBase;
public:

  //------------------------------
  //-- construction method
  inline Rgb(T red,T green,T blue) : Base(red, green, blue){}
  explicit inline Rgb(const Base& val) : Base(val) {}
  explicit inline Rgb(const T val=0) : Base(val,val,val) {}
  //-- construction method
  //------------------------------

  //------------------------------
  //-- accessors/getters methods
  inline const T& r() const { return (*this)(0); }
  inline T& r() { return (*this)(0); }
  inline const T& g() const { return (*this)(1); }
  inline T& g() { return (*this)(1); }
  inline const T& b() const { return (*this)(2); }
  inline T& b() { return (*this)(2); }
  /// Return gray
  inline operator T() const { return T(0.3*r()+0.59*g()+0.11*b());}
  //-- accessors/getters methods
  //------------------------------

  friend std::ostream& operator<<(std::ostream& os, const Rgb& col)
  {
    os << " {" ;
    for(int i=0; i<2; ++i)
      os << col(i) << ",";
    os << col(2) << "} ";
    return os;
  }

  template <typename Z>
  inline Rgb operator /(const Z& val) const
  {
    return Rgb(T((Z)((*this)(0)) / val),
               T((Z)((*this)(1)) / val),
               T((Z)((*this)(2)) / val) );
  }

  template <typename Z>
  inline Rgb operator *(const Z& val) const
  {
  	return Rgb(T((Z)(*this)(0) * val),
							 T((Z)(*this)(1) * val),
							 T((Z)(*this)(2) * val) );
  }
};
typedef Rgb<unsigned char> RGBColor;
typedef Rgb<float> RGBfColor;

/// RGBA templated pixel type
template <typename T>
class Rgba : public Eigen::Matrix<T, 4, 1, 0, 4, 1>
{
  typedef Eigen::Matrix<T, 4, 1, 0, 4, 1> Base;
public:

  //------------------------------
  //-- construction method
  inline Rgba(T red,T green,T blue, T alpha=1.0)
    : Base(red, green, blue, alpha){}
  explicit inline Rgba(const Base& val) : Base(val) {}
  explicit inline Rgba(const T val=0): Base(val,val,val,1.0) {}
  inline Rgba(const RGBColor val)
    : Base(val.r(),val.g(),val.b(),static_cast<T>(1)) {}
  //-- construction method
  //------------------------------

  //------------------------------
  //-- accessors/getters methods
  inline const T& r() const { return (*this)(0); }
  inline T& r() { return (*this)(0); }
  inline const T& g() const { return (*this)(1); }
  inline T& g() { return (*this)(1); }
  inline const T& b() const { return (*this)(2); }
  inline T& b() { return (*this)(2); }
  inline const T& a() const { return (*this)(3); }
  inline T& a() { return (*this)(3); }
  /// Return gray (not weighted by alpha)
  inline operator T() const { return T(0.3*r()+0.59*g()+0.11*b());}
  //-- accessors/getters methods
  //------------------------------

  friend std::ostream& operator<<(std::ostream& os, const Rgba& col) {
    os << " {" ;
    for(int i=0; i<3; ++i)
      os << col(i) << ",";
    os << col(3) << "} ";
    return os;
  }

  template<class Z>
  inline Rgba operator /(const Z& val) const
  {
		return Rgba(T((Z)(*this)(0) / val),
										 T((Z)(*this)(1) / val),
										 T((Z)(*this)(2) / val),
										 T((Z)(*this)(3) / val));
  }

  template<class Z>
  inline Rgba operator *(const Z& val) const
  {
  	return Rgba(T((Z)(*this)(0) * val),
										 T((Z)(*this)(1) * val),
										 T((Z)(*this)(2) * val),
										 T((Z)(*this)(3) * val));
  }
};
typedef Rgba<unsigned char> RGBAColor;

const RGBColor WHITE(255,255,255);
const RGBColor BLACK(0,0,0);
const RGBColor BLUE(0,0,255);
const RGBColor RED(255,0,0);
const RGBColor GREEN(0,255,0);
const RGBColor YELLOW(255,255,0);
const RGBColor CYAN(0,255,255);
const RGBColor MAGENTA(255,0,255);

} // namespace image
} // namespace openMVG

#endif // OPENMVG_IMAGE_PIXELTYPES_HPP

