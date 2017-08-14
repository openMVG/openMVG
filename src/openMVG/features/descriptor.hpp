// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_FEATURES_DESCRIPTOR_HPP
#define OPENMVG_FEATURES_DESCRIPTOR_HPP

#include <algorithm>
#include <array>
#include <iostream>
#include <iterator>
#include <fstream>
#include <string>

#include "openMVG/numeric/eigen_alias_definition.hpp"

namespace openMVG {
namespace features {

/**
 * Class that handle descriptor (a data container of N values of type T).
 * SiftDescriptor => <uchar,128> or <float,128>
 * Surf 64 => <float,64>
 * Brief 512 bits => <unsigned char,512/sizeof(unsigned char)>
 */
template <typename T, uint32_t N>
class Descriptor : public Eigen::Matrix<T, N, 1>
{
public:
  using bin_type = T;
  using size_type = uint32_t;

  /// Compile-time length of the descriptor
  static const uint32_t static_size = N;

  /// ostream interface
  std::ostream& print(std::ostream& os) const;
  /// istream interface
  std::istream& read(std::istream& in);

  template<class Archive>
  void save(Archive & archive) const
  {
    std::array<T,N> array;
    std::copy(this->data(), this->data()+N, array.begin());
    archive( array );
  }

  template<class Archive>
  void load(Archive & archive)
  {
    std::array<T, N> array;
    archive( array );
    std::memcpy(this->data(), array.data(), sizeof(T)*N);
  }
};

// Output stream definition
template <typename T, uint32_t N>
inline std::ostream& operator<<(std::ostream& out, const Descriptor<T, N>& obj)
{
  return obj.print(out); //simply call the print method.
}

// Input stream definition
template <typename T, uint32_t N>
inline std::istream& operator>>(std::istream& in, Descriptor<T, N>& obj)
{
  return obj.read(in); //simply call the read method.
}

//-- Use specialization to handle unsigned char case.
//-- We do not want confuse unsigned char value with the spaces written in the file

template<typename T>
inline std::ostream& printT(std::ostream& os, T *tab, uint32_t N)
{
  std::copy( tab, &tab[N], std::ostream_iterator<T>(os," "));
  return os;
}

template<>
inline std::ostream& printT<unsigned char>(std::ostream& os, unsigned char *tab, uint32_t N)
{
  for (uint32_t i=0; i < N; ++i)
    os << static_cast<int>(tab[i]) << " ";
  return os;
}

template<typename T>
inline std::istream& readT(std::istream& is, T *tab, uint32_t N)
{
  for (uint32_t i=0; i<N; ++i)
    is >> tab[i];
  return is;
}

template<>
inline std::istream& readT<unsigned char>(std::istream& is, unsigned char *tab, uint32_t N)
{
  int temp = -1;
  for (uint32_t i=0; i < N; ++i){
    is >> temp; tab[i] = static_cast<unsigned char>(temp);
  }
  return is;
}

template<typename T, uint32_t N>
std::ostream& Descriptor<T,N>::print(std::ostream& os) const
{
  return printT<T>(os, (T*) this->data(), N);
}

template<typename T, uint32_t N>
std::istream& Descriptor<T,N>::read(std::istream& in)
{
  return readT<T>(in, (T*) this->data(), N);
}

/// Read descriptors from file
template<typename DescriptorsT >
inline bool loadDescsFromFile(
  const std::string & sfileNameDescs,
  DescriptorsT & vec_desc)
{
  vec_desc.clear();

  std::ifstream fileIn(sfileNameDescs.c_str());
  if (!fileIn.is_open())
    return false;

  typename DescriptorsT::value_type value;
  while (fileIn >> value) {
    vec_desc.emplace_back(value);
  }

  const bool bOk = !fileIn.bad();
  fileIn.close();
  return bOk;
}

/// Write descriptors to file
template<typename DescriptorsT >
inline bool saveDescsToFile(
  const std::string & sfileNameDescs,
  DescriptorsT & vec_desc)
{
  std::ofstream file(sfileNameDescs.c_str());
  if (!file.is_open())
    return false;
  std::copy(vec_desc.begin(), vec_desc.end(),
            std::ostream_iterator<typename DescriptorsT::value_type >(file,"\n"));
  const bool bOk = file.good();
  file.close();
  return bOk;
}


/// Read descriptors from file (in binary mode)
template<typename DescriptorsT >
inline bool loadDescsFromBinFile(
  const std::string & sfileNameDescs,
  DescriptorsT & vec_desc)
{
  using VALUE = typename DescriptorsT::value_type;

  vec_desc.clear();
  std::ifstream fileIn(sfileNameDescs.c_str(), std::ios::in | std::ios::binary);
  if (!fileIn.is_open())
    return false;
  //Read the number of descriptor in the file
  std::size_t cardDesc = 0;
  fileIn.read(reinterpret_cast<char*>(&cardDesc), sizeof(std::size_t));
  vec_desc.resize(cardDesc);
  for (auto & it :vec_desc) {
    fileIn.read(reinterpret_cast<char*>(it.data()),
      VALUE::static_size*sizeof(typename VALUE::bin_type));
  }
  const bool bOk = !fileIn.bad();
  fileIn.close();
  return bOk;
}

/// Write descriptors to file (in binary mode)
template<typename DescriptorsT >
inline bool saveDescsToBinFile(
  const std::string & sfileNameDescs,
  DescriptorsT & vec_desc)
{
  using VALUE = typename DescriptorsT::value_type;

  std::ofstream file(sfileNameDescs.c_str(), std::ios::out | std::ios::binary);
  if (!file.is_open())
    return false;
  //Write the number of descriptor
  const std::size_t cardDesc = vec_desc.size();
  file.write((const char*) &cardDesc,  sizeof(std::size_t));
  for (typename DescriptorsT::const_iterator iter = vec_desc.begin();
    iter != vec_desc.end(); ++iter) {
    file.write((const char*) (*iter).data(),
      VALUE::static_size*sizeof(typename VALUE::bin_type));
  }
  const bool bOk = file.good();
  file.close();
  return bOk;
}

} // namespace features
} // namespace openMVG

#endif  // OPENMVG_FEATURES_DESCRIPTOR_HPP
