
// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_FEATURES_DESCRIPTOR_HPP
#define OPENMVG_FEATURES_DESCRIPTOR_HPP

#include "openMVG/numeric/numeric.h"
#include <iostream>
#include <iterator>
#include <fstream>
#include <string>
#include <vector>

namespace openMVG {
namespace features {

/**
 * Class that handle descriptor (a data container of N values of type T).
 * SiftDescriptor => <uchar,128> or <float,128>
 * Surf 64 => <float,64>
 * Brief 512 bits => <unsigned char,512/sizeof(unsigned char)>
 */
template <typename T, std::size_t N>
class Descriptor : public Eigen::Matrix<T, N, 1>
{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  typedef T bin_type;
  typedef std::size_t size_type;

  /// Compile-time length of the descriptor
  static const std::size_t static_size = N;

  /// Ostream interface
  std::ostream& print(std::ostream& os) const;
  /// Istream interface
  std::istream& read(std::istream& in);

  template<class Archive>
  void save(Archive & archive) const
  {
    const std::vector<T> array(this->data(), this->data()+N);
    archive( array );
  }

  template<class Archive>
  void load(Archive & archive)
  {
    std::vector<T> array(N);
    archive( array );
    std::memcpy(this->data(), array.data(), sizeof(T)*N);
  }
};

// Output stream definition
template <typename T, std::size_t N>
inline std::ostream& operator<<(std::ostream& out, const Descriptor<T, N>& obj)
{
  return obj.print(out); //simply call the print method.
}

// Input stream definition
template <typename T, std::size_t N>
inline std::istream& operator>>(std::istream& in, Descriptor<T, N>& obj)
{
  return obj.read(in); //simply call the read method.
}

//-- Use specialization to handle unsigned char case.
//-- We do not want confuse unsigned char value with the spaces written in the file

template<typename T>
inline std::ostream& printT(std::ostream& os, T *tab, size_t N)
{
  std::copy( tab, &tab[N], std::ostream_iterator<T>(os," "));
  return os;
}

template<>
inline std::ostream& printT<unsigned char>(std::ostream& os, unsigned char *tab, size_t N)
{
  for(size_t i=0; i < N; ++i)
    os << (int)tab[i] << " ";
  return os;
}

template<typename T>
inline std::istream& readT(std::istream& is, T *tab, size_t N)
{
  for(size_t i=0; i<N; ++i) is >> tab[i];
  return is;
}

template<>
inline std::istream& readT<unsigned char>(std::istream& is, unsigned char *tab, size_t N)
{
  int temp = -1;
  for(size_t i=0; i < N; ++i){
    is >> temp; tab[i] = (unsigned char)temp;
  }
  return is;
}

template<typename T, std::size_t N>
std::ostream& Descriptor<T,N>::print(std::ostream& os) const
{
  return printT<T>(os, (T*) this->data(), N);
}

template<typename T, std::size_t N>
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

  std::copy(
    std::istream_iterator<typename DescriptorsT::value_type >(fileIn),
    std::istream_iterator<typename DescriptorsT::value_type >(),
    std::back_inserter(vec_desc));
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
  typedef typename DescriptorsT::value_type VALUE;

  vec_desc.clear();
  std::ifstream fileIn(sfileNameDescs.c_str(), std::ios::in | std::ios::binary);
  if (!fileIn.is_open())
    return false;
  //Read the number of descriptor in the file
  std::size_t cardDesc = 0;
  fileIn.read((char*) &cardDesc,  sizeof(std::size_t));
  vec_desc.resize(cardDesc);
  for (typename DescriptorsT::const_iterator iter = vec_desc.begin();
    iter != vec_desc.end(); ++iter) {
    fileIn.read((char*) (*iter).data(),
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
  typedef typename DescriptorsT::value_type VALUE;

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
