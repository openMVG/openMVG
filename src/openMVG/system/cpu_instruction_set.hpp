// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (C) 2017 Pierre Moulon

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SYSTEM_CPUID_HPP
#define OPENMVG_SYSTEM_CPUID_HPP

#if defined _MSC_VER
  #include <array>
  #include <bitset>
  #include <intrin.h>
#endif

namespace openMVG
{
/**
* @brief Namespace handling various functions and classes about system programming
*/
namespace system
{

/**
* @brief Class to check runtime support of various CPU intrinsic set CPU capability
*/
class CpuInstructionSet
{
  bool m_SSE = false;
  bool m_SSE2 = false;
  bool m_SSE3 = false;
  bool m_SSE41 = false;
  bool m_SSE42 = false;
  bool m_AVX = false;
  bool m_AVX2 = false;
  bool m_POPCNT = false;

  public:

  CpuInstructionSet()
  {
#if defined _MSC_VER
    std::array<int, 4> cpui;
    __cpuid(cpui.data(),0);

    const int nIds = cpui[0];

    if (nIds<1){
     // possible old cpu
     return;
    }

    std::array<int, 4> cpui_ext;
    __cpuidex(cpui_ext.data(),1,0);

    const std::bitset<32> Edx (cpui_ext[3]);
    m_SSE = Edx[25];
    m_SSE2 = Edx[26];

    const std::bitset<32> Ecx (cpui_ext[2]);
    m_AVX = Ecx[28];
    m_SSE3 = Edx[0];
    m_SSE41 = Ecx[19];
    m_SSE42 = Ecx[20];
    m_POPCNT = Ecx[23];

    if (nIds > 6)
    {
      __cpuidex(cpui_ext.data(),7,0);
      const std::bitset<32> Ebx (cpui_ext[1]);
      m_AVX2 = Ebx[5];
    }
#endif
  }

  bool supportSSE() const
  {
  #if defined __GNUC__
    return __builtin_cpu_supports("sse");
  #else
    return m_SSE;
  #endif
  }

  bool supportSSE2() const
  {
  #if defined __GNUC__
    return __builtin_cpu_supports("sse2");
  #else
    return m_SSE2;
  #endif
  }

  bool supportSSE3() const
  {
  #if defined __GNUC__
    return __builtin_cpu_supports("sse3");
  #else
    return m_SSE3;
  #endif
  }

  bool supportSSE41() const
  {
  #if defined __GNUC__
    return __builtin_cpu_supports("sse4.1");
  #else
    return m_SSE41;
  #endif
  }

  bool supportSSE42() const
  {
  #if defined __GNUC__
    return __builtin_cpu_supports("sse4.2");
  #else
    return m_SSE42;
  #endif
  }

  bool supportAVX() const
  {
  #if defined __GNUC__
    return __builtin_cpu_supports("avx");
  #else
    return m_AVX;
  #endif
  }

  bool supportAVX2() const
  {
  #if defined __GNUC__
    return __builtin_cpu_supports("avx2");
  #else
    return m_AVX2;
  #endif
  }

  bool supportPOPCNT() const
  {
  #if defined __GNUC__
    return __builtin_cpu_supports("popcnt");
  #else
    return m_POPCNT;
  #endif
  }
};

} // namespace system
} // namespace openMVG

#endif // CPUID
