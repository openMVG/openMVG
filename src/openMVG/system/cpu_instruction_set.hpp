// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (C) 2017 Pierre Moulon

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SYSTEM_CPUID_HPP
#define OPENMVG_SYSTEM_CPUID_HPP

#include <array>
#include <bitset>

#if defined _MSC_VER
  #include <intrin.h>
#else
  #include <cpuid.h>
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
    std::array<int, 4> cpui;
    if (internal_cpuid(cpui.data(), 0))
    {
      const int nIds = cpui[0];
      if (nIds<1){
        // possible old cpu
        return;
      }

      internal_cpuid(cpui.data(), 1);

      const std::bitset<32> Edx (cpui[3]);
      m_SSE = Edx[25];
      m_SSE2 = Edx[26];

      const std::bitset<32> Ecx (cpui[2]);
      m_AVX = Ecx[28];
      m_SSE3 = Edx[0];
      m_SSE41 = Ecx[19];
      m_SSE42 = Ecx[20];
      m_POPCNT = Ecx[23];

      if (nIds > 6)
      {
        internal_cpuid(cpui.data(), 7);
        const std::bitset<32> Ebx (cpui[1]);
        m_AVX2 = Ebx[5];
      }
    }
  }

  bool supportSSE() const
  {
    return m_SSE;
  }

  bool supportSSE2() const
  {
    return m_SSE2;
  }

  bool supportSSE3() const
  {
    return m_SSE3;
  }

  bool supportSSE41() const
  {
    return m_SSE41;
  }

  bool supportSSE42() const
  {
    return m_SSE42;
  }

  bool supportAVX() const
  {
    return m_AVX;
  }

  bool supportAVX2() const
  {
    return m_AVX2;
  }

  bool supportPOPCNT() const
  {
    return m_POPCNT;
  }

private:
  static bool internal_cpuid(int32_t out[4], int32_t x)
  {
    #if defined __GNUC__
    __cpuid_count(x, 0, out[0], out[1], out[2], out[3]);
    return true;
    #endif
    #if defined _MSC_VER
    __cpuidex(out, x, 0);
    return true;
    #endif
    return false;
  }
};

} // namespace system
} // namespace openMVG

#endif // CPUID
