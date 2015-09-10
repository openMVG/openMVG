
// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "openMVG/sfm/sfm_data.hpp"
#include <string>

namespace openMVG {
namespace sfm {

/// Basic Reconstruction Engine.
/// Process Function handle the reconstruction.
class ReconstructionEngine
{
public:

  ReconstructionEngine(
    const SfM_Data & sfm_data,
    const std::string & soutDirectory)
    :_sOutDirectory(soutDirectory),
    _sfm_data(sfm_data),
    _bFixedIntrinsics(false)
  {
  }

  virtual ~ReconstructionEngine() {}

  virtual bool Process() = 0;

  bool Get_bFixedIntrinsics() const {return _bFixedIntrinsics;}
  void Set_bFixedIntrinsics(bool bVal) {_bFixedIntrinsics = bVal;}

  const SfM_Data & Get_SfM_Data() const {return _sfm_data;}

protected:
  std::string _sOutDirectory; // Output path where outputs will be stored

  //-----
  //-- Reconstruction data
  //-----
  SfM_Data _sfm_data; // internal SfM_Data

  //-----
  //-- Reconstruction parameters
  //-----
  bool _bFixedIntrinsics;
};

} // namespace sfm
} // namespace openMVG
