
// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_ENGINE_H
#define OPENMVG_SFM_ENGINE_H

#include <string>

namespace openMVG{

/// Basic Reconstruction Engine.
/// Provide essentially path to the data.
/// Process Function handle the reconstruction.
class OldReconstructionEngine
{
public:

  OldReconstructionEngine(const std::string & sSfM_Data_Path,
    const std::string & smatchesPath,
    const std::string & soutDirectory)
    :_sSfM_Data_Path(sSfM_Data_Path),
    _sMatchesPath(smatchesPath),
    _sOutDirectory(soutDirectory)
  {
  }

  virtual bool Process() = 0;

protected:
  std::string _sSfM_Data_Path;    // Path to the Sfm_Scene
  std::string _sMatchesPath;  // Path to correspondences and features
  std::string _sOutDirectory; // Output path where outputs will be stored
};
} // namespace openMVG

#endif // OPENMVG_SFM_ENGINE_H
