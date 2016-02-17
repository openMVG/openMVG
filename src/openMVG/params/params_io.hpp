
// Copyright (c) 2015 Pierre MOULON.
// Copyright (c) 2015 Klemen Istenic (CIRS - Universitat de Girona).

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


#ifndef OPENMVG_PARAMS_IO_HPP
#define OPENMVG_PARAMS_IO_HPP

#include "openMVG/params/params.hpp"

namespace openMVG {
namespace params {

// Camera data

/// Load camera parameters from file to paramsCamera
bool Load(paramsCamera & params_data, const std::string & filename);

/// Save camera parameters to file
bool Save(const paramsCamera & params_data, const std::string & filename);

// Detection data

/// Load detection parameters from file to paramsDetection
bool Load(paramsDetection & params_data, const std::string & filename);

/// Save detection parameters to file
bool Save(const paramsDetection & params_data, const std::string & filename);

// Matching data

/// Load matching parameters from file to paramsMatching
bool Load(paramsMatching & params_data, const std::string & filename);

/// Save matching parameters to file
bool Save(const paramsMatching & params_data, const std::string & filename);

// Incremental SfM data
/// Load incremental SfM parameters from file to paramsIncrementalSfM
bool Load(paramsIncrementalSfM & params_data, const std::string & filename);

/// Save incremental SfM parameters to file
bool Save(const paramsIncrementalSfM & params_data, const std::string & filename);

// Global SfM data

/// Load global SfM parameters from file to paramsGlobalSfM
bool Load(paramsGlobalSfM & params_data, const std::string & filename);

/// Save global SfM parameters to file
bool Save(const paramsGlobalSfM & params_data, const std::string & filename);


} // namespace params
} // namespace openMVG

#endif // OPENMVG_PARAMS_IO_HPP
