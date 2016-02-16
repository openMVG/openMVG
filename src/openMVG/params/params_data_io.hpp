
// Copyright (c) 2015 Pierre MOULON.
// Copyright (c) 2015 Klemen Istenic (CIRS - Universitat de Girona).

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


#ifndef OPENMVG_PARAMS_DATA_IO_HPP
#define OPENMVG_PARAMS_DATA_IO_HPP

#include "openMVG/params/params_data.hpp"

namespace openMVG {
namespace params {

/// Load parameters data from file to paramsData
bool Load(paramsData & params_data, const std::string & filename);

/// Save parameters data to file
bool Save(const paramsData & params_data, const std::string & filename);

} // namespace params
} // namespace openMVG

#endif // OPENMVG_PARAMS_DATA_IO_HPP
