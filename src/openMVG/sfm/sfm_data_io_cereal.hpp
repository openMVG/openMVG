
// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


#ifndef OPENMVG_SFM_SFM_DATA_IO_CEREAL_HPP
#define OPENMVG_SFM_SFM_DATA_IO_CEREAL_HPP

#include "openMVG/sfm/sfm_data_io.hpp"

namespace openMVG {
namespace sfm {

/// Load a SfM_Data SfM scene from a file using the Cereal Archive interface
template <
// JSONInputArchive/ ...
typename archiveType
>
bool Load_Cereal(
  SfM_Data & data,
  const std::string & filename,
  ESfM_Data flags_part);

/// Save a SfM_Data SfM scene to a file using the Cereal Archive interface
template <
// JSONOutputArchive/ ...
typename archiveType
>
bool Save_Cereal(
  const SfM_Data & data,
  const std::string & filename,
  ESfM_Data flags_part);

} // namespace sfm
} // namespace openMVG

//
// Explicit template instantiation
//

namespace cereal {
  class BinaryInputArchive;
  class BinaryOutputArchive;

 class PortableBinaryInputArchive;
 class PortableBinaryOutputArchive;

 class JSONInputArchive;
 class JSONOutputArchive;

 class XMLInputArchive;
 class XMLOutputArchive;
} // namespace cereal

extern template bool openMVG::sfm::Load_Cereal
<cereal::BinaryInputArchive>
(
  openMVG::sfm::SfM_Data & data,
  const std::string & filename,
  openMVG::sfm::ESfM_Data flags_part
);

extern template bool openMVG::sfm::Load_Cereal
<cereal::PortableBinaryInputArchive>
(
  openMVG::sfm::SfM_Data & data,
  const std::string & filename,
  openMVG::sfm::ESfM_Data flags_part
);

extern template bool openMVG::sfm::Load_Cereal
<cereal::JSONInputArchive>
(
  openMVG::sfm::SfM_Data & data,
  const std::string & filename,
  openMVG::sfm::ESfM_Data flags_part
);

extern template bool openMVG::sfm::Load_Cereal
<cereal::XMLInputArchive>
(
  openMVG::sfm::SfM_Data & data,
  const std::string & filename,
  openMVG::sfm::ESfM_Data flags_part
);

extern template bool openMVG::sfm::Save_Cereal
<cereal::BinaryOutputArchive>
(
  const openMVG::sfm::SfM_Data & data,
  const std::string & filename,
  openMVG::sfm::ESfM_Data flags_part
);

extern template bool openMVG::sfm::Save_Cereal
<cereal::PortableBinaryOutputArchive>
(
  const openMVG::sfm::SfM_Data & data,
  const std::string & filename,
  openMVG::sfm::ESfM_Data flags_part
);

extern template bool openMVG::sfm::Save_Cereal
<cereal::JSONOutputArchive>
(
  const openMVG::sfm::SfM_Data & data,
  const std::string & filename,
  openMVG::sfm::ESfM_Data flags_part
);

extern template bool openMVG::sfm::Save_Cereal
<cereal::XMLOutputArchive>
(
  const openMVG::sfm::SfM_Data & data,
  const std::string & filename,
  openMVG::sfm::ESfM_Data flags_part
);

#endif // OPENMVG_SFM_SFM_DATA_IO_CEREAL_HPP
