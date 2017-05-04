// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_SFM_REGIONS_PROVIDER_HPP
#define OPENMVG_SFM_SFM_REGIONS_PROVIDER_HPP

#include "openMVG/features/regions.hpp"
#include "openMVG/types.hpp"
#include "third_party/progress/progress.hpp"

#include <memory>
#include <string>

namespace openMVG {
namespace sfm {
struct SfM_Data;

/// Abstract Regions provider
/// Allow to load and return the regions related to a view
struct Regions_Provider
{
public:

  virtual ~Regions_Provider() = default;

  std::string Type_id();

  bool IsScalar();

  bool IsBinary();

  const openMVG::features::Regions* getRegionsType() const;

  virtual std::shared_ptr<features::Regions> get(const IndexT x) const;

  // Load Regions related to a provided SfM_Data View container
  virtual bool load(
    const SfM_Data & sfm_data,
    const std::string & feat_directory,
    std::unique_ptr<features::Regions>& region_type,
    C_Progress *  my_progress_bar = nullptr);
  
protected:
  /// Regions per ViewId of the considered SfM_Data container
  mutable Hash_Map<IndexT, std::shared_ptr<features::Regions> > cache_;
  std::unique_ptr<openMVG::features::Regions> region_type_;
}; // Regions_Provider

} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_SFM_REGIONS_PROVIDER_HPP
