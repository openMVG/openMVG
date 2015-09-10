// Copyright (c) 2015 Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

namespace openMVG {
namespace sfm {

struct SfM_Data;

// Group camera models that share common camera properties
// It modifies the intrinsic_id of the view field and change the sfm_data.intrinsics length
// Grouping is simplified by using a hash function over the camera intrinsics
// - it allow to merge camera model that share common camera parameters & image sizes
void GroupSharedIntrinsics(SfM_Data & sfm_data);

} // namespace sfm
} // namespace openMVG
