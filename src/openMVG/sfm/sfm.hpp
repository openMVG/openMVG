// Copyright (c) 2015 Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_HPP
#define OPENMVG_SFM_HPP

#include "openMVG/types.hpp"
#include "openMVG/numeric/numeric.h"

// Serialization
#include <cereal/cereal.hpp>

//-----------------
// SfM data
//-----------------
#include "openMVG/sfm/sfm_view.hpp"
#include "openMVG/sfm/sfm_landmark.hpp"

#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"
#include "openMVG/sfm/sfm_data_filters.hpp"
#include "openMVG/sfm/sfm_data_filters_frustum.hpp"
#include "openMVG/sfm/sfm_data_BA.hpp"
#include "openMVG/sfm/sfm_data_BA_ceres.hpp"

#include "openMVG/sfm/sfm_filters.hpp"
#include "openMVG/sfm/sfm_data_triangulation.hpp"

//-----------------
// SfM pipelines
//-----------------
#include "openMVG/sfm/pipelines/sfm_engine.hpp"
#include "openMVG/sfm/pipelines/sfm_features_provider.hpp"
#include "openMVG/sfm/pipelines/sfm_matches_provider.hpp"

#include "openMVG/sfm/pipelines/sfm_robust_model_estimation.hpp"

#include "openMVG/sfm/pipelines/global/mutexSet.hpp"
#include "openMVG/sfm/pipelines/global/sfm_global_reindex.hpp"
#include "openMVG/sfm/pipelines/global/sfm_global_engine_relative_motions.hpp"

#endif // OPENMVG_SFM_HPP
