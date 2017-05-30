// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <cstdint>
#include <vector>

#include <cuda.h>
#include <cuda_runtime.h>

#include "openMVG/matching/indMatch.hpp"

namespace openMVG {
namespace matching_image_collection {
namespace gpu {

class LatchBitMatcher {
  public:
    LatchBitMatcher(const unsigned int);
    void match(void*, void*, int, int);
    openMVG::matching::IndMatches retrieveMatches(float);
    ~LatchBitMatcher();
  private:
    const unsigned int m_maxKP;

    void* m_dQuery;
    void* m_dTraining;

    unsigned int m_numKPQuery;
    unsigned int m_numKPTraining;

    uint32_t* m_dMatches1;
    uint32_t* m_dMatches2;

    cudaTextureObject_t m_texQuery;
    cudaTextureObject_t m_texTraining;
    // CUDA stuff
    cudaStream_t m_stream1;
    cudaStream_t m_stream2;
};

} // namespace gpu
} // namespace matching_image_collection
} // namespace openMVG
