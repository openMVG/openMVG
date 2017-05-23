// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2014, 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_PIPELINES_MUTEX_SET_HPP
#define OPENMVG_SFM_PIPELINES_MUTEX_SET_HPP

#include <mutex>
#include <set>

using mutexT = std::mutex;
using lock_guardT = std::lock_guard<mutexT>;

namespace openMVG {
namespace sfm{

/// ThreadSafe Set thanks to a mutex
template <typename T>
class MutexSet {

public:
    void insert(const T & value) {
      lock_guardT guard(m_Mutex);
      m_Set.insert(value);
    }

    int count(const T & value) const {
      lock_guardT guard(m_Mutex);
      return m_Set.count(value);
    }

    size_t size() const {
      lock_guardT guard(m_Mutex);
      return m_Set.size();
    }

private:
    std::set<T> m_Set;
    mutable mutexT m_Mutex;
};

} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_PIPELINES_MUTEX_SET_HPP
