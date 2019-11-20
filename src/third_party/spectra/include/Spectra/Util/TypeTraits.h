// Copyright (C) 2018-2019 Yixuan Qiu <yixuan.qiu@cos.name>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef TYPE_TRAITS_H
#define TYPE_TRAITS_H

#include <Eigen/Core>
#include <limits>

/// \cond

namespace Spectra {


// For a real value type "Scalar", we want to know its smallest
// positive value, i.e., std::numeric_limits<Scalar>::min().
// However, we must take non-standard value types into account,
// so we rely on Eigen::NumTraits.
//
// Eigen::NumTraits has defined epsilon() and lowest(), but
// lowest() means negative highest(), which is a very small
// negative value.
//
// Therefore, we manually define this limit, and use eplison()^3
// to mimic it for non-standard types.

// Generic definition
template <typename Scalar>
struct TypeTraits
{
    static inline Scalar min()
    {
        return Eigen::numext::pow(Eigen::NumTraits<Scalar>::epsilon(), Scalar(3));
    }
};

// Full specialization
template <>
struct TypeTraits<float>
{
    static inline float min()
    {
        return std::numeric_limits<float>::min();
    }
};

template <>
struct TypeTraits<double>
{
    static inline double min()
    {
        return std::numeric_limits<double>::min();
    }
};

template <>
struct TypeTraits<long double>
{
    static inline long double min()
    {
        return std::numeric_limits<long double>::min();
    }
};


} // namespace Spectra

/// \endcond

#endif // TYPE_TRAITS_H
