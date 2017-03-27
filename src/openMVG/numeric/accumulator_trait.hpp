// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_NUMERIC_ACCUMULATOR_TRAIT_HPP
#define OPENMVG_NUMERIC_ACCUMULATOR_TRAIT_HPP

/// Accumulator trait to perform safe summation over a specified type
namespace openMVG {

template<typename T>
struct Accumulator { using Type = T; };
template<>
struct Accumulator<unsigned char>  { using Type = float; };
template<>
struct Accumulator<unsigned short> { using Type = float; };
template<>
struct Accumulator<unsigned int> { using Type = float; };
template<>
struct Accumulator<char>   { using Type = float; };
template<>
struct Accumulator<short>  { using Type = float; };
template<>
struct Accumulator<int> { using Type = float; };
template<>
struct Accumulator<bool>  { using Type = unsigned int; };

} // namespace openMVG

#endif //OPENMVG_NUMERIC_ACCUMULATOR_TRAIT_HPP
