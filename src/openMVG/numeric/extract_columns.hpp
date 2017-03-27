// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_NUMERIC_EXTRACT_COLUMNS_HPP
#define OPENMVG_NUMERIC_EXTRACT_COLUMNS_HPP

#include "openMVG/numeric/eigen_alias_definition.hpp"

namespace openMVG
{

/**
* @brief Extract a submatrix given a list of column
* @param A Input matrix
* @param columns A vector of columns index to extract
* @return Matrix containing a subset of input matrix columns
* @note columns index start at index 0
* @note Assuming columns contains a list of valid columns index
*/
template <typename TCols>
Mat ExtractColumns
(
  const Eigen::Ref<const Mat> &A,
  const TCols &columns
)
{
  Mat compressed( A.rows(), columns.size() );
  for ( size_t i = 0; i < static_cast<size_t>( columns.size() ); ++i )
  {
    compressed.col( i ) = A.col( columns[i] );
  }
  return compressed;
}

} // namespace openMVG

#endif  // OPENMVG_NUMERIC_EXTRACT_COLUMNS_HPP
