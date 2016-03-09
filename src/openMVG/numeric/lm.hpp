// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_NUMERIC_LM_H_
#define OPENMVG_NUMERIC_LM_H_

#include "openMVG/numeric/numeric.h"
// Levenberg Marquardt Non Linear Optimization
#include <unsupported/Eigen/NonLinearOptimization>
#include <unsupported/Eigen/NumericalDiff>

namespace openMVG
{
using namespace Eigen;

/**
 * @brief Generic functor Levenberg-Marquardt minimization
 * @tparam _Scalar Type of internal computation
 * @tparam NX Number of values per sample at compile time
 * @tparam NY Number of samples at compile time
 */
template<typename _Scalar, int NX = Dynamic, int NY = Dynamic>
struct Functor
{
  typedef _Scalar Scalar;
  enum
  {
    InputsAtCompileTime = NX,
    ValuesAtCompileTime = NY
  };
  typedef Matrix<Scalar, InputsAtCompileTime, 1> InputType;
  typedef Matrix<Scalar, ValuesAtCompileTime, 1> ValueType;
  typedef Matrix<Scalar, ValuesAtCompileTime, InputsAtCompileTime> JacobianType;


  /// Number of values per sample
  const int m_inputs ;

  // Number of sample
  const int m_values;

  /**
  * @brief Default constructor
  */
  Functor()
    : m_inputs( InputsAtCompileTime ),
      m_values( ValuesAtCompileTime )
  {

  }

  /**
  * @brief Constructor
  * @param inputs Number of column per sample
  * @param values Number of sample
  */
  Functor( int inputs, int values ) : m_inputs( inputs ), m_values( values ) {}

  /**
  * @brief Get number of samples
  * @return Number of samples
  */
  int inputs() const
  {
    return m_inputs;
  }

  /**
  * @brief Get number of samples
  * @return Number of samples
  */
  int values() const
  {
    return m_values;
  }

  // you should define that in the subclass :
  //  void operator() (const InputType& x, ValueType* v, JacobianType* _j=0) const;
};

}; // namespace openMVG

#endif // OPENMVG_NUMERIC_LM_H_
