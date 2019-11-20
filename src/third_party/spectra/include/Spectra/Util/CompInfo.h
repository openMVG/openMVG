// Copyright (C) 2016-2019 Yixuan Qiu <yixuan.qiu@cos.name>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef COMP_INFO_H
#define COMP_INFO_H

namespace Spectra {


///
/// \ingroup Enumerations
///
/// The enumeration to report the status of computation.
///
enum COMPUTATION_INFO
{
    SUCCESSFUL = 0,    ///< Computation was successful.

    NOT_COMPUTED,      ///< Used in eigen solvers, indicating that computation
                       ///< has not been conducted. Users should call
                       ///< the `compute()` member function of solvers.

    NOT_CONVERGING,    ///< Used in eigen solvers, indicating that some eigenvalues
                       ///< did not converge. The `compute()`
                       ///< function returns the number of converged eigenvalues.

    NUMERICAL_ISSUE    ///< Used in Cholesky decomposition, indicating that the
                       ///< matrix is not positive definite.
};


} // namespace Spectra

#endif // COMP_INFO_H
