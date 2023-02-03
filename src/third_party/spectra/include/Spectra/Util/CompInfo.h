// Copyright (C) 2016-2022 Yixuan Qiu <yixuan.qiu@cos.name>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SPECTRA_COMP_INFO_H
#define SPECTRA_COMP_INFO_H

namespace Spectra {

///
/// \ingroup Enumerations
///
/// The enumeration to report the status of computation.
///
enum class CompInfo
{
    Successful,  ///< Computation was successful.

    NotComputed,  ///< Used in eigen solvers, indicating that computation
                  ///< has not been conducted. Users should call
                  ///< the `compute()` member function of solvers.

    NotConverging,  ///< Used in eigen solvers, indicating that some eigenvalues
                    ///< did not converge. The `compute()`
                    ///< function returns the number of converged eigenvalues.

    NumericalIssue  ///< Used in various matrix factorization classes, for example in
                    ///< Cholesky decomposition it indicates that the
                    ///< matrix is not positive definite.
};

}  // namespace Spectra

#endif  // SPECTRA_COMP_INFO_H
