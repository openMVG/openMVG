// Copyright (C) 2016-2019 Yixuan Qiu <yixuan.qiu@cos.name>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef GEIGS_MODE_H
#define GEIGS_MODE_H

namespace Spectra {


///
/// \ingroup Enumerations
///
/// The enumeration to specify the mode of generalized eigenvalue solver.
///
enum GEIGS_MODE
{
    GEIGS_CHOLESKY = 0,     ///< Using Cholesky decomposition to solve generalized eigenvalues.

    GEIGS_REGULAR_INVERSE,  ///< Regular inverse mode for generalized eigenvalue solver.

    GEIGS_SHIFT_INVERT,     ///< Shift-and-invert mode for generalized eigenvalue solver.

    GEIGS_BUCKLING,         ///< Buckling mode for generalized eigenvalue solver.

    GEIGS_CAYLEY            ///< Cayley transformation mode for generalized eigenvalue solver.
};


} // namespace Spectra

#endif // GEIGS_MODE_H
