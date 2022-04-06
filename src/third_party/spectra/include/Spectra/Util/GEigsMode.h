// Copyright (C) 2016-2022 Yixuan Qiu <yixuan.qiu@cos.name>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SPECTRA_GEIGS_MODE_H
#define SPECTRA_GEIGS_MODE_H

namespace Spectra {

///
/// \ingroup Enumerations
///
/// The enumeration to specify the mode of generalized eigenvalue solver.
///
enum class GEigsMode
{
    Cholesky,        ///< Using Cholesky decomposition to solve generalized eigenvalues.
    RegularInverse,  ///< Regular inverse mode for generalized eigenvalue solver.
    ShiftInvert,     ///< Shift-and-invert mode for generalized eigenvalue solver.
    Buckling,        ///< Buckling mode for generalized eigenvalue solver.
    Cayley           ///< Cayley transformation mode for generalized eigenvalue solver.
};

}  // namespace Spectra

#endif  // SPECTRA_GEIGS_MODE_H
