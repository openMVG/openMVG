// Copyright (C) 2016-2019 Yixuan Qiu <yixuan.qiu@cos.name>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SELECTION_RULE_H
#define SELECTION_RULE_H

#include <vector>     // std::vector
#include <cmath>      // std::abs
#include <algorithm>  // std::sort
#include <complex>    // std::complex
#include <utility>    // std::pair
#include <stdexcept>  // std::invalid_argument

namespace Spectra {


///
/// \defgroup Enumerations
///
/// Enumeration types for the selection rule of eigenvalues.
///

///
/// \ingroup Enumerations
///
/// The enumeration of selection rules of desired eigenvalues.
///
enum SELECT_EIGENVALUE
{
    LARGEST_MAGN = 0,  ///< Select eigenvalues with largest magnitude. Magnitude
                       ///< means the absolute value for real numbers and norm for
                       ///< complex numbers. Applies to both symmetric and general
                       ///< eigen solvers.

    LARGEST_REAL,      ///< Select eigenvalues with largest real part. Only for general eigen solvers.

    LARGEST_IMAG,      ///< Select eigenvalues with largest imaginary part (in magnitude). Only for general eigen solvers.

    LARGEST_ALGE,      ///< Select eigenvalues with largest algebraic value, considering
                       ///< any negative sign. Only for symmetric eigen solvers.

    SMALLEST_MAGN,     ///< Select eigenvalues with smallest magnitude. Applies to both symmetric and general
                       ///< eigen solvers.

    SMALLEST_REAL,     ///< Select eigenvalues with smallest real part. Only for general eigen solvers.

    SMALLEST_IMAG,     ///< Select eigenvalues with smallest imaginary part (in magnitude). Only for general eigen solvers.

    SMALLEST_ALGE,     ///< Select eigenvalues with smallest algebraic value. Only for symmetric eigen solvers.

    BOTH_ENDS          ///< Select eigenvalues half from each end of the spectrum. When
                       ///< `nev` is odd, compute more from the high end. Only for symmetric eigen solvers.
};

///
/// \ingroup Enumerations
///
/// The enumeration of selection rules of desired eigenvalues. Alias for `SELECT_EIGENVALUE`.
///
enum SELECT_EIGENVALUE_ALIAS
{
    WHICH_LM = 0,  ///< Alias for `LARGEST_MAGN`
    WHICH_LR,      ///< Alias for `LARGEST_REAL`
    WHICH_LI,      ///< Alias for `LARGEST_IMAG`
    WHICH_LA,      ///< Alias for `LARGEST_ALGE`
    WHICH_SM,      ///< Alias for `SMALLEST_MAGN`
    WHICH_SR,      ///< Alias for `SMALLEST_REAL`
    WHICH_SI,      ///< Alias for `SMALLEST_IMAG`
    WHICH_SA,      ///< Alias for `SMALLEST_ALGE`
    WHICH_BE       ///< Alias for `BOTH_ENDS`
};

/// \cond

// Get the element type of a "scalar"
// ElemType<double>                   => double
// ElemType< std::complex<double> >   => double
template <typename T>
class ElemType
{
public:
    typedef T type;
};

template <typename T>
class ElemType< std::complex<T> >
{
public:
    typedef T type;
};

// When comparing eigenvalues, we first calculate the "target"
// to sort. For example, if we want to choose the eigenvalues with
// largest magnitude, the target will be -abs(x).
// The minus sign is due to the fact that std::sort() sorts in ascending order.

// Default target: throw an exception
template <typename Scalar, int SelectionRule>
class SortingTarget
{
public:
    static typename ElemType<Scalar>::type get(const Scalar& val)
    {
        using std::abs;
        throw std::invalid_argument("incompatible selection rule");
        return -abs(val);
    }
};

// Specialization for LARGEST_MAGN
// This covers [float, double, complex] x [LARGEST_MAGN]
template <typename Scalar>
class SortingTarget<Scalar, LARGEST_MAGN>
{
public:
    static typename ElemType<Scalar>::type get(const Scalar& val)
    {
        using std::abs;
        return -abs(val);
    }
};

// Specialization for LARGEST_REAL
// This covers [complex] x [LARGEST_REAL]
template <typename RealType>
class SortingTarget<std::complex<RealType>, LARGEST_REAL>
{
public:
    static RealType get(const std::complex<RealType>& val)
    {
        return -val.real();
    }
};

// Specialization for LARGEST_IMAG
// This covers [complex] x [LARGEST_IMAG]
template <typename RealType>
class SortingTarget<std::complex<RealType>, LARGEST_IMAG>
{
public:
    static RealType get(const std::complex<RealType>& val)
    {
        using std::abs;
        return -abs(val.imag());
    }
};

// Specialization for LARGEST_ALGE
// This covers [float, double] x [LARGEST_ALGE]
template <typename Scalar>
class SortingTarget<Scalar, LARGEST_ALGE>
{
public:
    static Scalar get(const Scalar& val)
    {
        return -val;
    }
};

// Here BOTH_ENDS is the same as LARGEST_ALGE, but
// we need some additional steps, which are done in
// SymEigsSolver.h => retrieve_ritzpair().
// There we move the smallest values to the proper locations.
template <typename Scalar>
class SortingTarget<Scalar, BOTH_ENDS>
{
public:
    static Scalar get(const Scalar& val)
    {
        return -val;
    }
};

// Specialization for SMALLEST_MAGN
// This covers [float, double, complex] x [SMALLEST_MAGN]
template <typename Scalar>
class SortingTarget<Scalar, SMALLEST_MAGN>
{
public:
    static typename ElemType<Scalar>::type get(const Scalar& val)
    {
        using std::abs;
        return abs(val);
    }
};

// Specialization for SMALLEST_REAL
// This covers [complex] x [SMALLEST_REAL]
template <typename RealType>
class SortingTarget<std::complex<RealType>, SMALLEST_REAL>
{
public:
    static RealType get(const std::complex<RealType>& val)
    {
        return val.real();
    }
};

// Specialization for SMALLEST_IMAG
// This covers [complex] x [SMALLEST_IMAG]
template <typename RealType>
class SortingTarget<std::complex<RealType>, SMALLEST_IMAG>
{
public:
    static RealType get(const std::complex<RealType>& val)
    {
        using std::abs;
        return abs(val.imag());
    }
};

// Specialization for SMALLEST_ALGE
// This covers [float, double] x [SMALLEST_ALGE]
template <typename Scalar>
class SortingTarget<Scalar, SMALLEST_ALGE>
{
public:
    static Scalar get(const Scalar& val)
    {
        return val;
    }
};

// Sort eigenvalues and return the order index
template <typename PairType>
class PairComparator
{
public:
    bool operator() (const PairType& v1, const PairType& v2)
    {
        return v1.first < v2.first;
    }
};

template <typename T, int SelectionRule>
class SortEigenvalue
{
private:
    typedef typename ElemType<T>::type TargetType; // Type of the sorting target, will be
                                                   // a floating number type, e.g. "double"
    typedef std::pair<TargetType, int> PairType;   // Type of the sorting pair, including
                                                   // the sorting target and the index

    std::vector<PairType> pair_sort;

public:
    SortEigenvalue(const T* start, int size) :
        pair_sort(size)
    {
        for(int i = 0; i < size; i++)
        {
            pair_sort[i].first = SortingTarget<T, SelectionRule>::get(start[i]);
            pair_sort[i].second = i;
        }
        PairComparator<PairType> comp;
        std::sort(pair_sort.begin(), pair_sort.end(), comp);
    }

    std::vector<int> index()
    {
        std::vector<int> ind(pair_sort.size());
        for(unsigned int i = 0; i < ind.size(); i++)
            ind[i] = pair_sort[i].second;

        return ind;
    }
};

/// \endcond


} // namespace Spectra

#endif // SELECTION_RULE_H
