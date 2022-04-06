// Copyright (C) 2016-2022 Yixuan Qiu <yixuan.qiu@cos.name>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SPECTRA_SELECTION_RULE_H
#define SPECTRA_SELECTION_RULE_H

#include <vector>     // std::vector
#include <cmath>      // std::abs
#include <algorithm>  // std::sort
#include <complex>    // std::complex
#include <utility>    // std::pair
#include <stdexcept>  // std::invalid_argument

#include <Eigen/Core>
#include "TypeTraits.h"

namespace Spectra {

///
/// \defgroup Enumerations Enumerations
///
/// Enumeration types for the selection rule of eigenvalues.
///

///
/// \ingroup Enumerations
///
/// The enumeration of selection rules of desired eigenvalues.
///
enum class SortRule
{
    LargestMagn,  ///< Select eigenvalues with largest magnitude. Magnitude
                  ///< means the absolute value for real numbers and norm for
                  ///< complex numbers. Applies to both symmetric and general
                  ///< eigen solvers.

    LargestReal,  ///< Select eigenvalues with largest real part. Only for general eigen solvers.

    LargestImag,  ///< Select eigenvalues with largest imaginary part (in magnitude). Only for general eigen solvers.

    LargestAlge,  ///< Select eigenvalues with largest algebraic value, considering
                  ///< any negative sign. Only for symmetric eigen solvers.

    SmallestMagn,  ///< Select eigenvalues with smallest magnitude. Applies to both symmetric and general
                   ///< eigen solvers.

    SmallestReal,  ///< Select eigenvalues with smallest real part. Only for general eigen solvers.

    SmallestImag,  ///< Select eigenvalues with smallest imaginary part (in magnitude). Only for general eigen solvers.

    SmallestAlge,  ///< Select eigenvalues with smallest algebraic value. Only for symmetric eigen solvers.

    BothEnds  ///< Select eigenvalues half from each end of the spectrum. When
              ///< `nev` is odd, compute more from the high end. Only for symmetric eigen solvers.
};

/// \cond

// When comparing eigenvalues, we first calculate the "target" to sort.
// For example, if we want to choose the eigenvalues with
// largest magnitude, the target will be -abs(x).
// The minus sign is due to the fact that std::sort() sorts in ascending order.

// Default target: throw an exception
template <typename Scalar, SortRule Rule>
class SortingTarget
{
public:
    static ElemType<Scalar> get(const Scalar& val)
    {
        using std::abs;
        throw std::invalid_argument("incompatible selection rule");
        return -abs(val);
    }
};

// Specialization for SortRule::LargestMagn
// This covers [float, double, complex] x [SortRule::LargestMagn]
template <typename Scalar>
class SortingTarget<Scalar, SortRule::LargestMagn>
{
public:
    static ElemType<Scalar> get(const Scalar& val)
    {
        using std::abs;
        return -abs(val);
    }
};

// Specialization for SortRule::LargestReal
// This covers [complex] x [SortRule::LargestReal]
template <typename RealType>
class SortingTarget<std::complex<RealType>, SortRule::LargestReal>
{
public:
    static RealType get(const std::complex<RealType>& val)
    {
        return -val.real();
    }
};

// Specialization for SortRule::LargestImag
// This covers [complex] x [SortRule::LargestImag]
template <typename RealType>
class SortingTarget<std::complex<RealType>, SortRule::LargestImag>
{
public:
    static RealType get(const std::complex<RealType>& val)
    {
        using std::abs;
        return -abs(val.imag());
    }
};

// Specialization for SortRule::LargestAlge
// This covers [float, double] x [SortRule::LargestAlge]
template <typename Scalar>
class SortingTarget<Scalar, SortRule::LargestAlge>
{
public:
    static Scalar get(const Scalar& val)
    {
        return -val;
    }
};

// Here SortRule::BothEnds is the same as SortRule::LargestAlge, but
// we need some additional steps, which are done in
// SymEigsSolver.h => retrieve_ritzpair().
// There we move the smallest values to the proper locations.
template <typename Scalar>
class SortingTarget<Scalar, SortRule::BothEnds>
{
public:
    static Scalar get(const Scalar& val)
    {
        return -val;
    }
};

// Specialization for SortRule::SmallestMagn
// This covers [float, double, complex] x [SortRule::SmallestMagn]
template <typename Scalar>
class SortingTarget<Scalar, SortRule::SmallestMagn>
{
public:
    static ElemType<Scalar> get(const Scalar& val)
    {
        using std::abs;
        return abs(val);
    }
};

// Specialization for SortRule::SmallestReal
// This covers [complex] x [SortRule::SmallestReal]
template <typename RealType>
class SortingTarget<std::complex<RealType>, SortRule::SmallestReal>
{
public:
    static RealType get(const std::complex<RealType>& val)
    {
        return val.real();
    }
};

// Specialization for SortRule::SmallestImag
// This covers [complex] x [SortRule::SmallestImag]
template <typename RealType>
class SortingTarget<std::complex<RealType>, SortRule::SmallestImag>
{
public:
    static RealType get(const std::complex<RealType>& val)
    {
        using std::abs;
        return abs(val.imag());
    }
};

// Specialization for SortRule::SmallestAlge
// This covers [float, double] x [SortRule::SmallestAlge]
template <typename Scalar>
class SortingTarget<Scalar, SortRule::SmallestAlge>
{
public:
    static Scalar get(const Scalar& val)
    {
        return val;
    }
};

// Sort eigenvalues
template <typename T, SortRule Rule>
class SortEigenvalue
{
private:
    using Index = Eigen::Index;
    using IndexArray = std::vector<Index>;

    const T* m_evals;
    IndexArray m_index;

public:
    // Sort indices according to the eigenvalues they point to
    inline bool operator()(Index i, Index j)
    {
        return SortingTarget<T, Rule>::get(m_evals[i]) < SortingTarget<T, Rule>::get(m_evals[j]);
    }

    SortEigenvalue(const T* start, Index size) :
        m_evals(start), m_index(size)
    {
        for (Index i = 0; i < size; i++)
        {
            m_index[i] = i;
        }
        std::sort(m_index.begin(), m_index.end(), *this);
    }

    inline IndexArray index() const { return m_index; }
    inline void swap(IndexArray& other) { m_index.swap(other); }
};

// Sort values[:len] according to the selection rule, and return the indices
template <typename Scalar>
std::vector<Eigen::Index> argsort(SortRule selection, const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& values, Eigen::Index len)
{
    using Index = Eigen::Index;

    // Sort Ritz values and put the wanted ones at the beginning
    std::vector<Index> ind;
    switch (selection)
    {
        case SortRule::LargestMagn:
        {
            SortEigenvalue<Scalar, SortRule::LargestMagn> sorting(values.data(), len);
            sorting.swap(ind);
            break;
        }
        case SortRule::BothEnds:
        case SortRule::LargestAlge:
        {
            SortEigenvalue<Scalar, SortRule::LargestAlge> sorting(values.data(), len);
            sorting.swap(ind);
            break;
        }
        case SortRule::SmallestMagn:
        {
            SortEigenvalue<Scalar, SortRule::SmallestMagn> sorting(values.data(), len);
            sorting.swap(ind);
            break;
        }
        case SortRule::SmallestAlge:
        {
            SortEigenvalue<Scalar, SortRule::SmallestAlge> sorting(values.data(), len);
            sorting.swap(ind);
            break;
        }
        default:
            throw std::invalid_argument("unsupported selection rule");
    }

    // For SortRule::BothEnds, the eigenvalues are sorted according to the
    // SortRule::LargestAlge rule, so we need to move those smallest values to the left
    // The order would be
    //     Largest => Smallest => 2nd largest => 2nd smallest => ...
    // We keep this order since the first k values will always be
    // the wanted collection, no matter k is nev_updated (used in SymEigsBase::restart())
    // or is nev (used in SymEigsBase::sort_ritzpair())
    if (selection == SortRule::BothEnds)
    {
        std::vector<Index> ind_copy(ind);
        for (Index i = 0; i < len; i++)
        {
            // If i is even, pick values from the left (large values)
            // If i is odd, pick values from the right (small values)
            if (i % 2 == 0)
                ind[i] = ind_copy[i / 2];
            else
                ind[i] = ind_copy[len - 1 - i / 2];
        }
    }

    return ind;
}

// Default vector length
template <typename Scalar>
std::vector<Eigen::Index> argsort(SortRule selection, const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& values)
{
    return argsort<Scalar>(selection, values, values.size());
}

/// \endcond

}  // namespace Spectra

#endif  // SPECTRA_SELECTION_RULE_H
