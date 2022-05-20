// Copyright (C) 2016-2022 Yixuan Qiu <yixuan.qiu@cos.name>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SPECTRA_SIMPLE_RANDOM_H
#define SPECTRA_SIMPLE_RANDOM_H

#include <Eigen/Core>

/// \cond

namespace Spectra {

// We need a simple pseudo random number generator here:
// 1. It is used to generate initial and restarted residual vector.
// 2. It is not necessary to be so "random" and advanced. All we hope
//    is that the residual vector is not in the space spanned by the
//    current Krylov space. This should be met almost surely.
// 3. We don't want to call RNG in C++, since we actually want the
//    algorithm to be deterministic. Also, calling RNG in C/C++ is not
//    allowed in R packages submitted to CRAN.
// 4. The method should be as simple as possible, so an LCG is enough.
// 5. Based on public domain code by Ray Gardner
//    http://stjarnhimlen.se/snippets/rg_rand.c

template <typename Scalar = double>
class SimpleRandom
{
private:
    using Index = Eigen::Index;
    using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

    static constexpr unsigned int m_a = 16807;           // multiplier
    static constexpr unsigned long m_max = 2147483647L;  // 2^31 - 1
    long m_rand;                                         // RNG state

    inline long next_long_rand(long seed) const
    {
        unsigned long lo, hi;

        lo = m_a * (long) (seed & 0xFFFF);
        hi = m_a * (long) ((unsigned long) seed >> 16);
        lo += (hi & 0x7FFF) << 16;
        if (lo > m_max)
        {
            lo &= m_max;
            ++lo;
        }
        lo += hi >> 15;
        if (lo > m_max)
        {
            lo &= m_max;
            ++lo;
        }
        return (long) lo;
    }

public:
    SimpleRandom(unsigned long init_seed) :
        m_rand(init_seed ? (init_seed & m_max) : 1)
    {}

    // Return a single random number, ranging from -0.5 to 0.5
    Scalar random()
    {
        m_rand = next_long_rand(m_rand);
        return Scalar(m_rand) / Scalar(m_max) - Scalar(0.5);
    }

    // Fill the given vector with random numbers
    // Ranging from -0.5 to 0.5
    void random_vec(Vector& vec)
    {
        const Index len = vec.size();
        for (Index i = 0; i < len; i++)
        {
            m_rand = next_long_rand(m_rand);
            vec[i] = Scalar(m_rand);
        }
        vec.array() = vec.array() / Scalar(m_max) - Scalar(0.5);
    }

    // Return a vector of random numbers
    // Ranging from -0.5 to 0.5
    Vector random_vec(const Index len)
    {
        Vector res(len);
        random_vec(res);
        return res;
    }
};

}  // namespace Spectra

/// \endcond

#endif  // SPECTRA_SIMPLE_RANDOM_H
