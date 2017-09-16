// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// Class derived from the Boost dynamic_bitset<> implementation.

// Distributed under the Boost Software License, Version 1.0.
//    (http://www.boost.org/LICENSE_1_0.txt)
//
// -----------------------------------------------------------

#ifndef OPENMVG_STL_DYNAMIC_BITSET_HPP
#define OPENMVG_STL_DYNAMIC_BITSET_HPP

#include <cassert>
#include <cstddef>
#include <limits>
#include <vector>

namespace stl
{

  struct dynamic_bitset
  {
    using BlockType =  unsigned char;
    static const int bits_per_block = (std::numeric_limits<BlockType>::digits);

    // Reference over a single bit of a sub-block
    class reference
    {
      BlockType & m_block;
      const BlockType m_mask;

    public:
      reference
      (
        BlockType & b,
        BlockType pos
      )
      : m_block(b),
        m_mask((assert(pos < bits_per_block), BlockType(1) << pos))
      { }

      operator bool() const { return (m_block & m_mask) != 0; }
      reference& operator=(bool x) { do_assign(x); return *this; } // for b[i] = x

      void do_set()   { m_block |= m_mask; }
      void do_reset() { m_block &= ~m_mask;}
      void do_flip()  { m_block ^= m_mask; }
      void do_assign(bool x) { x ? do_set() : do_reset(); }
    };

    size_t block_index(size_t pos) const
    {
      return pos / bits_per_block;
    }

    size_t bit_index(size_t pos) const
    {
      return static_cast<size_t>(pos % bits_per_block);
    }

    BlockType bit_mask(size_t pos) const
    {
      return BlockType(1) << bit_index(pos);
    }

    bool operator[](size_t pos) const
    {
      return (vec_bits[block_index(pos)] & bit_mask(pos)) != 0;
    }

    reference operator[](size_t pos)
    {
      return reference(vec_bits[block_index(pos)], bit_index(pos));
    }

    // return the number of stored bits
    size_t size() const { return m_num_bits; }

    // return the number of BlockType block used to manipulate the bits
    size_t num_blocks() const { return vec_bits.size(); }

    // Constructor
    dynamic_bitset(size_t num_bits = 0)
    {
      vec_bits.resize(calc_num_blocks(num_bits));
      m_num_bits = num_bits;
    }

    dynamic_bitset& reset()
    {
      std::fill(vec_bits.begin(), vec_bits.end(), BlockType(0));
      return *this;
    }

    const BlockType * data() const { return &vec_bits[0]; }

  private:
    inline size_t calc_num_blocks(size_t num_bits)
    {
      return num_bits / bits_per_block
        + static_cast<size_t>(num_bits % bits_per_block != 0);
    }

    // DATA
    std::vector<BlockType> vec_bits;
    size_t m_num_bits;
  };
} // namespace stl

#endif // OPENMVG_STL_DYNAMIC_BITSET_HPP
