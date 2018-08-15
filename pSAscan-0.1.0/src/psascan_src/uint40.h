/******************************************************************************
 *
 * Class representing a 40-bit unsigned integer encoded in five bytes.
 *
 ******************************************************************************
 * Copyright (C) 2012 Timo Bingmann <tb@panthema.net>
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 *****************************************************************************
 *
 * NOTE: This is slightly modified version of the file used in eSAIS-0.5.4
 * (https://panthema.net/2012/1119-eSAIS-Inducing-Suffix-and-LCP-Arrays-
 * in-External-Memory/). In particular, it contains a small bugfix in the
 * += operator.
 *
 * Modified by Dominik Kempa <dominik.kempa (at) gmail.com>
 *
 *****************************************************************************/


#ifndef __UINT40_H_INCLUDED
#define __UINT40_H_INCLUDED

#include <inttypes.h>
#include <stdint.h>
#include <cassert>
#include <iostream>
#include <limits>
#include <unistd.h>


class uint40
{
private:
    uint32_t    low;
    uint8_t     high;

public:
    inline uint40()
    {
    }

    inline uint40(uint32_t l, uint8_t h)
        : low(l), high(h)
    {
    }

    inline uint40(const uint40& a)
        : low(a.low), high(a.high)
    {
    }

    inline uint40(const int& a)
        : low(a), high(0)
    {
    }

    inline uint40(const unsigned int& a)
      : low(a), high(0)
    {
    }

    inline uint40(const uint64_t& a)
        : low(a & 0xFFFFFFFF), high((a >> 32) & 0xFF)
    {
        assert( a <= 0xFFFFFFFFFFLU );
    }

    inline uint40(const long& a)
      : low(a & 0xFFFFFFFFL), high((a >> 32) & 0xFF) {
      assert( a <= 0xFFFFFFFFFFL );
    }

    inline uint64_t ull() const {
        return ((uint64_t)high) << 32 | (uint64_t)low;
    }

    inline long ll() const
    {
        return (long)ull();
    }

    inline operator uint64_t() const
    {
        return ull();
    }

    inline uint64_t u64() const
    {
        return ((uint64_t)high) << 32 | (uint64_t)low;
    }

    inline uint40& operator++ ()
    {
        if (low == std::numeric_limits<uint32_t>::max())
            ++high, low = 0;
        else
            ++low;
        return *this;
    }

    inline uint40& operator-- ()
    {
        if (low == 0)
            --high, low = std::numeric_limits<uint32_t>::max();
        else
            --low;
        return *this;
    }

    inline uint40& operator+= (const uint40& b)
    {
        uint64_t add = (uint64_t)low + b.low;  // BUGFIX
        low = add & 0xFFFFFFFF;
        high += b.high + ((add >> 32) & 0xFF);
        return *this;
    }

    inline bool operator== (const uint40& b) const
    {
        return (low == b.low) && (high == b.high);
    }

    inline bool operator!= (const uint40& b) const
    {
        return (low != b.low) || (high != b.high);
    }

    inline bool operator< (const uint40& b) const
    {
        return (high < b.high) || (high == b.high && low < b.low);
    }

    inline bool operator<= (const uint40& b) const
    {
        return (high < b.high) || (high == b.high && low <= b.low);
    }

    inline bool operator> (const uint40& b) const
    {
        return (high > b.high) || (high == b.high && low > b.low);
    }

    inline bool operator>= (const uint40& b) const
    {
        return (high > b.high) || (high == b.high && low >= b.low);
    }

    friend std::ostream& operator<< (std::ostream& os, const uint40& a)
    {
        return os << a.ull();
    }

} __attribute__((packed));

namespace std {

template<>
class numeric_limits<uint40> {
public:
    static uint40 min() { return uint40(std::numeric_limits<uint32_t>::min(),
                                        std::numeric_limits<uint8_t>::min()); }

    static uint40 max() { return uint40(std::numeric_limits<uint32_t>::max(),
                                        std::numeric_limits<uint8_t>::max()); }
};

}

#endif  // __UINT40_H_INCLUDED
