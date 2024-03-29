/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018-2019 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline Foam::label Foam::bitSet::first_not_block() const
{
    if (empty())
    {
        return -1;
    }

    // Use complement to change 0 <-> 1 and check if any 1's now appear

    const label nblocks = num_blocks(size());

    // Extra bits in the final block?
    const unsigned int off = size() % elem_per_block;

    if (!off)
    {
        for (label blocki=0; blocki < nblocks; ++blocki)
        {
            if (~(blocks_[blocki]))
            {
                return blocki;
            }
        }
    }
    else
    {
        for (label blocki=0; blocki < nblocks-1; ++blocki)
        {
            if (~(blocks_[blocki]))
            {
                return blocki;
            }
        }

        // The final block needs masking
        if (~(blocks_[nblocks-1]) & mask_lower(off))
        {
            return nblocks-1;
        }
    }

    return -1;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::bitSet::bitSet() noexcept
:
    PackedList<1>()
{}


inline Foam::bitSet::bitSet(const label n)
:
    PackedList<1>(n)
{}


inline Foam::bitSet::bitSet(const label n, const bool val)
:
    bitSet(n)
{
    if (val) assign(val);
}


inline Foam::bitSet::bitSet(const bitSet& bitset)
:
    PackedList<1>(bitset)
{}


inline Foam::bitSet::bitSet(bitSet&& bitset)
:
    PackedList<1>(std::move(bitset))
{}


inline Foam::bitSet::bitSet(const UList<bool>& bools)
:
    bitSet()
{
    assign(bools);
}


inline Foam::bitSet::bitSet(const label n, const labelUList& locations)
:
    bitSet(n)
{
    setMany(locations.begin(), locations.end());
}


template<class Addr>
inline Foam::bitSet::bitSet
(
    const label n,
    const IndirectListBase<label, Addr>& locations
)
:
    bitSet(n)
{
    setMany(locations.begin(), locations.end());
}


inline Foam::bitSet::bitSet
(
    const label n,
    std::initializer_list<label> locations
)
:
    bitSet(n)
{
    setMany(locations.begin(), locations.end());
}


inline Foam::bitSet::bitSet(const labelUList& locations)
:
    bitSet()
{
    setMany(locations.begin(), locations.end());
}


template<class Addr>
inline Foam::bitSet::bitSet
(
    const IndirectListBase<label, Addr>& locations
)
:
    bitSet()
{
    setMany(locations.begin(), locations.end());
}


inline Foam::autoPtr<Foam::bitSet> Foam::bitSet::clone() const
{
    return autoPtr<bitSet>::New(*this);
}


// * * * * * * * * * * * * * * * * References * * * * * * * * * * * * * * * * //

inline Foam::bitSet::reference::reference
(
    bitSet* parent,
    const label index
)
:
    PackedList<1>::reference(parent, index)
{}


inline void Foam::bitSet::reference::flip()
{
    const unsigned int mask = (max_value << shift_);
    ref_ ^= mask;
}


inline void Foam::bitSet::reference::operator=
(
    const reference& other
)
{
    set(other.get());
}


inline void Foam::bitSet::reference::operator=
(
    const unsigned int val
)
{
    set(val);
}


inline Foam::bitSet::reference::operator unsigned int () const
{
    return get();
}


// * * * * * * * * * * * * * * * * Iterators * * * * * * * * * * * * * * * * //

inline Foam::bitSet::const_iterator::const_iterator() noexcept
:
    set_(nullptr),
    pos_(-1)
{}


inline Foam::bitSet::const_iterator::const_iterator(const bitSet* parent)
:
    set_(parent),
    pos_(set_->find_first())
{}


inline Foam::label Foam::bitSet::const_iterator::operator*() const noexcept
{
    return pos_;
}


inline Foam::bitSet::const_iterator& Foam::bitSet::const_iterator::operator++()
{
    pos_ = set_->find_next(pos_);
    return *this;
}


inline bool Foam::bitSet::const_iterator::operator==
(
    const const_iterator& iter
) const noexcept
{
    return (iter.pos_ == pos_);
}


inline bool Foam::bitSet::const_iterator::operator!=
(
    const const_iterator& iter
) const noexcept
{
    return (iter.pos_ != pos_);
}


inline Foam::bitSet::const_iterator Foam::bitSet::begin() const
{
    return const_iterator(this);
}


inline Foam::bitSet::const_iterator Foam::bitSet::cbegin() const
{
    return const_iterator(this);
}


inline Foam::bitSet::const_iterator Foam::bitSet::end() const noexcept
{
    return const_iterator();
}


inline Foam::bitSet::const_iterator Foam::bitSet::cend() const noexcept
{
    return const_iterator();
}


inline Foam::label Foam::bitSet::find_first() const
{
    // Process block-wise, detecting any '1' bits

    const label nblocks = num_blocks(size());
    for (label blocki = 0; blocki < nblocks; ++blocki)
    {
        label pos = (blocki * elem_per_block);

        for
        (
            unsigned int blockval = blocks_[blocki];
            blockval;
            blockval >>= 1u
        )
        {
            if (blockval & 1u)
            {
                return pos;
            }
            ++pos;
        }
    }

    return -1;
}


inline Foam::label Foam::bitSet::find_first_not() const
{
    const label blocki = first_not_block();

    if (blocki >= 0)
    {
        label pos = (blocki * elem_per_block);

        // Detect first '0' bit by checking the complement.

        // No special masking for the final block, that was already checked
        // in the first_not_block() call.

        for
        (
            unsigned int blockval = ~(blocks_[blocki]);
            blockval;
            blockval >>= 1u
        )
        {
            if (blockval & 1u)
            {
                return pos;
            }
            ++pos;
        }
    }

    return -1;
}


inline Foam::label Foam::bitSet::find_last() const
{
    // Process block-wise, detecting any '1' bits

    for (label blocki = num_blocks(size())-1; blocki >= 0; --blocki)
    {
        unsigned int blockval = blocks_[blocki];

        if (blockval)
        {
            label pos = (blocki * elem_per_block) - 1;

            while (blockval)
            {
                blockval >>= 1u;
                ++pos;
            }

            return pos;
        }
    }

    return -1;
}


inline Foam::label Foam::bitSet::find_next(label pos) const
{
    ++pos;
    if (pos < 0 || pos >= size())
    {
        return -1;
    }

    // The corresponding block/offset
    label blocki = pos / elem_per_block;
    unsigned int off = pos % elem_per_block;

    for
    (
        unsigned int blockval = (blocks_[blocki] >> off);
        blockval;
        blockval >>= 1u
    )
    {
        if (blockval & 1u)
        {
            return pos;
        }
        ++pos;
    }

    // Normal block-wise search. Starting at the next block

    const label nblocks = num_blocks(size());
    for (++blocki; blocki < nblocks; ++blocki)
    {
        label pos = (blocki * elem_per_block);

        for
        (
            unsigned int blockval = blocks_[blocki];
            blockval;
            blockval >>= 1u
        )
        {
            if (blockval & 1u)
            {
                return pos;
            }
            ++pos;
        }
    }

    return -1;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool Foam::bitSet::all() const
{
    return -1 == first_not_block();
}


inline bool Foam::bitSet::any() const
{
    if (size())
    {
        const label nblocks = num_blocks(size());

        for (label blocki=0; blocki < nblocks; ++blocki)
        {
            if (blocks_[blocki])
            {
                return true;
            }
        }
    }

    return false;
}


inline bool Foam::bitSet::none() const
{
    return !any();
}


inline bool Foam::bitSet::uniform() const
{
    return (size() == 1 || (size() > 1 && (test(0) ? all() : none())));
}


inline unsigned int Foam::bitSet::count(const bool on) const
{
    unsigned int total = 0;

    const label nblocks = num_blocks(size());

    for (label blocki = 0; blocki < nblocks; ++blocki)
    {
        total += BitOps::bit_count(blocks_[blocki]);
    }

    if (!on)
    {
        // Return the number of bits that are OFF.
        return (unsigned(size()) - total);
    }

    return total;
}


inline bool Foam::bitSet::test(const label pos) const
{
    return get(pos);
}


inline Foam::labelList Foam::bitSet::sortedToc() const
{
    return toc();
}


inline void Foam::bitSet::swap(bitSet& bitset)
{
    PackedList<1>::swap(bitset);
}


inline void Foam::bitSet::transfer(bitSet& bitset)
{
    PackedList<1>::transfer(bitset);
}


inline void Foam::bitSet::assign(const bool val)
{
    if (empty())
    {
        return;  // Trivial case
    }

    const label nblocks = num_blocks(size());

    if (val)
    {
        for (label blocki=0; blocki < nblocks; ++blocki)
        {
            blocks_[blocki] = (~0u);
        }
        clear_trailing_bits();
    }
    else
    {
        for (label blocki=0; blocki < nblocks; ++blocki)
        {
            blocks_[blocki] = (0u);
        }
    }
}


inline void Foam::bitSet::set(const bitSet& bitset)
{
    orEq(bitset, false);  // Non-strict: Lets the set size grow.
}


inline Foam::label Foam::bitSet::set(const labelUList& locations)
{
    return setMany(locations.begin(), locations.end());
}


template<class Addr>
inline Foam::label Foam::bitSet::set
(
    const IndirectListBase<label, Addr>& locations
)
{
    return setMany(locations.begin(), locations.end());
}


inline Foam::label Foam::bitSet::unset(const labelUList& locations)
{
    return unset(locations.begin(), locations.end());
}


template<class Addr>
inline Foam::label Foam::bitSet::unset
(
    const IndirectListBase<label, Addr>& locations
)
{
    return unset(locations.begin(), locations.end());
}


inline Foam::bitSet& Foam::bitSet::unset(const bitSet& other)
{
    return minusEq(other);
}


inline void Foam::bitSet::flip()
{
    if (size())
    {
        const label nblocks = num_blocks(size());

        for (label blocki=0; blocki < nblocks; ++blocki)
        {
            blocks_[blocki] = ~(blocks_[blocki]);
        }
        clear_trailing_bits();
    }
}


inline void Foam::bitSet::flip(const label i)
{
    if (i >= 0 && i < size())
    {
        reference(this, i).flip();
    }
}


inline Foam::bitSet& Foam::bitSet::bound(const label maxSize)
{
    if (maxSize < size())
    {
        resize(maxSize);
    }

    return *this;
}


inline Foam::bitSet& Foam::bitSet::bound(const bitSet& other)
{
    return bound(other.size());
}


inline Foam::bitSet& Foam::bitSet::extend(const label minSize)
{
    if (size() < minSize)
    {
        resize(minSize);
    }

    return *this;
}


inline Foam::bitSet& Foam::bitSet::extend(const bitSet& other)
{
    return extend(other.size());
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline bool Foam::bitSet::operator()(const label pos) const
{
    return test(pos);
}


inline unsigned int Foam::bitSet::operator[](const label i) const
{
    return get(i);
}


inline Foam::bitSet::reference Foam::bitSet::operator[](const label i)
{
    #ifdef FULLDEBUG
    checkIndex(i);
    #endif
    return reference(this, i);
}


inline Foam::bitSet& Foam::bitSet::operator=(const bool val)
{
    PackedList<1>::operator=(val);
    return *this;
}


inline Foam::bitSet& Foam::bitSet::operator=(const bitSet& bitset)
{
    PackedList<1>::operator=(bitset);
    return *this;
}


inline Foam::bitSet& Foam::bitSet::operator=(bitSet&& bitset)
{
    transfer(bitset);
    return *this;
}


inline Foam::bitSet& Foam::bitSet::operator&=(const bitSet& other)
{
    return andEq(other);
}


inline Foam::bitSet& Foam::bitSet::operator|=(const bitSet& other)
{
    return orEq(other);
}


inline Foam::bitSet& Foam::bitSet::operator^=(const bitSet& other)
{
    return xorEq(other);
}


inline Foam::bitSet& Foam::bitSet::operator-=(const bitSet& other)
{
    return minusEq(other);
}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

inline Foam::bitSet Foam::operator~(const bitSet& bitset)
{
    bitSet result(bitset);
    result.flip();
    return result;
}


inline Foam::bitSet Foam::operator&(const bitSet& a, const bitSet& b)
{
    bitSet result(a);
    return (result &= b);
}


inline Foam::bitSet Foam::operator|(const bitSet& a, const bitSet& b)
{
    bitSet result(a);
    return (result |= b);
}


inline Foam::bitSet Foam::operator^(const bitSet& a, const bitSet& b)
{
    bitSet result(a);
    return (result ^= b);
}


inline Foam::bitSet Foam::operator-(const bitSet& a, const bitSet& b)
{
    bitSet result(a);
    return (result -= b);
}


// ************************************************************************* //
