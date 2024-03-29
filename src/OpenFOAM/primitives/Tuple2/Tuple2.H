/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2010, 2019 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2011-2016 OpenFOAM Foundation
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

Class
    Foam::Tuple2

Description
    A 2-tuple for storing two objects of dissimilar types.
    The container is similar in purpose to std::pair, but does not expose
    its members directly.

See also
    Foam::Pair for storing two objects of identical types.

\*---------------------------------------------------------------------------*/

#ifndef Tuple2_H
#define Tuple2_H

#include "Istream.H"
#include "Ostream.H"
#include <utility>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class Tuple2 Declaration
\*---------------------------------------------------------------------------*/

template<class T1, class T2>
class Tuple2
{
    // Private Data

        T1 f_;
        T2 s_;

public:

    // Typedefs (cf. std::pair)

        //- Type of member first, the first template parameter (T1)
        typedef T1 first_type;

        //- Type of member second, the second template parameter (T2)
        typedef T2 second_type;


    // Constructors

        //- Construct null
        inline Tuple2() = default;

        //- Copy construct from components
        inline Tuple2(const T1& f, const T2& s)
        :
            f_(f),
            s_(s)
        {}

        //- Move construct from components
        inline Tuple2(T1&& f, T2&& s)
        :
            f_(std::move(f)),
            s_(std::move(s))
        {}

        //- Copy construct from std::pair
        inline Tuple2(const std::pair<T1,T2>& vals)
        :
            f_(vals.first),
            s_(vals.second)
        {}

        //- Move construct from std::pair
        inline Tuple2(std::pair<T1,T2>&& vals)
        :
            f_(std::move(vals.first)),
            s_(std::move(vals.second))
        {}

        //- Construct from Istream
        inline explicit Tuple2(Istream& is)
        {
            is >> *this;
        }


    // Member Functions

        //- Return first
        inline const T1& first() const noexcept
        {
            return f_;
        }

        //- Return first
        inline T1& first() noexcept
        {
            return f_;
        }

        //- Return second
        inline const T2& second() const noexcept
        {
            return s_;
        }

        //- Return second
        inline T2& second() noexcept
        {
            return s_;
        }
};


//- Return reverse of a Tuple2
template<class T1, class T2>
inline Tuple2<T2, T1> reverse(const Tuple2<T1,T2>& t)
{
    return Tuple2<T2, T1>(t.second(), t.first());
}


template<class T1, class T2>
inline bool operator==(const Tuple2<T1,T2>& a, const Tuple2<T1,T2>& b)
{
    return (a.first() == b.first() && a.second() == b.second());
}


template<class T1, class T2>
inline bool operator!=(const Tuple2<T1,T2>& a, const Tuple2<T1,T2>& b)
{
    return !(a == b);
}


template<class T1, class T2>
inline bool operator<(const Tuple2<T1,T2>& a, const Tuple2<T1,T2>& b)
{
    return
    (
        a.first() < b.first()
     || (!(b.first() < a.first()) && a.second() < b.second())
    );
}



template<class T1, class T2>
inline bool operator<=(const Tuple2<T1,T2>& a, const Tuple2<T1,T2>& b)
{
    return !(b < a);
}


template<class T1, class T2>
inline bool operator>(const Tuple2<T1,T2>& a, const Tuple2<T1,T2>& b)
{
    return (b < a);
}


template<class T1, class T2>
inline bool operator>=(const Tuple2<T1,T2>& a, const Tuple2<T1,T2>& b)
{
    return !(a < b);
}


// IOstream Operators

//- Read Tuple2 from Istream, discarding contents of existing Tuple2.
template<class T1, class T2>
inline Istream& operator>>(Istream& is, Tuple2<T1,T2>& t)
{
    is.readBegin("Tuple2");
    is >> t.first() >> t.second();
    is.readEnd("Tuple2");

    is.check(FUNCTION_NAME);
    return is;
}


//- Read std::pair from Istream
template<class T1, class T2>
inline Istream& operator>>(Istream& is, std::pair<T1,T2>& t)
{
    is.readBegin("std::pair");
    is >> t.first >> t.second;
    is.readEnd("std::pair");

    is.check(FUNCTION_NAME);
    return is;
}


//- Write Tuple2 to Ostream.
template<class T1, class T2>
inline Ostream& operator<<(Ostream& os, const Tuple2<T1,T2>& t)
{
    os  << token::BEGIN_LIST
        << t.first() << token::SPACE << t.second()
        << token::END_LIST;

    return os;
}


//- Write std::pair to Ostream.
template<class T1, class T2>
inline Ostream& operator<<(Ostream& os, const std::pair<T1,T2>& t)
{
    os  << token::BEGIN_LIST
        << t.first << token::SPACE << t.second
        << token::END_LIST;

    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
