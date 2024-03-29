/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2007-2010, 2018-2019 OpenCFD Ltd.
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

Class
    Foam::scalarRange

Description
    Scalar bounds to be used as a unary predicate.

    The bound can be specified as an "MIN:MAX" range, as a "MIN:" or ":MAX"
    bound or simply as a single "VALUE".

    When defined via the parse() method, the special string "none" can be
    used to define an empty (inverse) range.

SeeAlso
    Foam::MinMax
    Foam::predicates::scalars

SourceFiles
    scalarRange.C

\*---------------------------------------------------------------------------*/
#ifndef scalarRange_H
#define scalarRange_H

#include "scalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class scalarRange;
class Ostream;
template<class T> class MinMax;

Ostream& operator<<(Ostream& os, const scalarRange& range);


/*---------------------------------------------------------------------------*\
                         Class scalarRange Declaration
\*---------------------------------------------------------------------------*/

class scalarRange
{
    //- Enumeration defining type of range test to use
    enum testType
    {
        NONE = 0,   //!< Empty or invalid (inverse) range. Never matches.
        EQ,         //!< A single value test.
        GE,         //!< Greater-than-equal.
        GT,         //!< Greater-than.
        LE,         //!< Less-than-equal.
        LT,         //!< Less-than.
        GE_LE,      //!< An inclusive range.
        ALWAYS      //!< Placeholder range that always matches
    };


    // Private Member Data

        //- The min value of the range
        scalar min_;

        //- The max value of the range
        scalar max_;

        //- The type of range test
        enum testType type_;

        //- Construct from components without sanity checking
        inline scalarRange
        (
            const testType type,
            const scalar minVal,
            const scalar maxVal
        ) noexcept;


public:

    // STL type definitions

        //- Type of values the range contains
        typedef scalar value_type;


    // Static Data Members

        //- An empty/inverted range that never matches
        static const scalarRange null;

        //- A range that always matches
        static const scalarRange always;


    // Constructors

        //- Construct an empty (inverse) range
        inline scalarRange() noexcept;

        //- Construct an exact value matcher
        inline explicit scalarRange(const scalar value) noexcept;

        //- Construct a range from min-value to max-value
        inline scalarRange(const scalar minVal, const scalar maxVal) noexcept;

        //- Copy construct from a min/max range.
        //  Automatically decides if this is a GE_LE or NONE range
        explicit scalarRange(const MinMax<label>& range) noexcept;

        //- Copy construct from a min/max range.
        //  Automatically decides if this is a GE_LE or NONE range
        explicit scalarRange(const MinMax<scalar>& range) noexcept;


    // Static Constructors

        //- Construct by parsing string content.
        //  A colon (:) is used as a range marker or when specifying
        //  greater-than or less-than bounds.
        //
        //  \note The special string "none" can be used define an empty
        //      (inverse) range
        //
        //  \return True if no parse problems were encountered.
        static bool parse(const std::string& str, scalarRange& range);

        //- Construct by parsing string content.
        //  Any parse problems are emitted as information and the returned
        //  range is of type empty().
        //  \return The parsed range, which is empty() on any problems
        static scalarRange parse(const std::string& str);


        //- Construct a greater-equals bound
        inline static scalarRange ge(const scalar minVal) noexcept;

        //- Construct a greater-than bound
        inline static scalarRange gt(const scalar minVal) noexcept;

        //- Construct a greater-equals zero bound
        inline static scalarRange ge0() noexcept;

        //- Construct a greater-than zero bound
        inline static scalarRange gt0() noexcept;

        //- Construct a less-equals bound
        inline static scalarRange le(const scalar maxVal) noexcept;

        //- Construct a less-than bound
        inline static scalarRange lt(const scalar maxVal) noexcept;

        //- Construct a greater-equals 0, less-equals 1 bound
        inline static scalarRange zero_one() noexcept;


    // Member Functions

        //- Reset to an empty (inverse) range.
        inline void clear() noexcept;

        //- True if range is empty (eg, inverted)
        inline bool empty() const noexcept;

        //- True if range is non-empty.
        inline bool valid() const noexcept;

        //- True if the range bounds represent a single value.
        inline bool single() const noexcept;

        //- The min value of the range.
        inline scalar min() const noexcept;

        //- The max value of the range.
        inline scalar max() const noexcept;

        //- A representative (average) value for the range.
        //  For GE, LE bounds it is the min/max value, respectively.
        inline scalar value() const;

        //- True if the value matches the condition.
        inline bool match(const scalar& val) const;


    // Member Operators

        //- Identical to match(), for use as a predicate.
        inline bool operator()(const scalar& val) const;

        inline bool operator==(const scalarRange& rhs) const noexcept;
        inline bool operator!=(const scalarRange& rhs) const noexcept;


    // IOstream Operators

        //- Print information about the range.
        friend Ostream& operator<<(Ostream& os, const scalarRange& range);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "scalarRangeI.H"

#endif

// ************************************************************************* //
