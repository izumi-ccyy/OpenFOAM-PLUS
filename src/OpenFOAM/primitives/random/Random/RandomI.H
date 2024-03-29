/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenCFD Ltd.
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

inline Foam::scalar Foam::Random::scalar01()
{
    return uniform01_(generator_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline int Foam::Random::bit()
{
    return generator_() & 0x1;
}


inline Foam::label Foam::Random::seed() const
{
    return seed_;
}


inline void Foam::Random::reset(const label seedValue)
{
    seed_ = seedValue;
    generator_.seed(seed_);
    uniform01_.reset();   // A no-op, but for completeness
}


// ************************************************************************* //
