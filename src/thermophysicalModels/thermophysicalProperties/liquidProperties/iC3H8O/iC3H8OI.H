/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2011 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2011-2017 OpenFOAM Foundation
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

inline Foam::scalar Foam::iC3H8O::rho(scalar p, scalar T) const
{
    return rho_.f(p, T);
}


inline Foam::scalar Foam::iC3H8O::pv(scalar p, scalar T) const
{
    return pv_.f(p, T);
}


inline Foam::scalar Foam::iC3H8O::hl(scalar p, scalar T) const
{
    return hl_.f(p, T);
}


inline Foam::scalar Foam::iC3H8O::Cp(scalar p, scalar T) const
{
    return Cp_.f(p, T);
}


inline Foam::scalar Foam::iC3H8O::h(scalar p, scalar T) const
{
    return h_.f(p, T);
}


inline Foam::scalar Foam::iC3H8O::Cpg(scalar p, scalar T) const
{
    return Cpg_.f(p, T);
}


inline Foam::scalar Foam::iC3H8O::B(scalar p, scalar T) const
{
    return B_.f(p, T);
}


inline Foam::scalar Foam::iC3H8O::mu(scalar p, scalar T) const
{
    return mu_.f(p, T);
}


inline Foam::scalar Foam::iC3H8O::mug(scalar p, scalar T) const
{
    return mug_.f(p, T);
}


inline Foam::scalar Foam::iC3H8O::kappa(scalar p, scalar T) const
{
    return kappa_.f(p, T);
}


inline Foam::scalar Foam::iC3H8O::kappag(scalar p, scalar T) const
{
    return kappag_.f(p, T);
}


inline Foam::scalar Foam::iC3H8O::sigma(scalar p, scalar T) const
{
    return sigma_.f(p, T);
}


inline Foam::scalar Foam::iC3H8O::D(scalar p, scalar T) const
{
    return D_.f(p, T);
}


inline Foam::scalar Foam::iC3H8O::D(scalar p, scalar T, scalar) const
{
    // Use 'standard' NSRDS D(p, T) function
    return D(p, T);
}


// ************************************************************************* //
