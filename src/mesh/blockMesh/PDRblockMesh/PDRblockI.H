/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::PDRblock::PDRblock()
:
    ijkMesh(),
    verbose_(false),
    grid_(),
    bounds_(),
    patches_(),
    minEdgeLen_(Zero)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool Foam::PDRblock::location::valid() const
{
    return (scalarList::size() > 1);
}


inline Foam::label Foam::PDRblock::location::nCells() const
{
    return (scalarList::size()-1);
}


inline Foam::label Foam::PDRblock::location::nPoints() const
{
    return scalarList::size();
}


inline bool Foam::PDRblock::location::contains(const scalar p) const
{
    return (scalarList::size() > 1 && first() <= p && p <= last());
}


inline const Foam::scalar& Foam::PDRblock::location::min() const
{
    return scalarList::empty() ? pTraits<scalar>::rootMax : first();
}


inline const Foam::scalar& Foam::PDRblock::location::max() const
{
    return scalarList::empty() ? pTraits<scalar>::rootMin : last();
}


inline Foam::scalar Foam::PDRblock::location::centre() const
{
    return scalarList::empty() ? 0 : (0.5*first() + 0.5*last());
}


inline void Foam::PDRblock::location::checkIndex(const label i) const
{
    if (i < 0 || i >= nCells())
    {
        FatalErrorInFunction
            << "The index " << i
            << " is out of range [0," << nCells() << ']' << nl
            << abort(FatalError);
    }
}


inline Foam::scalar Foam::PDRblock::location::width(const label i) const
{
    #ifdef FULLDEBUG
    checkIndex(i);
    #endif

    return (operator[](i+1) - operator[](i));
}


inline Foam::scalar Foam::PDRblock::location::C(const label i) const
{
    if (i == -1)
    {
        #ifdef FULLDEBUG
        checkIndex(0);
        #endif

        // "Halo" centre [-1] == x0 - 1/2 (x1 - x0)
        return first() - 0.5*(operator[](1) - first());
    }
    else if (i > 1 && i == scalarList::size()-1)
    {
        // "Halo" centre [nCells] == xN + 1/2 (xN - xN1)
        return last() + 0.5*(last() - operator[](scalarList::size()-2));
    }

    #ifdef FULLDEBUG
    checkIndex(i);
    #endif

    return 0.5*(operator[](i+1) + operator[](i));
}


inline const Foam::scalar&
Foam::PDRblock::location::clip(const scalar& val) const
{
    if (scalarList::size())
    {
        if (val < first())
        {
            return first();
        }
        else if (last() < val)
        {
            return last();
        }
    }

    return val; // Pass-through
}


inline const Foam::Vector<Foam::PDRblock::location>&
Foam::PDRblock::grid() const
{
    return grid_;
}


inline const Foam::scalar& Foam::PDRblock::minEdgeLen() const
{
    return minEdgeLen_;
}


inline const Foam::boundBox& Foam::PDRblock::bounds() const
{
    return bounds_;
}


inline Foam::scalar Foam::PDRblock::dx(const label i) const
{
    return grid_.x().width(i);
}


inline Foam::scalar Foam::PDRblock::dx(const labelVector& ijk) const
{
    return grid_.x().width(ijk.x());
}


inline Foam::scalar Foam::PDRblock::dy(const label j) const
{
    return grid_.y().width(j);
}


inline Foam::scalar Foam::PDRblock::dy(const labelVector& ijk) const
{
    return grid_.y().width(ijk.y());
}


inline Foam::scalar Foam::PDRblock::dz(const label k) const
{
    return grid_.z().width(k);
}


inline Foam::scalar Foam::PDRblock::dz(const labelVector& ijk) const
{
    return grid_.z().width(ijk.z());
}


inline Foam::vector Foam::PDRblock::span
(
    const label i,
    const label j,
    const label k
) const
{
    return vector(dx(i), dy(j), dz(k));
}


inline Foam::vector Foam::PDRblock::span(const labelVector& ijk) const
{
    return vector(dx(ijk), dy(ijk), dz(ijk));
}


inline Foam::point Foam::PDRblock::grid
(
    const label i,
    const label j,
    const label k
) const
{
    return point(grid_.x()[i], grid_.y()[j], grid_.z()[k]);
}


inline Foam::point Foam::PDRblock::grid(const labelVector& ijk) const
{
    return
        point
        (
            grid_.x()[ijk.x()],
            grid_.y()[ijk.y()],
            grid_.z()[ijk.z()]
        );
}


inline Foam::point Foam::PDRblock::C
(
    const label i,
    const label j,
    const label k
) const
{
    return point(grid_.x().C(i), grid_.y().C(j), grid_.z().C(k));
}


inline Foam::point Foam::PDRblock::C(const labelVector& ijk) const
{
    return
        point
        (
            grid_.x().C(ijk.x()),
            grid_.y().C(ijk.y()),
            grid_.z().C(ijk.z())
        );
}


inline Foam::scalar Foam::PDRblock::V
(
    const label i,
    const label j,
    const label k
) const
{
    return dx(i)*dy(j)*dz(k);
}


inline Foam::scalar Foam::PDRblock::V(const labelVector& ijk) const
{
    return dx(ijk.x())*dy(ijk.y())*dz(ijk.z());
}


inline Foam::scalar Foam::PDRblock::width
(
    const label i,
    const label j,
    const label k
) const
{
    return Foam::cbrt(V(i, j, k));
}


inline Foam::scalar Foam::PDRblock::width(const labelVector& ijk) const
{
    return Foam::cbrt(V(ijk));
}


// ************************************************************************* //
