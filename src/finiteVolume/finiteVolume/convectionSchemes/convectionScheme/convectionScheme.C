/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2010 OpenCFD Ltd.
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

Description
    Abstract base class for finite volume calculus convection schemes.

\*---------------------------------------------------------------------------*/

#include "fv.H"
#include "HashTable.H"
#include "linear.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
convectionScheme<Type>::convectionScheme(const convectionScheme& cs)
:
    refCount(),
    mesh_(cs.mesh_)
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

template<class Type>
tmp<convectionScheme<Type>> convectionScheme<Type>::New
(
    const fvMesh& mesh,
    const surfaceScalarField& faceFlux,
    Istream& schemeData
)
{
    if (fv::debug)
    {
        InfoInFunction << "Constructing convectionScheme<Type>" << endl;
    }

    if (schemeData.eof())
    {
        FatalIOErrorInFunction(schemeData)
            << "Convection scheme not specified" << endl << endl
            << "Valid convection schemes are :" << endl
            << IstreamConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    const word schemeName(schemeData);

    if (fv::debug)
    {
        InfoInFunction << "schemeName:" << schemeName << endl;
    }

    auto cstrIter = IstreamConstructorTablePtr_->cfind(schemeName);

    if (!cstrIter.found())
    {
        FatalIOErrorInFunction(schemeData)
            << "Unknown convection scheme "
            << schemeName << nl << nl
            << "Valid convection schemes are :" << endl
            << IstreamConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    return cstrIter()(mesh, faceFlux, schemeData);
}


template<class Type>
tmp<convectionScheme<Type>> convectionScheme<Type>::New
(
    const fvMesh& mesh,
    const typename multivariateSurfaceInterpolationScheme<Type>::
        fieldTable& fields,
    const surfaceScalarField& faceFlux,
    Istream& schemeData
)
{
    if (fv::debug)
    {
        InfoInFunction << "Constructing convectionScheme<Type>" << endl;
    }

    if (schemeData.eof())
    {
        FatalIOErrorInFunction(schemeData)
            << "Convection scheme not specified" << endl << endl
            << "Valid convection schemes are :" << endl
            << MultivariateConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    const word schemeName(schemeData);

    auto cstrIter = MultivariateConstructorTablePtr_->cfind(schemeName);

    if (!cstrIter.found())
    {
        FatalIOErrorInFunction(schemeData)
            << "Unknown convection scheme "
            << schemeName << nl << nl
            << "Valid convection schemes are :" << endl
            << MultivariateConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    return cstrIter()(mesh, fields, faceFlux, schemeData);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
convectionScheme<Type>::~convectionScheme()
{}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void convectionScheme<Type>::operator=(const convectionScheme<Type>& cs)
{
    if (this == &cs)
    {
        FatalErrorInFunction
            << "attempted assignment to self"
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
