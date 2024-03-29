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

\*---------------------------------------------------------------------------*/

#include "fv.H"
#include "multivariateSurfaceInterpolationScheme.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::multivariateSurfaceInterpolationScheme<Type>::
multivariateSurfaceInterpolationScheme
(
    const fvMesh& mesh,
    const multivariateSurfaceInterpolationScheme<Type>::fieldTable& vtfs,
    const surfaceScalarField&,
    Istream&
)
:
    mesh_(mesh),
    fields_(vtfs)
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::multivariateSurfaceInterpolationScheme<Type>>
Foam::multivariateSurfaceInterpolationScheme<Type>::New
(
    const fvMesh& mesh,
    const multivariateSurfaceInterpolationScheme<Type>::fieldTable& vtfs,
    const surfaceScalarField& faceFlux,
    Istream& schemeData
)
{
    if (fv::debug)
    {
        InfoInFunction
            << "Constructing surfaceInterpolationScheme<Type>" << endl;
    }

    const word schemeName(schemeData);

    auto cstrIter = IstreamConstructorTablePtr_->cfind(schemeName);

    if (!cstrIter.found())
    {
        FatalIOErrorInFunction(schemeData)
            << "Unknown discretisation scheme " << schemeName << nl << nl
            << "Valid schemes are :" << endl
            << IstreamConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    return cstrIter()(mesh, vtfs, faceFlux, schemeData);
}


// ************************************************************************* //
