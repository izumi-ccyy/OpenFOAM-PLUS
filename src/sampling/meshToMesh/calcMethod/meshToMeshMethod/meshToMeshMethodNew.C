/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2013-2015 OpenFOAM Foundation
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

#include "meshToMeshMethod.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::meshToMeshMethod> Foam::meshToMeshMethod::New
(
    const word& methodName,
    const polyMesh& src,
    const polyMesh& tgt
)
{
    if (debug)
    {
        Info<< "Selecting AMIMethod " << methodName << endl;
    }

    auto cstrIter = componentsConstructorTablePtr_->cfind(methodName);

    if (!cstrIter.found())
    {
        FatalErrorInFunction
            << "Unknown meshToMesh type "
            << methodName << nl << nl
            << "Valid meshToMesh types :" << nl
            << componentsConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<meshToMeshMethod>(cstrIter()(src, tgt));
}


// ************************************************************************* //
