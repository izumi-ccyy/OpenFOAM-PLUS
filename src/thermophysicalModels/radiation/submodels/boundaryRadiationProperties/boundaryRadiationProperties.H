/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2018 OpenCFD Ltd.
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
    Foam::radiation::boundaryRadiationProperties

Description
    Boundary radiation properties holder

SourceFiles
    boundaryRadiationProperties.C

\*---------------------------------------------------------------------------*/

#ifndef boundaryRadiationProperties_H
#define boundaryRadiationProperties_H

#include "MeshObject.H"
#include "boundaryRadiationPropertiesPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;

namespace radiation
{

/*---------------------------------------------------------------------------*\
                 Class boundaryRadiationProperties Declaration
\*---------------------------------------------------------------------------*/

class boundaryRadiationProperties
:
    public MeshObject
    <
        fvMesh,
        Foam::GeometricMeshObject,
        boundaryRadiationProperties
    >
{
    // Private data

        //- Ptr list of boundaryRadiationProperties
        List<autoPtr<boundaryRadiationPropertiesPatch>>
            radBoundaryPropertiesPtrList_;


public:

    // Declare name of the class and its debug switch
    TypeName("boundaryRadiationProperties");


    // Constructors

        //- Construct given fvMesh
        explicit boundaryRadiationProperties(const fvMesh&);


    // Member Functions

        //- Access boundary emissivity on patch
        tmp<scalarField> emissivity
        (
            const label patchI,
            const label bandI = 0,
            vectorField* incomingDirection = nullptr,
            scalarField* T = nullptr
        ) const;

        //- Access boundary emissivity on face
        scalar faceEmissivity
        (
            const label patchI,
            const label faceI,
            const label bandI = 0,
            vector incomingDirection = Zero,
            scalar T = 0
        ) const;

        //- Access boundary absorptivity on patch
        tmp<scalarField> absorptivity
        (
            const label patchI,
            const label bandI = 0,
            vectorField* incomingDirection = nullptr,
            scalarField* T = nullptr
        ) const;

        //- Access boundary absorptivity on face
        scalar faceAbsorptivity
        (
            const label patchI,
            const label faceI,
            const label bandI = 0,
            vector incomingDirection = Zero,
            scalar T = 0
        ) const;

        //- Access boundary transmissivity on patch
        tmp<scalarField> transmissivity
        (
            const label patchI,
            const label bandI = 0,
            vectorField* incomingDirection = nullptr,
            scalarField* T = nullptr
        ) const;

        //- Access boundary transmissivity on face
        scalar faceTransmissivity
        (
            const label patchI,
            const label faceI,
            const label bandI = 0,
            vector incomingDirection = Zero,
            scalar T = 0
        ) const;

        //- Access boundary diffuse reflectivity on patch
        tmp<scalarField> diffReflectivity
        (
            const label patchI,
            const label bandI = 0,
            vectorField* incomingDirection = nullptr,
            scalarField* T = nullptr
        ) const;

        //- Access boundary diffuse reflectivity on face
        scalar faceDiffReflectivity
        (
            const label patchI,
            const label faceI,
            const label bandI = 0,
            vector incomingDirection = Zero,
            scalar T = 0
        ) const;

        //- Access boundary specular reflectivity on patch
        tmp<scalarField> specReflectivity
        (
            const label patchI,
            const label bandI = 0,
            vectorField* incomingDirection = nullptr,
            scalarField* T = nullptr
        ) const;

        //- Access boundary specular reflectivity on face
        scalar faceSpecReflectivity
        (
            const label patchI,
            const label faceI,
            const label bandI = 0,
            vector incomingDirection = Zero,
            scalar T = 0
        ) const;


    //- Destructor
    ~boundaryRadiationProperties() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
