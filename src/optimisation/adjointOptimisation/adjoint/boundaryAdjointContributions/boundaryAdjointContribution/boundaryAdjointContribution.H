/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2007-2019 PCOpt/NTUA
                            | Copyright (C) 2013-2019 FOSS GP
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
    Foam::boundaryAdjointContribution

Description
    Abstract base class for computing contributions of the objective functions
    to the adjoint boundary conditions

SourceFiles
    boundaryAdjointContribution.C

\*---------------------------------------------------------------------------*/

#ifndef boundaryAdjointContribution_H
#define boundaryAdjointContribution_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "fvPatchFields.H"
#include "fvsPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class boundaryAdjointContribution Declaration
\*---------------------------------------------------------------------------*/

class boundaryAdjointContribution
{
private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        boundaryAdjointContribution
        (
            const boundaryAdjointContribution&
        ) = delete;

        //- Disallow default bitwise assignment
        void operator=(const boundaryAdjointContribution&) = delete;


protected:

    // Protected data

        const fvPatch& patch_;


public:

    //- Runtime type information
    TypeName("boundaryAdjointContribution");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            boundaryAdjointContribution,
            dictionary,
            (
                const word& managerName,
                const word& adjointSolverName,
                const word& simulationType,
                const fvPatch& patch
            ),
            (managerName, adjointSolverName, simulationType, patch)
        );

    // Constructors

        //- Construct from components
        boundaryAdjointContribution
        (
            const word& managerName,
            const word& adjointSolverName,
            const word& simulationType,
            const fvPatch& patch
        );


    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<boundaryAdjointContribution> New
        (
            const word& managerName,
            const word& adjointSolverName,
            const word& simulationType,
            const fvPatch& patch
        );


    //- Destructor
    virtual ~boundaryAdjointContribution() = default;


    // Member Functions

        // Contribution to surface sensitivities for a specific patch
        virtual tmp<scalarField> pressureSource() = 0;
        virtual tmp<vectorField> velocitySource() = 0;
        virtual tmp<vectorField> tangentVelocitySource() = 0;
        virtual tmp<vectorField> normalVelocitySource() = 0;
        virtual tmp<scalarField> adjointTMVariable1Source();
        virtual tmp<scalarField> adjointTMVariable2Source();
        virtual tmp<scalarField> energySource() = 0;

        virtual tmp<scalarField> momentumDiffusion() = 0;
        virtual tmp<scalarField> laminarDiffusivity() = 0;
        virtual tmp<scalarField> thermalDiffusion() = 0;
        virtual tmp<scalarField> wallDistance() = 0;

        virtual tmp<scalarField> TMVariable1Diffusion();
        virtual tmp<scalarField> TMVariable2Diffusion();
        virtual tmp<scalarField> TMVariable1();
        virtual tmp<scalarField> TMVariable2();

        // References to primal and adjoint fields for the specific patch
        virtual const fvPatchVectorField& Ub() const = 0;
        virtual const fvPatchScalarField& pb() const = 0;
        virtual const fvsPatchScalarField& phib() const = 0;
        virtual const fvPatchScalarField& turbulentDiffusivity() const = 0;
        virtual const fvPatchVectorField& Uab() const = 0;
        virtual const fvPatchScalarField& pab() const = 0;
        virtual const fvsPatchScalarField& phiab() const = 0;

        // Field suffixes for primal and adjoint fields
        virtual const word primalSolverName() const = 0;
        virtual const word adjointSolverName() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
