/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenCFD Ltd
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
    Foam::MHFModels::MHFModel

Description
    Base class for nucleation site density models

SourceFiles
    MHFModel.C
    newMHFModel.C

\*---------------------------------------------------------------------------*/

#ifndef MHFModel_H
#define MHFModel_H

#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"

#include "phaseModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace wallBoilingModels
{

/*---------------------------------------------------------------------------*\
                         Class MHFModel Declaration
\*---------------------------------------------------------------------------*/

class MHFModel
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        MHFModel(const MHFModel&);

        //- Disallow default bitwise assignment
        void operator=(const MHFModel&);


public:

    //- Runtime type information
    TypeName("MHFModel");


    //- Declare runtime construction
    declareRunTimeSelectionTable
    (
        autoPtr,
        MHFModel,
        dictionary,
        (
            const dictionary& dict
        ),
        (dict)
    );


    // Constructors

        //- Construct null
        MHFModel();


    // Selectors

        //- Select null constructed
        static autoPtr<MHFModel> New(const dictionary& dict);


    //- Destructor
    virtual ~MHFModel();


    // Member Functions

        //- Calculate temperature
        virtual tmp<scalarField> MHF
        (
            const phaseModel& liquid,
            const phaseModel& vapor,
            const label patchi,
            const scalarField& Tl,
            const scalarField& Tsatw,
            const scalarField& L
        ) const = 0;

        virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace wallBoilingModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
