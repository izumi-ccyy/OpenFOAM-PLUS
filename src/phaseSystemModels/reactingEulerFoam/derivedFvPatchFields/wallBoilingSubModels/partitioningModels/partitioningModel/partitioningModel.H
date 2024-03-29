/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2016-2018 OpenFOAM Foundation
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
    Foam::wallBoilingModels::partitioningModel

Description
    Base class for wall heat flux partitioning models

SourceFiles
    partitioningModel.C
    newpartitioningModel.C

\*---------------------------------------------------------------------------*/

#ifndef partitioningModel_H
#define partitioningModel_H

#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace wallBoilingModels
{

/*---------------------------------------------------------------------------*\
                         Class partitioningModel Declaration
\*---------------------------------------------------------------------------*/

class partitioningModel
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        partitioningModel(const partitioningModel&);

        //- Disallow default bitwise assignment
        void operator=(const partitioningModel&);


public:

    //- Runtime type information
    TypeName("partitioningModel");


    //- Declare runtime construction
    declareRunTimeSelectionTable
    (
        autoPtr,
        partitioningModel,
        dictionary,
        (
            const dictionary& dict
        ),
        (dict)
    );


    // Constructors

        //- Construct null
        partitioningModel();


    // Selectors

        //- Select null constructed
        static autoPtr<partitioningModel> New(const dictionary& dict);


    //- Destructor
    virtual ~partitioningModel();


    // Member Functions

        //- Calculate and return the wall heat-flux partitioning
        virtual tmp<scalarField> fLiquid
        (
            const scalarField& alphaLiquid
        ) const = 0;

        virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace wallBoilingModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
