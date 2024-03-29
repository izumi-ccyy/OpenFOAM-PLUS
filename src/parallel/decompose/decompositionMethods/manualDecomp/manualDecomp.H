/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2010, 2017-2018 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2011-2015 OpenFOAM Foundation
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
    Foam::manualDecomp

Description
    Decompose based on cell-to-processor association in a file

    Method coefficients:
    \table
        Property  | Description                           | Required | Default
        dataFile  | filename of cell decomposition data   | yes      |
    \endtable

SourceFiles
    manualDecomp.C

\*---------------------------------------------------------------------------*/

#ifndef manualDecomp_H
#define manualDecomp_H

#include "decompositionMethod.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class manualDecomp Declaration
\*---------------------------------------------------------------------------*/

class manualDecomp
:
    public decompositionMethod
{
    // Private data

        fileName dataFile_;


    // Private Member Functions

        //- No copy construct
        manualDecomp(const manualDecomp&) = delete;

        //- No copy assignment
        void operator=(const manualDecomp&) = delete;


public:

    //- Runtime type information
    TypeName("manual");


    // Constructors

        //- Construct given the decomposition dictionary
        manualDecomp(const dictionary& decompDict);

        //- Construct given the decomposition dictionary and region name
        manualDecomp
        (
            const dictionary& decompDict,
            const word& regionName
        );


    //- Destructor
    virtual ~manualDecomp() = default;


    // Member Functions

        //- Manual decompose does not care about proc boundaries.
        //- Is all up to the user.
        virtual bool parallelAware() const
        {
            return true;
        }

        //- Return for every coordinate the wanted processor number.
        //  Use the mesh connectivity (if needed)
        virtual labelList decompose
        (
            const polyMesh& mesh,
            const pointField& cc,
            const scalarField& cWeights
        ) const;

        //- Return for every coordinate the wanted processor number.
        //  Explicitly provided connectivity - does not use mesh_.
        virtual labelList decompose
        (
            const labelListList& globalCellCells,
            const pointField& cc,
            const scalarField& cWeights
        ) const
        {
            NotImplemented;
            return labelList();
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
