/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2010, 2019 OpenCFD Ltd.
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

Class
    Foam::nutUWallFunctionFvPatchScalarField

Group
    grpWallFunctions

Description
    This boundary condition provides a turbulent kinematic viscosity condition
    when using wall functions, based on velocity.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            nutUWallFunction;
    }
    \endverbatim

See also
    Foam::nutWallFunctionFvPatchScalarField

Note
    Suffers from non-exact restart since correctNut() (called through
    turbulence->validate) returns a slightly different value every time
    it is called. See nutUSpaldingWallFunctionFvPatchScalarField.C

SourceFiles
    nutUWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef nutUWallFunctionFvPatchScalarField_H
#define nutUWallFunctionFvPatchScalarField_H

#include "nutWallFunctionFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class nutUWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class nutUWallFunctionFvPatchScalarField
:
    public nutWallFunctionFvPatchScalarField
{
protected:

    // Protected Member Functions

        //- Calculate yPlus
        virtual tmp<scalarField> calcYPlus(const scalarField& magUp) const;

        //- Calculate the turbulence viscosity
        virtual tmp<scalarField> calcNut() const;


public:

    //- Runtime type information
    TypeName("nutUWallFunction");


    // Constructors

        //- Construct from patch and internal field
        nutUWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        nutUWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  nutUWallFunctionFvPatchScalarField
        //  onto a new patch
        nutUWallFunctionFvPatchScalarField
        (
            const nutUWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        nutUWallFunctionFvPatchScalarField
        (
            const nutUWallFunctionFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new nutUWallFunctionFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        nutUWallFunctionFvPatchScalarField
        (
            const nutUWallFunctionFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new nutUWallFunctionFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

            //- Calculate and return the yPlus at the boundary
            virtual tmp<scalarField> yPlus() const;


        // I-O

            //- Write
            virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
