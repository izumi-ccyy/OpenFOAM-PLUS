/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2015-2016 OpenFOAM Foundation
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

#include "temperatureDependentAlphaContactAngleFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::temperatureDependentAlphaContactAngleFvPatchScalarField::
temperatureDependentAlphaContactAngleFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    alphaContactAngleTwoPhaseFvPatchScalarField(p, iF),
    TName_("T"),
    theta0_()
{}


Foam::temperatureDependentAlphaContactAngleFvPatchScalarField::
temperatureDependentAlphaContactAngleFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    alphaContactAngleTwoPhaseFvPatchScalarField(p, iF, dict),
    TName_(dict.lookupOrDefault<word>("T", "T")),
    theta0_(Function1<scalar>::New("theta0", dict))
{
    evaluate();
}


Foam::temperatureDependentAlphaContactAngleFvPatchScalarField::
temperatureDependentAlphaContactAngleFvPatchScalarField
(
    const temperatureDependentAlphaContactAngleFvPatchScalarField& psf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    alphaContactAngleTwoPhaseFvPatchScalarField(psf, p, iF, mapper),
    TName_(psf.TName_),
    theta0_(psf.theta0_.clone())
{}


Foam::temperatureDependentAlphaContactAngleFvPatchScalarField::
temperatureDependentAlphaContactAngleFvPatchScalarField
(
    const temperatureDependentAlphaContactAngleFvPatchScalarField& psf
)
:
    alphaContactAngleTwoPhaseFvPatchScalarField(psf),
    TName_(psf.TName_),
    theta0_(psf.theta0_.clone())
{}


Foam::temperatureDependentAlphaContactAngleFvPatchScalarField::
temperatureDependentAlphaContactAngleFvPatchScalarField
(
    const temperatureDependentAlphaContactAngleFvPatchScalarField& psf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    alphaContactAngleTwoPhaseFvPatchScalarField(psf, iF),
    TName_(psf.TName_),
    theta0_(psf.theta0_.clone())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::temperatureDependentAlphaContactAngleFvPatchScalarField::theta
(
    const fvPatchVectorField&,
    const fvsPatchVectorField&
) const
{
    return theta0_->value
    (
        patch().lookupPatchField<volScalarField, scalar>(TName_)
    );
}


void Foam::temperatureDependentAlphaContactAngleFvPatchScalarField::write
(
    Ostream& os
) const
{
    alphaContactAngleTwoPhaseFvPatchScalarField::write(os);
    os.writeEntryIfDifferent<word>("T", "T", TName_);
    theta0_->writeData(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        temperatureDependentAlphaContactAngleFvPatchScalarField
    );
}

// ************************************************************************* //
