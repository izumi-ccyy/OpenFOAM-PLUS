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

\*---------------------------------------------------------------------------*/

#include "adjointOutletVelocityFluxFvPatchVectorField.H"
#include "emptyFvPatch.H"
#include "fvMatrix.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::adjointOutletVelocityFluxFvPatchVectorField::
adjointOutletVelocityFluxFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    adjointBoundaryCondition(p, iF, word::null)
{}


Foam::adjointOutletVelocityFluxFvPatchVectorField::
adjointOutletVelocityFluxFvPatchVectorField
(
    const adjointOutletVelocityFluxFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    adjointBoundaryCondition(p, iF, ptf.adjointSolverName_)
{}


Foam::adjointOutletVelocityFluxFvPatchVectorField::
adjointOutletVelocityFluxFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    adjointBoundaryCondition(p, iF, dict.get<word>("solverName"))
{
    fvPatchField<vector>::operator=
    (
         vectorField("value", dict, p.size())
    );
}


Foam::adjointOutletVelocityFluxFvPatchVectorField::
adjointOutletVelocityFluxFvPatchVectorField
(
    const adjointOutletVelocityFluxFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pivpvf, iF),
    adjointBoundaryCondition(pivpvf)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::adjointOutletVelocityFluxFvPatchVectorField::manipulateMatrix
(
    fvMatrix<vector>& matrix
)
{
    vectorField& source = matrix.source();
    const vectorField& Sf = patch().Sf();
    const labelList& faceCells = patch().faceCells();
    const scalarField& magSf = patch().magSf();
    tmp<vectorField> tvelocitySource(boundaryContrPtr_->velocitySource());
    const vectorField& velocitySource = tvelocitySource();
    const fvPatchScalarField& pab = boundaryContrPtr_->pab();
    const word& fieldName = internalField().name();
    tmp<tensorField> tgradUab(computePatchGrad<vector>(fieldName));
    const tensorField& gradUab = tgradUab();

    // Momentum diffusion coefficient
    tmp<scalarField> tmomentumDiffusion(boundaryContrPtr_->momentumDiffusion());
    const scalarField& momentumDiffusion = tmomentumDiffusion();

    vectorField explDiffusiveFlux
    (
        -momentumDiffusion*(gradUab - sphericalTensor::oneThirdI*tr(gradUab))
       & Sf
    );

//  const fvPatchVectorField& Ub = boundaryContrPtr_->Ub();
//  const fvPatchVectorField& Uab = boundaryContrPtr_->Uab();
//  vectorField cmFormTerm = (Ub & Uab)*Sf;

    forAll(faceCells, fI)
    {
        const label cI = faceCells[fI];
        // Contributions from the convection and diffusion term (except from
        // the transpose part) will be canceled out through the value and
        // gradient coeffs. The pressure flux will be inserted later through
        // grad(pa), so it must be canceled out here. Once the typical fluxes
        // have been canceled out, add the objective flux. velocitySource
        // includes also fluxes from the adjoint turbulence-dependent terms
        // found in the adjoint momentum equations.
        source[cI] +=
            pab[fI]*Sf[fI]
//        - cmFormTerm[fI]
          + explDiffusiveFlux[fI]
          - velocitySource[fI]*magSf[fI];
    }
}


void Foam::adjointOutletVelocityFluxFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    tmp<vectorField> tnf = patch().nf();
    const vectorField& nf = tnf();

    // vectorField Ua = (patchInternalField() & nf) * nf;
    const fvsPatchScalarField& phia = boundaryContrPtr_->phiab();
    vectorField Ua((phia/patch().magSf())*nf);

    operator==(Ua);

    fixedValueFvPatchVectorField::updateCoeffs();
}


Foam::tmp<Foam::Field<Foam::vector>>
Foam::adjointOutletVelocityFluxFvPatchVectorField::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    return tmp<Field<vector>>::New(this->size(), pTraits<vector>::zero);
}


Foam::tmp<Foam::Field<Foam::vector>>
Foam::adjointOutletVelocityFluxFvPatchVectorField::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    return tmp<Field<vector>>::New(this->size(), pTraits<vector>::zero);
}


Foam::tmp<Foam::Field<Foam::vector>>
Foam::adjointOutletVelocityFluxFvPatchVectorField::
gradientBoundaryCoeffs() const
{
    return tmp<Field<vector>>::New(this->size(), pTraits<vector>::zero);
}


Foam::tmp<Foam::Field<Foam::vector>>
Foam::adjointOutletVelocityFluxFvPatchVectorField::
gradientInternalCoeffs() const
{
    return tmp<Field<vector>>::New(this->size(), pTraits<vector>::zero);
}


void Foam::adjointOutletVelocityFluxFvPatchVectorField::write
(
    Ostream& os
) const
{
    fvPatchVectorField::write(os);
    writeEntry("value", os);
    os.writeEntry("solverName", adjointSolverName_);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::adjointOutletVelocityFluxFvPatchVectorField::operator=
(
    const fvPatchField<vector>& pvf
)
{
    fvPatchField<vector>::operator=(patch().nf()*(patch().nf() & pvf));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        adjointOutletVelocityFluxFvPatchVectorField
    );
}

// ************************************************************************* //
