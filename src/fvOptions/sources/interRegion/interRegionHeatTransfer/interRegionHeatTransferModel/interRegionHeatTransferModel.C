/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
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

#include "interRegionHeatTransferModel.H"
#include "basicThermo.H"
#include "fvmSup.H"
#include "zeroGradientFvPatchFields.H"
#include "fvcVolumeIntegrate.H"
#include "fvOptionList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(interRegionHeatTransferModel, 0);
}
}


// * * * * * * * * * * * *  Protected member functions * * * * * * * * * * * //

void Foam::fv::interRegionHeatTransferModel::setNbrModel()
{
    if (!firstIter_)
    {
        return;
    }

    const fvMesh& nbrMesh = mesh_.time().lookupObject<fvMesh>(nbrRegionName_);

    const optionList& fvOptions = nbrMesh.lookupObject<optionList>("fvOptions");

    bool nbrModelFound = false;

    forAll(fvOptions, i)
    {
        if (fvOptions[i].name() == nbrModelName_)
        {
            nbrModel_ = &const_cast<interRegionHeatTransferModel&>
            (
                refCast<const interRegionHeatTransferModel>(fvOptions[i])
            );
            nbrModelFound = true;
            break;
        }
    }

    if (!nbrModelFound)
    {
        FatalErrorInFunction
            << "Neighbour model not found" << nbrModelName_
            << " in region " << nbrMesh.name() << nl
            << exit(FatalError);
    }

    firstIter_ = false;

    // Set nbr model's nbr model to avoid construction order problems
    nbrModel_->setNbrModel();
}


void Foam::fv::interRegionHeatTransferModel::correct()
{
    if (master_)
    {
        if (mesh_.time().timeIndex() != timeIndex_)
        {
            calculateHtc();
            timeIndex_ = mesh_.time().timeIndex();
        }
    }
    else
    {
        nbrModel().correct();
        interpolate(nbrModel().htc(), htc_);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::interRegionHeatTransferModel::interRegionHeatTransferModel
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    interRegionOption
    (
        name,
        modelType,
        dict,
        mesh
    ),
    nbrModelName_(coeffs_.get<word>("nbrModel")),
    nbrModel_(nullptr),
    firstIter_(true),
    semiImplicit_(false),
    timeIndex_(-1),
    htc_
    (
        IOobject
        (
            type() + ":htc",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimEnergy/dimTime/dimTemperature/dimVolume, Zero),
        zeroGradientFvPatchScalarField::typeName
    ),
    TName_(coeffs_.lookupOrDefault<word>("T", "T")),
    TNbrName_(coeffs_.lookupOrDefault<word>("TNbr", "T"))
{
    if (active())
    {
        coeffs_.readEntry("fields", fieldNames_);
        applied_.setSize(fieldNames_.size(), false);

        coeffs_.readEntry("semiImplicit", semiImplicit_);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::interRegionHeatTransferModel::addSup
(
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    setNbrModel();

    correct();

    const volScalarField& he = eqn.psi();

    const volScalarField& T = mesh_.lookupObject<volScalarField>(TName_);

    auto tTmapped = tmp<volScalarField>::New
    (
        IOobject
        (
            type() + ":Tmapped",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        T
    );

    auto& Tmapped = tTmapped.ref();

    const fvMesh& nbrMesh = mesh_.time().lookupObject<fvMesh>(nbrRegionName_);

    const volScalarField& Tnbr =
        nbrMesh.lookupObject<volScalarField>(TNbrName_);

    interpolate(Tnbr, Tmapped.primitiveFieldRef());

    if (debug)
    {
        Info<< "Volumetric integral of htc: "
            << fvc::domainIntegrate(htc_).value()
            << endl;

        if (mesh_.time().writeTime())
        {
            Tmapped.write();
            htc_.write();
        }
    }

    if (semiImplicit_)
    {
        if (he.dimensions() == dimEnergy/dimMass)
        {
            const basicThermo* thermoPtr =
                mesh_.findObject<basicThermo>(basicThermo::dictName);

            if (thermoPtr)
            {
                const basicThermo& thermo = *thermoPtr;

                volScalarField htcByCpv(htc_/thermo.Cpv());

                eqn += htc_*(Tmapped - T) + htcByCpv*he - fvm::Sp(htcByCpv, he);

                if (debug)
                {
                    const dimensionedScalar energy =
                        fvc::domainIntegrate(htc_*(he/thermo.Cp() - Tmapped));

                    Info<< "Energy exchange from region " << nbrMesh.name()
                        << " To " << mesh_.name() << " : " <<  energy.value()
                        << endl;
                }
            }
            else
            {
                FatalErrorInFunction
                    << " on mesh " << mesh_.name()
                    << " could not find object basicThermo."
                    << " The available objects are: "
                    << mesh_.names()
                    << exit(FatalError);
            }
        }
        else if (he.dimensions() == dimTemperature)
        {
            eqn += htc_*Tmapped - fvm::Sp(htc_, he);
        }
    }
    else
    {
        eqn += htc_*(Tmapped - T);
    }
}


void Foam::fv::interRegionHeatTransferModel::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    addSup(eqn, fieldi);
}


bool Foam::fv::interRegionHeatTransferModel::read(const dictionary& dict)
{
    if (interRegionOption::read(dict))
    {
        return true;
    }

    return false;
}


// ************************************************************************* //
