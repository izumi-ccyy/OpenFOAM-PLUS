/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2018 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2013-2016 OpenFOAM Foundation
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

#include "turbulenceFields.H"
#include "turbulentTransportModel.H"
#include "turbulentFluidThermoModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(turbulenceFields, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        turbulenceFields,
        dictionary
    );
}
}

const Foam::Enum
<
    Foam::functionObjects::turbulenceFields::compressibleField
>
Foam::functionObjects::turbulenceFields::compressibleFieldNames_
({
    { compressibleField::cfK, "k" },
    { compressibleField::cfEpsilon, "epsilon" },
    { compressibleField::cfOmega, "omega" },
    { compressibleField::cfNuTilda, "nuTilda" },
    { compressibleField::cfMut, "mut" },
    { compressibleField::cfMuEff, "muEff" },
    { compressibleField::cfAlphat, "alphat" },
    { compressibleField::cfAlphaEff, "alphaEff" },
    { compressibleField::cfR, "R" },
    { compressibleField::cfDevRhoReff, "devRhoReff" },
    { compressibleField::cfL, "L" },
    { compressibleField::cfI, "I" },
});


const Foam::Enum
<
    Foam::functionObjects::turbulenceFields::incompressibleField
>
Foam::functionObjects::turbulenceFields::incompressibleFieldNames_
({
    { incompressibleField::ifK, "k" },
    { incompressibleField::ifEpsilon, "epsilon" },
    { incompressibleField::ifOmega, "omega" },
    { incompressibleField::ifNuTilda, "nuTilda" },
    { incompressibleField::ifNut, "nut" },
    { incompressibleField::ifNuEff, "nuEff" },
    { incompressibleField::ifR, "R" },
    { incompressibleField::ifDevReff, "devReff" },
    { incompressibleField::ifL, "L" },
    { incompressibleField::ifI, "I" },
});


const Foam::word Foam::functionObjects::turbulenceFields::modelName
(
    Foam::turbulenceModel::propertiesName
);


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::turbulenceFields::compressible()
{
    if (obr_.foundObject<compressible::turbulenceModel>(modelName))
    {
        return true;
    }
    else if (obr_.foundObject<incompressible::turbulenceModel>(modelName))
    {
        return false;
    }

    FatalErrorInFunction
        << "Turbulence model not found in database, deactivating"
        << exit(FatalError);

    return false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::turbulenceFields::turbulenceFields
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    fieldSet_()
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::turbulenceFields::~turbulenceFields()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::turbulenceFields::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    if (dict.found("field"))
    {
        fieldSet_.insert(dict.get<word>("field"));
    }
    else
    {
        fieldSet_.insert(dict.get<wordList>("fields"));
    }

    Info<< type() << " " << name() << ": ";
    if (fieldSet_.size())
    {
        Info<< "storing fields:" << nl;
        for (const word& f : fieldSet_)
        {
            Info<< "    " << modelName << ':' << f << nl;
        }
        Info<< endl;
    }
    else
    {
        Info<< "no fields requested to be stored" << nl << endl;
    }

    return true;
}


bool Foam::functionObjects::turbulenceFields::execute()
{
    bool comp = compressible();

    if (comp)
    {
        const compressible::turbulenceModel& model =
            obr_.lookupObject<compressible::turbulenceModel>(modelName);

        for (const word& f : fieldSet_)
        {
            switch (compressibleFieldNames_[f])
            {
                case cfK:
                {
                    processField<scalar>(f, model.k());
                    break;
                }
                case cfEpsilon:
                {
                    processField<scalar>(f, model.epsilon());
                    break;
                }
                case cfOmega:
                {
                    processField<scalar>(f, omega(model));
                    break;
                }
                case cfNuTilda:
                {
                    processField<scalar>(f, nuTilda(model));
                    break;
                }
                case cfMut:
                {
                    processField<scalar>(f, model.mut());
                    break;
                }
                case cfMuEff:
                {
                    processField<scalar>(f, model.muEff());
                    break;
                }
                case cfAlphat:
                {
                    processField<scalar>(f, model.alphat());
                    break;
                }
                case cfAlphaEff:
                {
                    processField<scalar>(f, model.alphaEff());
                    break;
                }
                case cfR:
                {
                    processField<symmTensor>(f, model.R());
                    break;
                }
                case cfDevRhoReff:
                {
                    processField<symmTensor>(f, model.devRhoReff());
                    break;
                }
                case cfL:
                {
                    processField<scalar>(f, L(model));
                    break;
                }
                case cfI:
                {
                    processField<scalar>(f, I(model));
                    break;
                }
                default:
                {
                    FatalErrorInFunction
                        << "Invalid field selection" << abort(FatalError);
                }
            }
        }
    }
    else
    {
        const incompressible::turbulenceModel& model =
            obr_.lookupObject<incompressible::turbulenceModel>(modelName);

        for (const word& f : fieldSet_)
        {
            switch (incompressibleFieldNames_[f])
            {
                case ifK:
                {
                    processField<scalar>(f, model.k());
                    break;
                }
                case ifEpsilon:
                {
                    processField<scalar>(f, model.epsilon());
                    break;
                }
                case ifOmega:
                {
                    processField<scalar>(f, omega(model));
                    break;
                }
                case ifNuTilda:
                {
                    processField<scalar>(f, nuTilda(model));
                    break;
                }
                case ifNut:
                {
                    processField<scalar>(f, model.nut());
                    break;
                }
                case ifNuEff:
                {
                    processField<scalar>(f, model.nuEff());
                    break;
                }
                case ifR:
                {
                    processField<symmTensor>(f, model.R());
                    break;
                }
                case ifDevReff:
                {
                    processField<symmTensor>(f, model.devReff());
                    break;
                }
                case ifL:
                {
                    processField<scalar>(f, L(model));
                    break;
                }
                case ifI:
                {
                    processField<scalar>(f, I(model));
                    break;
                }
                default:
                {
                    FatalErrorInFunction
                        << "Invalid field selection" << abort(FatalError);
                }
            }
        }
    }

    return true;
}


bool Foam::functionObjects::turbulenceFields::write()
{
    for (const word& f : fieldSet_)
    {
        const word fieldName = modelName + ':' + f;
        writeObject(fieldName);
    }

    return true;
}


// ************************************************************************* //
