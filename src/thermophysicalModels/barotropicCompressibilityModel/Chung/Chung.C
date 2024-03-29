/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2010 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "Chung.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace compressibilityModels
    {
        defineTypeNameAndDebug(Chung, 0);
        addToRunTimeSelectionTable
        (
            barotropicCompressibilityModel,
            Chung,
            dictionary
        );
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::compressibilityModels::Chung::Chung
(
    const dictionary& compressibilityProperties,
    const volScalarField& gamma,
    const word& psiName
)
:
    barotropicCompressibilityModel(compressibilityProperties, gamma, psiName),
    psiv_
    (
        "psiv",
        dimCompressibility,
        compressibilityProperties_
    ),
    psil_
    (
        "psil",
        dimCompressibility,
        compressibilityProperties_
    ),
    rhovSat_
    (
        "rhovSat",
        dimDensity,
        compressibilityProperties_
    ),
    rholSat_
    (
        "rholSat",
        dimDensity,
        compressibilityProperties_
    )
{
    correct();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::compressibilityModels::Chung::correct()
{
    volScalarField sfa
    (
        sqrt
        (
            (rhovSat_/psiv_)
           /((scalar(1) - gamma_)*rhovSat_/psiv_ + gamma_*rholSat_/psil_)
        )
    );

    psi_ = sqr
    (
        ((scalar(1) - gamma_)/sqrt(psiv_) + gamma_*sfa/sqrt(psil_))
       *sqrt(psiv_*psil_)/sfa
    );
}


bool Foam::compressibilityModels::Chung::read
(
    const dictionary& compressibilityProperties
)
{
    barotropicCompressibilityModel::read(compressibilityProperties);

    compressibilityProperties_.readEntry("psiv", psiv_);
    compressibilityProperties_.readEntry("psil", psil_);
    compressibilityProperties_.readEntry("rhovSat", rhovSat_);
    compressibilityProperties_.readEntry("rholSat", rholSat_);

    return true;
}


// ************************************************************************* //
