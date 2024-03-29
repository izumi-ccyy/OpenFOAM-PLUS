/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "Lee.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo, class OtherThermo>
Foam::meltingEvaporationModels::Lee<Thermo, OtherThermo>::Lee
(
    const dictionary& dict,
    const phasePair& pair
)
:
    InterfaceCompositionModel<Thermo, OtherThermo>(dict, pair),
    C_("C", inv(dimTime), dict),
    Tactivate_("Tactivate", dimTemperature, dict),
    alphaMin_(dict.lookupOrDefault<scalar>("alphaMin", 0))
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Thermo, class OtherThermo>
Foam::tmp<Foam::volScalarField>
Foam::meltingEvaporationModels::Lee<Thermo, OtherThermo>::Kexp
(
    label variable,
    const volScalarField& refValue
)
{
    if (this->modelVariable_ == variable)
    {
        volScalarField from
        (
            min(max(this->pair().from(), scalar(0)), scalar(1))
        );

        if (sign(C_.value()) > 0)
        {
            return
            (
                C_
              * from
              * this->pair().from().rho()
              * (refValue.oldTime() - Tactivate_)
              * pos(from - alphaMin_)
              * pos(refValue.oldTime() - Tactivate_)/Tactivate_
            );
        }
        else
        {
            return
            (
               -C_
              * from
              * this->pair().from().rho()
              * pos(from - alphaMin_)
              * (Tactivate_ - refValue.oldTime())
              * pos(Tactivate_ - refValue.oldTime())/Tactivate_
            );

        }
    }
    else
    {
        return tmp<volScalarField> ();
    }
}


template<class Thermo, class OtherThermo>
const Foam::dimensionedScalar&
Foam::meltingEvaporationModels::Lee<Thermo, OtherThermo>::Tactivate() const
{
    return Tactivate_;
}


// ************************************************************************* //
