/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017-2018 OpenFOAM Foundation
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

// * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

inline Foam::label Foam::diameterModels::populationBalanceModel::nCorr() const
{
    return mesh_.solverDict(name_).get<label>("nCorr");
}


inline Foam::label
Foam::diameterModels::populationBalanceModel::sourceUpdateInterval() const
{
    return
        mesh_.solverDict(name_)
       .lookupOrDefault<label>("sourceUpdateInterval", 1);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::phaseSystem&
Foam::diameterModels::populationBalanceModel::fluid() const
{
    return fluid_;
}


inline const Foam::fvMesh&
Foam::diameterModels::populationBalanceModel::mesh() const
{
    return mesh_;
}


inline const Foam::dictionary&
Foam::diameterModels::populationBalanceModel::dict() const
{
    return dict_;
}


inline const Foam::phaseModel&
Foam::diameterModels::populationBalanceModel::continuousPhase() const
{
    return continuousPhase_;
}


inline const Foam::UPtrList<Foam::diameterModels::velocityGroup>&
Foam::diameterModels::populationBalanceModel::velocityGroups() const
{
    return velocityGroups_;
}


inline const Foam::UPtrList<Foam::diameterModels::sizeGroup>&
Foam::diameterModels::populationBalanceModel::sizeGroups() const
{
    return sizeGroups_;
}


inline const Foam::diameterModels::populationBalanceModel::phasePairTable&
Foam::diameterModels::populationBalanceModel::phasePairs() const
{
    return phasePairs_;
}


inline const Foam::PtrList<Foam::dimensionedScalar>&
Foam::diameterModels::populationBalanceModel::v() const
{
    return v_;
}


inline const Foam::volScalarField&
Foam::diameterModels::populationBalanceModel::alphas() const
{
    if (velocityGroups_.size() > 1)
    {
        return alphas_();
    }
    else
    {
        return velocityGroups_.first().phase();
    }
}


inline const Foam::volVectorField&
Foam::diameterModels::populationBalanceModel::U() const
{
    if (velocityGroups_.size() > 1)
    {
        return U_();
    }
    else
    {
        return velocityGroups_.first().phase().U();
    }
}


// ************************************************************************* //
