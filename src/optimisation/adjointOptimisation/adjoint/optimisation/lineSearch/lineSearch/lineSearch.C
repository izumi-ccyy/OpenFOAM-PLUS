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

#include "lineSearch.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(lineSearch, 0);
defineRunTimeSelectionTable(lineSearch, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

lineSearch::lineSearch(const dictionary& dict, const Time& time)
:
    dict_(dict),
    lineSearchDict_
    (
        IOobject
        (
            "lineSearch",
            time.timeName(),
            "uniform",
            time,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        )
    ),
    directionalDeriv_(Zero),
    direction_(0),
    oldMeritValue_(Zero),
    newMeritValue_(Zero),
    prevMeritDeriv_
    (
        lineSearchDict_.lookupOrDefault<scalar>("prevMeritDeriv", Zero)
    ),
    initialStep_(dict.lookupOrDefault<scalar>("initialStep", 1)),
    minStep_(dict.lookupOrDefault<scalar>("minStep", 0.3)),
    step_(Zero),
    iter_(lineSearchDict_.lookupOrDefault<label>("iter", 0)),
    maxIters_(dict.lookupOrDefault<scalar>("maxIters", 10)),
    extrapolateInitialStep_
    (
        dict.lookupOrDefault<bool>
        (
            "extrapolateInitialStep",
            false
        )
    ),
    stepUpdate_(stepUpdate::New(dict))
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<lineSearch> lineSearch::New
(
    const dictionary& dict,
    const Time& time
)
{
    autoPtr<lineSearch> lineSrch(nullptr);

    const word modelType =
        dict.lookupOrDefault<word>("lineSearchType", "none");

    Info<< "lineSearch type : " << modelType << endl;

    if (modelType != "none")
    {
        auto cstrIter = dictionaryConstructorTablePtr_->cfind(modelType);

        if (!cstrIter.found())
        {
            FatalErrorInFunction
                << "Unknown lineSearch type " << modelType
                << nl << nl
                << "Valid lineSearch types are :" << nl
                << dictionaryConstructorTablePtr_->sortedToc()
                << exit(FatalError);
        }

        lineSrch.reset((cstrIter()(dict, time)).ptr());
    }
    else
    {
        Info<< "No line search method specified. "
            << "Proceeding with constant step" << endl;
    }

    return lineSrch;
}


// * * * * * * * * * * * * * * *  Member Functions   * * * * * * * * * * * * //

void lineSearch::setDeriv(const scalar deriv)
{
    directionalDeriv_ = deriv;
    stepUpdate_->setDeriv(deriv);
}


void lineSearch::setDirection(const scalarField& direction)
{
    direction_ = direction;
}


void lineSearch::setNewMeritValue(const scalar value)
{
    newMeritValue_ = value;
    stepUpdate_->setNewMeritValue(value);
}


void lineSearch::setOldMeritValue(const scalar value)
{
    oldMeritValue_ = value;
    stepUpdate_->setOldMeritValue(value);
}


void lineSearch::reset()
{
    if (extrapolateInitialStep_ && iter_ != 0)
    {
        // step_ = 2*(oldMeritValue_-prevMeritValue_)/directionalDeriv_;
        // Interpolate in order to get same improvement with the previous
        // optimisation cycle
        step_ = max(min(step_*prevMeritDeriv_/directionalDeriv_, 1.), minStep_);
        Info<< "\n------- Computing initial step-------" << endl;
        Info<< "old dphi(0) "  << prevMeritDeriv_ << endl;
        Info<< "dphi(0) "      << directionalDeriv_ << endl;
        Info<< "Setting initial step value " << step_ << endl << endl;
    }
    else
    {
        step_ = initialStep_;
    }
}


label lineSearch::maxIters() const
{
    return maxIters_;
}


scalar lineSearch::step() const
{
    return step_;
}


void lineSearch::updateStep(const scalar newStep)
{
    step_ = newStep;
}


lineSearch& lineSearch::operator++()
{
    iter_++;
    prevMeritDeriv_ = directionalDeriv_;
    lineSearchDict_.add<scalar>("prevMeritDeriv_", prevMeritDeriv_, true);
    lineSearchDict_.add<label>("iter", iter_, true);
    lineSearchDict_.regIOobject::write();

    return *this;
}


lineSearch& lineSearch::operator++(int)
{
    return operator++();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
