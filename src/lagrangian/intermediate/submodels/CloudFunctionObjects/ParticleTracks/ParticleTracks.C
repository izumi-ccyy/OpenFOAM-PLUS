/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2010-2011, 2019 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "ParticleTracks.H"
#include "Pstream.H"
#include "ListListOps.H"
#include "IOPtrList.H"

// * * * * * * * * * * * * * protected Member Functions  * * * * * * * * * * //

template<class CloudType>
void Foam::ParticleTracks<CloudType>::write()
{
    if (cloudPtr_.valid())
    {
        cloudPtr_->write();

        if (resetOnWrite_)
        {
            cloudPtr_->clear();
        }
    }
    else
    {
        DebugInFunction << "invalid cloud pointer" << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ParticleTracks<CloudType>::ParticleTracks
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName
)
:
    CloudFunctionObject<CloudType>(dict, owner, modelName, typeName),
    trackInterval_(this->coeffDict().getLabel("trackInterval")),
    maxSamples_(this->coeffDict().getLabel("maxSamples")),
    resetOnWrite_(this->coeffDict().getBool("resetOnWrite")),
    faceHitCounter_(),
    cloudPtr_(nullptr)
{}


template<class CloudType>
Foam::ParticleTracks<CloudType>::ParticleTracks
(
    const ParticleTracks<CloudType>& ppm
)
:
    CloudFunctionObject<CloudType>(ppm),
    trackInterval_(ppm.trackInterval_),
    maxSamples_(ppm.maxSamples_),
    resetOnWrite_(ppm.resetOnWrite_),
    faceHitCounter_(ppm.faceHitCounter_),
    cloudPtr_(ppm.cloudPtr_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ParticleTracks<CloudType>::~ParticleTracks()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::ParticleTracks<CloudType>::preEvolve()
{
    if (!cloudPtr_.valid())
    {
        cloudPtr_.reset
        (
            this->owner().cloneBare(this->owner().name() + "Tracks").ptr()
        );
    }
}


template<class CloudType>
void Foam::ParticleTracks<CloudType>::postFace(const parcelType& p, bool&)
{
    if
    (
        this->owner().solution().output()
     || this->owner().solution().transient()
    )
    {
        if (!cloudPtr_.valid())
        {
            FatalErrorInFunction
             << "Cloud storage not allocated" << abort(FatalError);
        }

        const label count =
            ++(faceHitCounter_(labelPair(p.origProc(), p.origId()), 0));

        const label nSamples = floor(count/trackInterval_);

        if ((count % trackInterval_) == 0 && nSamples < maxSamples_)
        {
            cloudPtr_->append
            (
                static_cast<parcelType*>(p.clone(this->owner().mesh()).ptr())
            );
        }
    }
}


// ************************************************************************* //
