/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenCFD Ltd.
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

#include "faceAreaWeightAMI.H"
#include "profiling.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class SourcePatch, class TargetPatch>
void Foam::faceAreaWeightAMI<SourcePatch, TargetPatch>::calcAddressing
(
    List<DynamicList<label>>& srcAddr,
    List<DynamicList<scalar>>& srcWght,
    List<DynamicList<label>>& tgtAddr,
    List<DynamicList<scalar>>& tgtWght,
    label srcFacei,
    label tgtFacei
)
{
    addProfiling(ami, "faceAreaWeightAMI::calcAddressing");

    // construct weights and addressing
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    label nFacesRemaining = srcAddr.size();

    // list of tgt face neighbour faces
    DynamicList<label> nbrFaces(10);

    // list of faces currently visited for srcFacei to avoid multiple hits
    DynamicList<label> visitedFaces(10);

    // list to keep track of tgt faces used to seed src faces
    labelList seedFaces(nFacesRemaining, -1);
    seedFaces[srcFacei] = tgtFacei;

    // list to keep track of whether src face can be mapped
    boolList mapFlag(nFacesRemaining, true);

    // reset starting seed
    label startSeedi = 0;

    DynamicList<label> nonOverlapFaces;
    do
    {
        // Do advancing front starting from srcFacei,tgtFacei
        bool faceProcessed = processSourceFace
        (
            srcFacei,
            tgtFacei,

            nbrFaces,
            visitedFaces,

            srcAddr,
            srcWght,
            tgtAddr,
            tgtWght
        );

        mapFlag[srcFacei] = false;

        nFacesRemaining--;

        if (!faceProcessed)
        {
            nonOverlapFaces.append(srcFacei);
        }

        // choose new src face from current src face neighbour
        if (nFacesRemaining > 0)
        {
            setNextFaces
            (
                startSeedi,
                srcFacei,
                tgtFacei,
                mapFlag,
                seedFaces,
                visitedFaces
            );
        }
    } while (nFacesRemaining > 0);

    this->srcNonOverlap_.transfer(nonOverlapFaces);
}


template<class SourcePatch, class TargetPatch>
bool Foam::faceAreaWeightAMI<SourcePatch, TargetPatch>::processSourceFace
(
    const label srcFacei,
    const label tgtStartFacei,

    // list of tgt face neighbour faces
    DynamicList<label>& nbrFaces,
    // list of faces currently visited for srcFacei to avoid multiple hits
    DynamicList<label>& visitedFaces,

    // temporary storage for addressing and weights
    List<DynamicList<label>>& srcAddr,
    List<DynamicList<scalar>>& srcWght,
    List<DynamicList<label>>& tgtAddr,
    List<DynamicList<scalar>>& tgtWght
)
{
    addProfiling(ami, "faceAreaWeightAMI::processSourceFace");

    if (tgtStartFacei == -1)
    {
        return false;
    }

    nbrFaces.clear();
    visitedFaces.clear();

    // append initial target face and neighbours
    nbrFaces.append(tgtStartFacei);
    this->appendNbrFaces
    (
        tgtStartFacei,
        this->tgtPatch_,
        visitedFaces,
        nbrFaces
    );

    bool faceProcessed = false;

    label maxNeighbourFaces = nbrFaces.size();

    do
    {
        // process new target face
        label tgtFacei = nbrFaces.remove();
        visitedFaces.append(tgtFacei);
        scalar area = interArea(srcFacei, tgtFacei);

        // store when intersection fractional area > tolerance
        if (area/this->srcMagSf_[srcFacei] > faceAreaIntersect::tolerance())
        {
            srcAddr[srcFacei].append(tgtFacei);
            srcWght[srcFacei].append(area);

            tgtAddr[tgtFacei].append(srcFacei);
            tgtWght[tgtFacei].append(area);

            this->appendNbrFaces
            (
                tgtFacei,
                this->tgtPatch_,
                visitedFaces,
                nbrFaces
            );

            faceProcessed = true;

            maxNeighbourFaces = max(maxNeighbourFaces, nbrFaces.size());
        }

    } while (nbrFaces.size() > 0);

    if (debug > 1)
    {
        DebugVar(maxNeighbourFaces);
    }

    return faceProcessed;
}


template<class SourcePatch, class TargetPatch>
void Foam::faceAreaWeightAMI<SourcePatch, TargetPatch>::setNextFaces
(
    label& startSeedi,
    label& srcFacei,
    label& tgtFacei,
    const boolList& mapFlag,
    labelList& seedFaces,
    const DynamicList<label>& visitedFaces,
    bool errorOnNotFound
) const
{
    addProfiling(ami, "faceAreaWeightAMI::setNextFaces");

    const labelList& srcNbrFaces = this->srcPatch_.faceFaces()[srcFacei];

    // initialise tgtFacei
    tgtFacei = -1;

    // set possible seeds for later use
    bool valuesSet = false;
    for (label faceS: srcNbrFaces)
    {
        if (mapFlag[faceS] && seedFaces[faceS] == -1)
        {
            for (label faceT : visitedFaces)
            {
                const scalar threshold =
                    this->srcMagSf_[faceS]*faceAreaIntersect::tolerance();

                // store when intersection fractional area > tolerance
                // Check that faces have enough overlap for robust walking
                if (overlaps(faceS, faceT, threshold))
                {
                    seedFaces[faceS] = faceT;

                    if (!valuesSet)
                    {
                        srcFacei = faceS;
                        tgtFacei = faceT;
                        valuesSet = true;
                    }
                }
            }
        }
    }

    // set next src and tgt faces if not set above
    if (valuesSet)
    {
        return;
    }
    else
    {
        // try to use existing seed
        bool foundNextSeed = false;
        for (label facei = startSeedi; facei < mapFlag.size(); ++facei)
        {
            if (mapFlag[facei])
            {
                if (!foundNextSeed)
                {
                    startSeedi = facei;
                    foundNextSeed = true;
                }

                if (seedFaces[facei] != -1)
                {
                    srcFacei = facei;
                    tgtFacei = seedFaces[facei];

                    return;
                }
            }
        }

        // perform new search to find match
        if (debug)
        {
            Pout<< "Advancing front stalled: searching for new "
                << "target face" << endl;
        }

        foundNextSeed = false;
        for (label facei = startSeedi; facei < mapFlag.size(); ++facei)
        {
            if (mapFlag[facei])
            {
                if (!foundNextSeed)
                {
                    startSeedi = facei + 1;
                    foundNextSeed = true;
                }

                srcFacei = facei;
                tgtFacei = this->findTargetFace(srcFacei);

                if (tgtFacei >= 0)
                {
                    return;
                }
            }
        }

        if (errorOnNotFound)
        {
            FatalErrorInFunction
               << "Unable to set source and target faces" << abort(FatalError);
        }
    }
}


template<class SourcePatch, class TargetPatch>
Foam::scalar Foam::faceAreaWeightAMI<SourcePatch, TargetPatch>::interArea
(
    const label srcFacei,
    const label tgtFacei
) const
{
    addProfiling(ami, "faceAreaWeightAMI::interArea");

    scalar area = 0;

    const pointField& srcPoints = this->srcPatch_.points();
    const pointField& tgtPoints = this->tgtPatch_.points();

    // references to candidate faces
    const face& src = this->srcPatch_[srcFacei];
    const face& tgt = this->tgtPatch_[tgtFacei];

    // quick reject if either face has zero area
    if
    (
        (this->srcMagSf_[srcFacei] < ROOTVSMALL)
     || (this->tgtMagSf_[tgtFacei] < ROOTVSMALL)
    )
    {
        return area;
    }

    // create intersection object
    faceAreaIntersect inter
    (
        srcPoints,
        tgtPoints,
        this->srcTris_[srcFacei],
        this->tgtTris_[tgtFacei],
        this->reverseTarget_,
        AMIInterpolation<SourcePatch, TargetPatch>::cacheIntersections_
    );

    // crude resultant norm
    vector n(-this->srcPatch_.faceNormals()[srcFacei]);
    if (this->reverseTarget_)
    {
        n -= this->tgtPatch_.faceNormals()[tgtFacei];
    }
    else
    {
        n += this->tgtPatch_.faceNormals()[tgtFacei];
    }
    scalar magN = mag(n);

    if (magN > ROOTVSMALL)
    {
        area = inter.calc(src, tgt, n/magN);
    }
    else
    {
        WarningInFunction
            << "Invalid normal for source face " << srcFacei
            << " points " << UIndirectList<point>(srcPoints, src)
            << " target face " << tgtFacei
            << " points " << UIndirectList<point>(tgtPoints, tgt)
            << endl;
    }


    if ((debug > 1) && (area > 0))
    {
        this->writeIntersectionOBJ(area, src, tgt, srcPoints, tgtPoints);
    }

    return area;
}


template<class SourcePatch, class TargetPatch>
bool Foam::faceAreaWeightAMI<SourcePatch, TargetPatch>::overlaps
(
    const label srcFacei,
    const label tgtFacei,
    const scalar threshold
) const
{
    const pointField& srcPoints = this->srcPatch_.points();
    const pointField& tgtPoints = this->tgtPatch_.points();

    // references to candidate faces
    const face& src = this->srcPatch_[srcFacei];
    const face& tgt = this->tgtPatch_[tgtFacei];

    // quick reject if either face has zero area
    if
    (
        (this->srcMagSf_[srcFacei] < ROOTVSMALL)
     || (this->tgtMagSf_[tgtFacei] < ROOTVSMALL)
    )
    {
        return false;
    }

    faceAreaIntersect inter
    (
        srcPoints,
        tgtPoints,
        this->srcTris_[srcFacei],
        this->tgtTris_[tgtFacei],
        this->reverseTarget_,
        AMIInterpolation<SourcePatch, TargetPatch>::cacheIntersections_
    );

    // crude resultant norm
    vector n(-this->srcPatch_.faceNormals()[srcFacei]);
    if (this->reverseTarget_)
    {
        n -= this->tgtPatch_.faceNormals()[tgtFacei];
    }
    else
    {
        n += this->tgtPatch_.faceNormals()[tgtFacei];
    }
    scalar magN = mag(n);

    if (magN > ROOTVSMALL)
    {
        return inter.overlaps(src, tgt, n/magN, threshold);
    }
    else
    {
        WarningInFunction
            << "Invalid normal for source face " << srcFacei
            << " points " << UIndirectList<point>(srcPoints, src)
            << " target face " << tgtFacei
            << " points " << UIndirectList<point>(tgtPoints, tgt)
            << endl;
    }

    return false;
}


template<class SourcePatch, class TargetPatch>
void Foam::faceAreaWeightAMI<SourcePatch, TargetPatch>::
restartUncoveredSourceFace
(
    List<DynamicList<label>>& srcAddr,
    List<DynamicList<scalar>>& srcWght,
    List<DynamicList<label>>& tgtAddr,
    List<DynamicList<scalar>>& tgtWght
)
{
    addProfiling(ami, "faceAreaWeightAMI::restartUncoveredSourceFace");

    // Collect all src faces with a low weight
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    labelHashSet lowWeightFaces(100);
    forAll(srcWght, srcFacei)
    {
        scalar s = sum(srcWght[srcFacei]);
        scalar t = s/this->srcMagSf_[srcFacei];

        if (t < 0.5)
        {
            lowWeightFaces.insert(srcFacei);
        }
    }

    if (debug)
    {
        Pout<< "faceAreaWeightAMI: restarting search on "
            << lowWeightFaces.size() << " faces since sum of weights < 0.5"
            << endl;
    }

    if (lowWeightFaces.size() > 0)
    {
        // Erase all the lowWeight source faces from the target
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        DynamicList<label> okSrcFaces(10);
        DynamicList<scalar> okSrcWeights(10);
        forAll(tgtAddr, tgtFacei)
        {
            okSrcFaces.clear();
            okSrcWeights.clear();
            DynamicList<label>& srcFaces = tgtAddr[tgtFacei];
            DynamicList<scalar>& srcWeights = tgtWght[tgtFacei];
            forAll(srcFaces, i)
            {
                if (!lowWeightFaces.found(srcFaces[i]))
                {
                    okSrcFaces.append(srcFaces[i]);
                    okSrcWeights.append(srcWeights[i]);
                }
            }
            if (okSrcFaces.size() < srcFaces.size())
            {
                srcFaces.transfer(okSrcFaces);
                srcWeights.transfer(okSrcWeights);
            }
        }



        // Restart search from best hit
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        // list of tgt face neighbour faces
        DynamicList<label> nbrFaces(10);

        // list of faces currently visited for srcFacei to avoid multiple hits
        DynamicList<label> visitedFaces(10);

        for (const label srcFacei : lowWeightFaces)
        {
            const label tgtFacei = this->findTargetFace(srcFacei);
            if (tgtFacei != -1)
            {
                //bool faceProcessed =
                processSourceFace
                (
                    srcFacei,
                    tgtFacei,

                    nbrFaces,
                    visitedFaces,

                    srcAddr,
                    srcWght,
                    tgtAddr,
                    tgtWght
                );
                // ? Check faceProcessed to see if restarting has worked.
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class SourcePatch, class TargetPatch>
Foam::faceAreaWeightAMI<SourcePatch, TargetPatch>::faceAreaWeightAMI
(
    const SourcePatch& srcPatch,
    const TargetPatch& tgtPatch,
    const faceAreaIntersect::triangulationMode& triMode,
    const bool reverseTarget,
    const bool requireMatch,
    const bool restartUncoveredSourceFace
)
:
    AMIMethod<SourcePatch, TargetPatch>
    (
        srcPatch,
        tgtPatch,
        triMode,
        reverseTarget,
        requireMatch
    ),
    restartUncoveredSourceFace_(restartUncoveredSourceFace),
    srcTris_(),
    tgtTris_()
{
    this->triangulatePatch(srcPatch, srcTris_, this->srcMagSf_);
    this->triangulatePatch(tgtPatch, tgtTris_, this->tgtMagSf_);
}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

template<class SourcePatch, class TargetPatch>
Foam::faceAreaWeightAMI<SourcePatch, TargetPatch>::~faceAreaWeightAMI()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class SourcePatch, class TargetPatch>
void Foam::faceAreaWeightAMI<SourcePatch, TargetPatch>::calculate
(
    labelListList& srcAddress,
    scalarListList& srcWeights,
    labelListList& tgtAddress,
    scalarListList& tgtWeights,
    label srcFacei,
    label tgtFacei
)
{
    addProfiling(ami, "faceAreaWeightAMI::calculate");

    bool ok =
        this->initialise
        (
            srcAddress,
            srcWeights,
            tgtAddress,
            tgtWeights,
            srcFacei,
            tgtFacei
        );

    if (!ok)
    {
        return;
    }

    // temporary storage for addressing and weights
    List<DynamicList<label>> srcAddr(this->srcPatch_.size());
    List<DynamicList<scalar>> srcWght(srcAddr.size());
    List<DynamicList<label>> tgtAddr(this->tgtPatch_.size());
    List<DynamicList<scalar>> tgtWght(tgtAddr.size());

    calcAddressing
    (
        srcAddr,
        srcWght,
        tgtAddr,
        tgtWght,
        srcFacei,
        tgtFacei
    );

    if (debug && !this->srcNonOverlap_.empty())
    {
        Pout<< "    AMI: " << this->srcNonOverlap_.size()
            << " non-overlap faces identified"
            << endl;
    }


    // Check for badly covered faces
    if (restartUncoveredSourceFace_)
    {
        restartUncoveredSourceFace
        (
            srcAddr,
            srcWght,
            tgtAddr,
            tgtWght
        );
    }


    // transfer data to persistent storage
    forAll(srcAddr, i)
    {
        srcAddress[i].transfer(srcAddr[i]);
        srcWeights[i].transfer(srcWght[i]);
    }
    forAll(tgtAddr, i)
    {
        tgtAddress[i].transfer(tgtAddr[i]);
        tgtWeights[i].transfer(tgtWght[i]);
    }
}


template<class SourcePatch, class TargetPatch>
void Foam::faceAreaWeightAMI<SourcePatch, TargetPatch>::setMagSf
(
    const TargetPatch& tgtPatch,
    const mapDistribute& map,
    scalarList& srcMagSf,
    scalarList& tgtMagSf
) const
{
    srcMagSf = std::move(this->srcMagSf_);
    tgtMagSf = std::move(this->tgtMagSf_);
    map.reverseDistribute(tgtPatch.size(), tgtMagSf);
}


template<class SourcePatch, class TargetPatch>
void Foam::faceAreaWeightAMI<SourcePatch, TargetPatch>::normaliseWeights
(
    const bool verbose,
    AMIInterpolation<SourcePatch, TargetPatch>& inter
)
{
    inter.normaliseWeights(this->conformal(), verbose);
}


// ************************************************************************* //
