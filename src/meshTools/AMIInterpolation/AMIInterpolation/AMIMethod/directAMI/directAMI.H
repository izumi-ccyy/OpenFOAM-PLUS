/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
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

Class
    Foam::directAMI

Description
    Direct mapped Arbitrary Mesh Interface (AMI) method

SourceFiles
    directAMI.C

\*---------------------------------------------------------------------------*/

#ifndef directAMI_H
#define directAMI_H

#include "AMIMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class directAMI Declaration
\*---------------------------------------------------------------------------*/

template<class SourcePatch, class TargetPatch>
class directAMI
:
    public AMIMethod<SourcePatch, TargetPatch>
{

private:

    // Private Member Functions

        //- No copy construct
        directAMI(const directAMI&) = delete;

        //- No copy assignment
        void operator=(const directAMI&) = delete;

        // Marching front

            //- Append to list of src face seed indices
            void appendToDirectSeeds
            (
                labelList& mapFlag,
                labelList& srcTgtSeed,
                DynamicList<label>& srcSeeds,
                DynamicList<label>& nonOverlapFaces,
                label& srcFacei,
                label& tgtFacei
            ) const;

            //- Restart the advancing front - typically happens for
            //  disconnected regions
            void restartAdvancingFront
            (
                labelList& mapFlag,
                DynamicList<label>& nonOverlapFaces,
                label& srcFacei,
                label& tgtFacei
            ) const;


        // Evaluation

            //- Area of intersection between source and target faces
            scalar interArea
            (
                const label srcFacei,
                const label tgtFacei
            ) const;


public:

    //- Runtime type information
    TypeName("directAMI");


    // Constructors

        //- Construct from components
        directAMI
        (
            const SourcePatch& srcPatch,
            const TargetPatch& tgtPatch,
            const faceAreaIntersect::triangulationMode& triMode,
            const bool reverseTarget = false,
            const bool requireMatch = true
        );


    //- Destructor
    virtual ~directAMI();


    // Member Functions

        // Manipulation

            //- Update addressing and weights
            virtual void calculate
            (
                labelListList& srcAddress,
                scalarListList& srcWeights,
                labelListList& tgtAddress,
                scalarListList& tgtWeights,
                label srcFacei = -1,
                label tgtFacei = -1
            );

            //- Set the face areas for parallel runs
            virtual void setMagSf
            (
                const TargetPatch& tgtPatch,
                const mapDistribute& map,
                scalarList& srcMagSf,
                scalarList& tgtMagSf
            ) const;

            //- Normalise the weight. Can optionally subset addressing
            //  (e.g. for mapNearest)
            virtual void normaliseWeights
            (
                const bool verbose,
                AMIInterpolation<SourcePatch, TargetPatch>& inter
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "directAMI.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
