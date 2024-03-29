/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
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

// OpenFOAM includes
#include "geometryPatches.H"
#include "fvMesh.H"
#include "volFields.H"
#include "emptyPolyPatch.H"
#include "processorPolyPatch.H"
#include "foamVtkTools.H"
#include "runTimePostProcessing.H"

// VTK includes
#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkCompositeDataSet.h"
#include "vtkMultiPieceDataSet.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

vtkSmartPointer<vtkMultiPieceDataSet>
Foam::functionObjects::runTimePostPro::geometryPatches::gatherPatchPieces
(
    const labelListList& patchIds
) const
{
    const polyBoundaryMesh& patches = parent().mesh().boundaryMesh();

    label nPieces = 0;
    for (const labelList& ids : patchIds)
    {
        nPieces += ids.size();
    }

    auto multiPiece = vtkSmartPointer<vtkMultiPieceDataSet>::New();
    multiPiece->SetNumberOfPieces(nPieces);

    label pieceId = 0;

    if (!needsCollective())
    {
        // Simple case

        for (int proci=0; proci < Pstream::myProcNo(); ++proci)
        {
            pieceId += patchIds[proci].size();
        }

        for (const label patchId : patchIds[Pstream::myProcNo()])
        {
            const polyPatch& pp = patches[patchId];

            multiPiece->SetPiece
            (
                pieceId,
                Foam::vtk::Tools::Patch::mesh(pp)
            );

            ++pieceId;
        }
    }
    else if (Pstream::master())
    {
        // Gather pieces on master

        // Add myself
        for (const label patchId : patchIds[Pstream::myProcNo()])
        {
            const polyPatch& pp = patches[patchId];

            multiPiece->SetPiece
            (
                pieceId,
                Foam::vtk::Tools::Patch::mesh(pp)
            );

            ++pieceId;
        }

        // Receive surfaces
        pointField points;
        faceList faces;

        for
        (
            int slave=Pstream::firstSlave();
            slave<=Pstream::lastSlave();
            ++slave
        )
        {
            const label nSlavePatches = patchIds[slave].size();

            if (!nSlavePatches)
            {
                continue;
            }

            IPstream fromSlave(Pstream::commsTypes::scheduled, slave);

            for (label recvi=0; recvi < nSlavePatches; ++recvi)
            {
                points.clear();
                faces.clear();

                fromSlave
                    >> points >> faces;

                multiPiece->SetPiece
                (
                    pieceId,
                    Foam::vtk::Tools::Patch::mesh(points, faces)
                );

                ++pieceId;
            }
        }
    }
    else
    {
        // Slave - send surfaces
        const labelList& slavePatchIds = patchIds[Pstream::myProcNo()];

        if (slavePatchIds.size())
        {
            OPstream toMaster
            (
                Pstream::commsTypes::scheduled,
                Pstream::masterNo()
            );

            for (const label patchId : patchIds[Pstream::myProcNo()])
            {
                const polyPatch& pp = patches[patchId];

                toMaster
                    << pp.localPoints() << pp.localFaces();
            }
        }
    }

    return multiPiece;
}


vtkSmartPointer<vtkMultiPieceDataSet>
Foam::functionObjects::runTimePostPro::geometryPatches::gatherPatchFaceCentres
(
    const labelListList& patchIds
) const
{
    const polyBoundaryMesh& patches = parent().mesh().boundaryMesh();

    label nPieces = 0;
    for (const labelList& ids : patchIds)
    {
        nPieces += ids.size();
    }

    auto multiPiece = vtkSmartPointer<vtkMultiPieceDataSet>::New();
    multiPiece->SetNumberOfPieces(nPieces);

    label pieceId = 0;

    if (!needsCollective())
    {
        // Simple case

        for (int proci=0; proci < Pstream::myProcNo(); ++proci)
        {
            pieceId += patchIds[proci].size();
        }

        for (const label patchId : patchIds[Pstream::myProcNo()])
        {
            const polyPatch& pp = patches[patchId];

            auto geom = vtkSmartPointer<vtkPolyData>::New();

            geom->SetPoints(Foam::vtk::Tools::Patch::faceCentres(pp));
            geom->SetVerts(Foam::vtk::Tools::identityVertices(pp.size()));

            multiPiece->SetPiece(pieceId, geom);

            ++pieceId;
        }
    }
    else if (Pstream::master())
    {
        // Gather pieces (face centres) on master

        // Add myself
        for (const label patchId : patchIds[Pstream::myProcNo()])
        {
            const polyPatch& pp = patches[patchId];

            auto geom = vtkSmartPointer<vtkPolyData>::New();

            geom->SetPoints(Foam::vtk::Tools::Patch::faceCentres(pp));
            geom->SetVerts(Foam::vtk::Tools::identityVertices(pp.size()));

            multiPiece->SetPiece(pieceId, geom);

            ++pieceId;
        }

        // Receive points
        pointField points;

        for
        (
            int slave=Pstream::firstSlave();
            slave<=Pstream::lastSlave();
            ++slave
        )
        {
            const label nSlavePatches = patchIds[slave].size();

            if (!nSlavePatches)
            {
                continue;
            }

            IPstream fromSlave(Pstream::commsTypes::scheduled, slave);

            for (label recvi=0; recvi < nSlavePatches; ++recvi)
            {
                points.clear();

                fromSlave >> points;

                multiPiece->SetPiece
                (
                    pieceId,
                    Foam::vtk::Tools::Vertices(points)
                );

                ++pieceId;
            }
        }
    }
    else
    {
        // Slave - send face centres

        const labelList& slavePatchIds = patchIds[Pstream::myProcNo()];

        if (slavePatchIds.size())
        {
            OPstream toMaster
            (
                Pstream::commsTypes::scheduled,
                Pstream::masterNo()
            );

            for (const label patchId : patchIds[Pstream::myProcNo()])
            {
                const polyPatch& pp = patches[patchId];

                toMaster << pp.faceCentres();
            }
        }
    }

    return multiPiece;
}


// ************************************************************************* //
