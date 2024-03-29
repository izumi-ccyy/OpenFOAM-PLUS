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


Class
    Foam::objectiveManager

Description
    class for managing incompressible objective functions.

SourceFiles
    objectiveManager.C

\*---------------------------------------------------------------------------*/

#ifndef objectiveManager_H
#define objectiveManager_H

#include "fvMesh.H"
#include "objective.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class objectiveManager Declaration
\*---------------------------------------------------------------------------*/

class objectiveManager
:
    public regIOobject
{
protected:

    // Protected data

        const fvMesh& mesh_;
        const dictionary& dict_;
        const word adjointSolverName_;
        const word primalSolverName_;
        PtrList<objective> objectives_;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        objectiveManager(const objectiveManager&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const objectiveManager&) = delete;


public:

    TypeName("objectiveManager");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            objectiveManager,
            dictionary,
            (
                const fvMesh& mesh,
                const dictionary& dict,
                const word& adjointSolverName,
                const word& primalSolverName
            ),
            (mesh, dict, adjointSolverName, primalSolverName)
        );

    // Constructors

        //- Construct from components
        objectiveManager
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const word& adjointSolverName,
            const word& primalSolverName
        );

    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<objectiveManager> New
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const word& adjointSolverName,
            const word& primalSolverName
        );


    //- Destructor
    virtual ~objectiveManager() = default;


    // Member Functions

        virtual bool readDict(const dictionary& dict);

        //- Update objective function related values
        void updateNormalizationFactor();

        //- Update objective function related values
        void update();

        //- Print to screen
        scalar print();

        //- Write objective function history
        virtual bool write(const bool valid = true) const;

        //- Call all functions required prior to the solution of the adjoint
        //- equations
        void updateAndWrite();

        //- Return reference to objective functions
        PtrList<objective>& getObjectiveFunctions();

        //- Return constant reference to objective functions
        const PtrList<objective>& getObjectiveFunctions() const;

        //- Return name of adjointSolverManager
        const word& adjointSolverManagerName() const;

        //- Return name of the adjointSolver
        const word& adjointSolverName() const;

        //- Return name of the primalSolver
        const word& primalSolverName() const;

        //- Add contribution to adjoint momentum PDEs
        virtual void addUaEqnSource(fvVectorMatrix& UaEqn) = 0;

        //- Add contribution to adjoint momentum PDEs
        virtual void addPaEqnSource(fvScalarMatrix& paEqn) = 0;

        //- Add contribution to first adjoint turbulence model PDE
        virtual void addTMEqn1Source(fvScalarMatrix& adjTMEqn1) = 0;

        //- Add contribution to second adjoint turbulence model PDE
        virtual void addTMEqn2Source(fvScalarMatrix& adjTMEqn2) = 0;


    // IO

        virtual bool writeData(Ostream&) const
        {
            return true;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
