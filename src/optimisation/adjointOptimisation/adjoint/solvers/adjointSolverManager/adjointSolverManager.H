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
    Foam::adjointSolverManager

Description
    Class for managing adjoint solvers, which may be more than one per
    operating point

SourceFiles
    adjointSolverManager.C

\*---------------------------------------------------------------------------*/

#ifndef adjointSolverManager_H
#define adjointSolverManager_H

#include "adjointSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class adjointSolverManager Declaration
\*---------------------------------------------------------------------------*/

class adjointSolverManager
:
    public regIOobject
{
private:

    // Private Member Functions

        //- No  copy construct
        adjointSolverManager(const adjointSolverManager&) = delete;

        //- No copy assignment
        void operator=(const adjointSolverManager&) = delete;


protected:

    // Protected Data

        fvMesh& mesh_;

        dictionary dict_;

        const word managerName_;

        const word primalSolverName_;

        PtrList<adjointSolver> adjointSolvers_;

        labelList objectiveSolverIDs_;

        labelList constraintSolverIDs_;

        scalar operatingPointWeight_;


public:

    TypeName("adjointSolverManager");


    // Constructors

        //- Construct from components
        adjointSolverManager
        (
            fvMesh& mesh,
            const word& managerType,
            const dictionary& dict
        );


    //- Destructor
    virtual ~adjointSolverManager() = default;


    // Member Functions

        virtual bool readDict(const dictionary& dict);


        // Access

            //- Const access to the manager name
            const word& managerName() const;

            //- Const access to the construction dictionary
            const dictionary& dict() const;

            //- Const access to adjoint solvers
            const PtrList<adjointSolver>& adjointSolvers() const;

            //- Non-const access to adjoint solvers
            PtrList<adjointSolver>& adjointSolvers();

            //- Const access to adjoint solvers
            scalar operatingPointWeight() const;

            //- Number of adjoint solvers corresponding to contraints
            label nConstraints() const;

            //- Number of adjoint solvers corresponding to objectives
            label nObjectives() const;

            //- Total number of adjoint solvers
            label nAdjointSolvers() const;


        // Evolution

            //- Update objective function-related values and solve adjoint
            //- equations
            virtual void solveAdjointEquations();

            //- Aggregate sensitivities from various adjoint solvers
            virtual tmp<scalarField> aggregateSensitivities();

            //- Get constraint sensitivities. One scalarField per constraint
            virtual PtrList<scalarField> constraintSensitivities();

            //- Compute sensitivities for all adjoint solvers
            //- (both objective- and constraint-related ones)
            //- Clears previous sensitivity fields
            void computeAllSensitivities();

            //- Get objective value
            scalar objectiveValue();

            //- Get constraint values
            virtual tmp<scalarField> constraintValues();

            //- Update fields related to primal solution.
            //  For instance, primal fields of adjoint turbulence models
            void updatePrimalBasedQuantities(const word& name);


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
