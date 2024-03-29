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
    Foam::primalSolver

Description
    Base class for primal solvers

\*---------------------------------------------------------------------------*/

#ifndef primalSolver_H
#define primalSolver_H

#include "fvMesh.H"
#include "solver.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                            Class primalSolver Declaration
\*---------------------------------------------------------------------------*/

class primalSolver
:
    public solver
{
private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        primalSolver(const primalSolver&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const primalSolver&) = delete;


public:


    // Static Data Members

        //- Run-time type information
        TypeName("primalSolver");


    // Declare run-time constructor selection table

        declareRunTimeNewSelectionTable
        (
            autoPtr,
            primalSolver,
            primalSolver,
            (
                fvMesh& mesh,
                const word& managerType,
                const dictionary& dict
            ),
            (mesh, managerType, dict)
        );


    // Constructors

        //- Construct from mesh and dictionary
        primalSolver
        (
            fvMesh& mesh,
            const word& managerType,
            const dictionary& dict
        );


    // Selectors

        //- Return a reference to the selected primal solver
        static autoPtr<primalSolver> New
        (
            fvMesh& mesh,
            const word& managerType,
            const dictionary& dict
        );


    //- Destructor
    virtual ~primalSolver() = default;


    // Member Functions

        virtual bool readDict(const dictionary& dict);


        // Evolution

            //- update boundary conditions
            virtual void correctBoundaryConditions();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
