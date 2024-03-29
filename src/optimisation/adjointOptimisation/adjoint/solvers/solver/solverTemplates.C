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

#include "fvmSup.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::solver::solveWithArgs
(
    Type& type,
    List<void (Type::*)()>& funcs
)
{
    // Iterate
    if (active_)
    {
        restoreInitValues();
        while(loop())
        {
            solveIter();
            forAll(funcs, fI)
            {
                (type.*funcs[fI])();
            }
        }

        mesh_.time().printExecutionTime(Info);
    }
}


template<class Type>
void Foam::solver::addOptimisationTypeSource
(
    fvMatrix<Type>& matrix
) const
{
    // If source has been allocated, add source * variable
    if (optTypeSource_)
    {
        const GeometricField<Type, fvPatchField, volMesh>& psi = matrix.psi();

        matrix += fvm::Sp(*optTypeSource_, psi);
    }
}


// ************************************************************************* //
