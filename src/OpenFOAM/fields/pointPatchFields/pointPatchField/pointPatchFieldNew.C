/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         |2011 OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2011-2016 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::autoPtr<Foam::pointPatchField<Type>> Foam::pointPatchField<Type>::New
(
    const word& patchFieldType,
    const word& actualPatchType,
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF
)
{
    if (debug)
    {
        InfoInFunction << "Constructing pointPatchField<Type>" << endl;
    }

    auto cstrIter = pointPatchConstructorTablePtr_->cfind(patchFieldType);

    if (!cstrIter.found())
    {
        FatalErrorInFunction
            << "Unknown patchFieldType type "
            << patchFieldType << nl << nl
            << "Valid patchField types :" << endl
            << pointPatchConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    autoPtr<pointPatchField<Type>> pfPtr(cstrIter()(p, iF));

    if
    (
        actualPatchType == word::null
     || actualPatchType != p.type()
    )
    {
        if (pfPtr().constraintType() != p.constraintType())
        {
            // Use default constraint type
            auto patchTypeCstrIter =
                pointPatchConstructorTablePtr_->cfind(p.type());

            if (!patchTypeCstrIter.found())
            {
                FatalErrorInFunction
                    << "inconsistent patch and patchField types for \n"
                    << "    patch type " << p.type()
                    << " and patchField type " << patchFieldType
                    << exit(FatalError);
            }

            return patchTypeCstrIter()(p, iF);
        }
    }
    else
    {
        if (pointPatchConstructorTablePtr_->found(p.type()))
        {
            pfPtr().patchType() = actualPatchType;
        }
    }

    return pfPtr;
}


template<class Type>
Foam::autoPtr<Foam::pointPatchField<Type>> Foam::pointPatchField<Type>::New
(
    const word& patchFieldType,
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF
)
{
    return New(patchFieldType, word::null, p, iF);
}


template<class Type>
Foam::autoPtr<Foam::pointPatchField<Type>> Foam::pointPatchField<Type>::New
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const dictionary& dict
)
{
    if (debug)
    {
        InfoInFunction << "Constructing pointPatchField<Type>" << endl;
    }

    const word patchFieldType(dict.get<word>("type"));

    auto cstrIter = dictionaryConstructorTablePtr_->cfind(patchFieldType);

    if (!cstrIter.found())
    {
        if (!disallowGenericPointPatchField)
        {
            cstrIter = dictionaryConstructorTablePtr_->cfind("generic");
        }

        if (!cstrIter.found())
        {
            FatalIOErrorInFunction(dict)
                << "Unknown patchField type " << patchFieldType
                << " for patch type " << p.type() << nl << nl
                << "Valid patchField types :" << endl
                << dictionaryConstructorTablePtr_->sortedToc()
                << exit(FatalIOError);
        }
    }

    // Construct (but not necessarily returned)
    autoPtr<pointPatchField<Type>> pfPtr(cstrIter()(p, iF, dict));

    if
    (
       !dict.found("patchType")
     || dict.get<word>("patchType") != p.type()
    )
    {
        if (pfPtr().constraintType() == p.constraintType())
        {
            // Compatible (constraint-wise) with the patch type
            return pfPtr;
        }
        else
        {
            // Use default constraint type
            auto patchTypeCstrIter
                = dictionaryConstructorTablePtr_->cfind(p.type());

            if (!patchTypeCstrIter.found())
            {
                FatalIOErrorInFunction(dict)
                    << "inconsistent patch and patchField types for \n"
                    << "    patch type " << p.type()
                    << " and patchField type " << patchFieldType
                    << exit(FatalIOError);
            }

            return patchTypeCstrIter()(p, iF, dict);
        }
    }

    return cstrIter()(p, iF, dict);
}


template<class Type>
Foam::autoPtr<Foam::pointPatchField<Type>> Foam::pointPatchField<Type>::New
(
    const pointPatchField<Type>& ptf,
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const pointPatchFieldMapper& pfMapper
)
{
    if (debug)
    {
        InfoInFunction << "Constructing pointPatchField<Type>" << endl;
    }

    auto cstrIter = patchMapperConstructorTablePtr_->cfind(ptf.type());

    if (!cstrIter.found())
    {
        FatalErrorInFunction
            << "Unknown patchField type "
            << ptf.type() << nl << nl
            << "Valid patchField types :" << endl
            << patchMapperConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return cstrIter()(ptf, p, iF, pfMapper);
}


// ************************************************************************* //
