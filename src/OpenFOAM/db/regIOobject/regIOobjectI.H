/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2010, 2018-2019 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2011-2015 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool Foam::regIOobject::ownedByRegistry() const
{
    return ownedByRegistry_;
}


inline void Foam::regIOobject::store()
{
    ownedByRegistry_ = true;
}


template<class Type>
inline Type& Foam::regIOobject::store(Type* p)
{
    if (!p)
    {
        FatalErrorInFunction
            << "object deallocated"
            << abort(FatalError);
    }

    p->regIOobject::ownedByRegistry_ = true;

    return *p;
}


template<class Type>
inline Type& Foam::regIOobject::store(autoPtr<Type>& aptr)
{
    return store(aptr.ptr()); // release, pass management to regIOobject
}


template<class Type>
inline Type& Foam::regIOobject::store(autoPtr<Type>&& aptr)
{
    return store(aptr.ptr()); // release, pass management to regIOobject
}


inline void Foam::regIOobject::release(const bool unregister)
{
    ownedByRegistry_ = false;
    if (unregister)
    {
        registered_ = false;
    }
}


inline Foam::label Foam::regIOobject::eventNo() const
{
    return eventNo_;
}

inline Foam::label& Foam::regIOobject::eventNo()
{
    return eventNo_;
}


inline const Foam::labelList& Foam::regIOobject::watchIndices() const
{
    return watchIndices_;
}


inline Foam::labelList& Foam::regIOobject::watchIndices()
{
    return watchIndices_;
}


// ************************************************************************* //
