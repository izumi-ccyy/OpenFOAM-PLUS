/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2019 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2013-2017 OpenFOAM Foundation
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

#include "MPPICParcel.H"
#include "IOstreams.H"
#include "IOField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class ParcelType>
Foam::string Foam::MPPICParcel<ParcelType>::propertyList_ =
    Foam::MPPICParcel<ParcelType>::propertyList();

template<class ParcelType>
Foam::string Foam::MPPICParcel<ParcelType>::propertyTypes_ =
    Foam::MPPICParcel<ParcelType>::propertyTypes();

template<class ParcelType>
const std::size_t Foam::MPPICParcel<ParcelType>::sizeofFields
(
    sizeof(MPPICParcel<ParcelType>) - sizeof(ParcelType)
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::MPPICParcel<ParcelType>::MPPICParcel
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields,
    bool newFormat
)
:
    ParcelType(mesh, is, readFields, newFormat),
    UCorrect_(Zero)
{
    if (readFields)
    {
        if (is.format() == IOstream::ASCII)
        {
            is >> UCorrect_;
        }
        else
        {
            is.read(reinterpret_cast<char*>(&UCorrect_), sizeofFields);
        }
    }

    is.check(FUNCTION_NAME);
}


template<class ParcelType>
template<class CloudType>
void Foam::MPPICParcel<ParcelType>::readFields(CloudType& c)
{
    bool valid = c.size();

    ParcelType::readFields(c);

    IOField<vector> UCorrect
    (
        c.fieldIOobject("UCorrect", IOobject::MUST_READ),
        valid
    );
    c.checkFieldIOobject(c, UCorrect);

    label i = 0;
    for (MPPICParcel<ParcelType>& p : c)
    {
        p.UCorrect_ = UCorrect[i];

        ++i;
    }
}


template<class ParcelType>
template<class CloudType>
void Foam::MPPICParcel<ParcelType>::writeFields(const CloudType& c)
{
    ParcelType::writeFields(c);

    label np = c.size();

    IOField<vector>
        UCorrect(c.fieldIOobject("UCorrect", IOobject::NO_READ), np);

    label i = 0;

    for (const MPPICParcel<ParcelType>& p : c)
    {
        UCorrect[i] = p.UCorrect();

        ++i;
    }

    UCorrect.write(np > 0);
}


template<class ParcelType>
template<class CloudType>
void Foam::MPPICParcel<ParcelType>::writeObjects
(
    const CloudType& c,
    objectRegistry& obr
)
{
    ParcelType::writeObjects(c, obr);

    label np = c.size();

    IOField<vector>&
        UCorrect(cloud::createIOField<vector>("UCorrect", np, obr));

    label i = 0;

    for (const MPPICParcel<ParcelType>& p : c)
    {
        UCorrect[i] = p.UCorrect();

        ++i;
    }
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class ParcelType>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const MPPICParcel<ParcelType>& p
)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const ParcelType&>(p)
            << token::SPACE << p.UCorrect();
    }
    else
    {
        os  << static_cast<const ParcelType&>(p);
        os.write
        (
            reinterpret_cast<const char*>(&p.UCorrect_),
            MPPICParcel<ParcelType>::sizeofFields
        );
    }

    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //
