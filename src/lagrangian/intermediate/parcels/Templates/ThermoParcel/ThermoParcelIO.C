/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2011, 2016-2019 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "ThermoParcel.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class ParcelType>
Foam::string Foam::ThermoParcel<ParcelType>::propertyList_ =
    Foam::ThermoParcel<ParcelType>::propertyList();

template<class ParcelType>
Foam::string Foam::ThermoParcel<ParcelType>::propertyTypes_ =
    Foam::ThermoParcel<ParcelType>::propertyTypes();

template<class ParcelType>
const std::size_t Foam::ThermoParcel<ParcelType>::sizeofFields
(
    sizeof(ThermoParcel<ParcelType>)
  - offsetof(ThermoParcel<ParcelType>, T_)
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::ThermoParcel<ParcelType>::ThermoParcel
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields,
    bool newFormat
)
:
    ParcelType(mesh, is, readFields, newFormat),
    T_(0.0),
    Cp_(0.0)
{
    if (readFields)
    {
        if (is.format() == IOstream::ASCII)
        {
            T_ = readScalar(is);
            Cp_ = readScalar(is);
        }
        else
        {
            is.read(reinterpret_cast<char*>(&T_), sizeofFields);
        }
    }

    is.check(FUNCTION_NAME);
}


template<class ParcelType>
template<class CloudType>
void Foam::ThermoParcel<ParcelType>::readFields(CloudType& c)
{
    bool valid = c.size();

    ParcelType::readFields(c);

    IOField<scalar> T(c.fieldIOobject("T", IOobject::MUST_READ), valid);
    c.checkFieldIOobject(c, T);

    IOField<scalar> Cp(c.fieldIOobject("Cp", IOobject::MUST_READ), valid);
    c.checkFieldIOobject(c, Cp);


    label i = 0;
    for (ThermoParcel<ParcelType>& p : c)
    {
        p.T_ = T[i];
        p.Cp_ = Cp[i];
        ++i;
    }
}


template<class ParcelType>
template<class CloudType>
void Foam::ThermoParcel<ParcelType>::writeFields(const CloudType& c)
{
    ParcelType::writeFields(c);

    label np = c.size();

    IOField<scalar> T(c.fieldIOobject("T", IOobject::NO_READ), np);
    IOField<scalar> Cp(c.fieldIOobject("Cp", IOobject::NO_READ), np);

    label i = 0;
    for (const ThermoParcel<ParcelType>& p : c)
    {
        T[i] = p.T_;
        Cp[i] = p.Cp_;
        ++i;
    }

    T.write(np > 0);
    Cp.write(np > 0);
}


template<class ParcelType>
template<class CloudType>
void Foam::ThermoParcel<ParcelType>::writeObjects
(
    const CloudType& c,
    objectRegistry& obr
)
{
    ParcelType::writeObjects(c, obr);

    label np = c.size();

    IOField<scalar>& T(cloud::createIOField<scalar>("T", np, obr));
    IOField<scalar>& Cp(cloud::createIOField<scalar>("Cp", np, obr));

    label i = 0;
    for (const ThermoParcel<ParcelType>& p : c)
    {
        T[i] = p.T_;
        Cp[i] = p.Cp_;
        ++i;
    }
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class ParcelType>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const ThermoParcel<ParcelType>& p
)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const ParcelType&>(p)
            << token::SPACE << p.T()
            << token::SPACE << p.Cp();
    }
    else
    {
        os  << static_cast<const ParcelType&>(p);
        os.write
        (
            reinterpret_cast<const char*>(&p.T_),
            ThermoParcel<ParcelType>::sizeofFields
        );
    }

    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //
