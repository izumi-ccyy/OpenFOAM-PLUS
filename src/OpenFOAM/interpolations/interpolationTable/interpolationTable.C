/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2011 OpenCFD Ltd.
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

#include "interpolationTable.H"
#include "openFoamTableReader.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class Type>
void Foam::interpolationTable<Type>::readTable()
{
    // preserve the original (unexpanded) fileName to avoid absolute paths
    // appearing subsequently in the write() method
    fileName fName(fileName_);

    fName.expand();

    // Read data from file
    reader_()(fName, *this);

    if (this->empty())
    {
        FatalErrorInFunction
            << "table read from " << fName << " is empty" << nl
            << exit(FatalError);
    }

    // Check that the data are okay
    check();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::interpolationTable<Type>::interpolationTable()
:
    List<Tuple2<scalar, Type>>(),
    bounding_(bounds::repeatableBounding::WARN),
    fileName_("fileNameIsUndefined"),
    reader_(nullptr)
{}


template<class Type>
Foam::interpolationTable<Type>::interpolationTable
(
    const List<Tuple2<scalar, Type>>& values,
    const bounds::repeatableBounding bounding,
    const fileName& fName
)
:
    List<Tuple2<scalar, Type>>(values),
    bounding_(bounding),
    fileName_(fName),
    reader_(nullptr)
{}


template<class Type>
Foam::interpolationTable<Type>::interpolationTable(const fileName& fName)
:
    List<Tuple2<scalar, Type>>(),
    bounding_(bounds::repeatableBounding::WARN),
    fileName_(fName),
    reader_(new openFoamTableReader<Type>(dictionary()))
{
    readTable();
}


template<class Type>
Foam::interpolationTable<Type>::interpolationTable(const dictionary& dict)
:
    List<Tuple2<scalar, Type>>(),
    bounding_
    (
        bounds::repeatableBoundingNames.lookupOrDefault
        (
            "outOfBounds",
            dict,
            bounds::repeatableBounding::WARN,
            true  // Failsafe behaviour
        )
    ),
    fileName_(dict.lookup("file")),
    reader_(tableReader<Type>::New(dict))
{
    readTable();
}


template<class Type>
Foam::interpolationTable<Type>::interpolationTable
(
     const interpolationTable& interpTable
)
:
    List<Tuple2<scalar, Type>>(interpTable),
    bounding_(interpTable.bounding_),
    fileName_(interpTable.fileName_),
    reader_(interpTable.reader_)    // note: steals reader. Used in write().
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::interpolationTable<Type>::check() const
{
    const label n = this->size();
    scalar prevValue = this->first().first();

    for (label i=1; i<n; ++i)
    {
        const scalar currValue =
            List<Tuple2<scalar, Type>>::operator[](i).first();

        // avoid duplicate values (divide-by-zero error)
        if (currValue <= prevValue)
        {
            FatalErrorInFunction
                << "out-of-order value: "
                << currValue << " at index " << i << nl
                << exit(FatalError);
        }
        prevValue = currValue;
    }
}


template<class Type>
void Foam::interpolationTable<Type>::write(Ostream& os) const
{
    os.writeEntry("file", fileName_);
    os.writeEntry("outOfBounds", bounds::repeatableBoundingNames[bounding_]);
    if (reader_.valid())
    {
        reader_->write(os);
    }
}


template<class Type>
Type Foam::interpolationTable<Type>::rateOfChange(const scalar value) const
{
    label n = this->size();

    if (n <= 1)
    {
        // There are not enough entries to provide a rate of change
        return 0;
    }

    const scalar minLimit = this->first().first();
    const scalar maxLimit = this->last().first();
    scalar lookupValue = value;

    if (lookupValue < minLimit)
    {
        switch (bounding_)
        {
            case bounds::repeatableBounding::ERROR:
            {
                FatalErrorInFunction
                    << "value (" << lookupValue << ") underflow" << nl
                    << exit(FatalError);
                break;
            }
            case bounds::repeatableBounding::WARN:
            {
                WarningInFunction
                    << "value (" << lookupValue << ") underflow" << nl
                    << "    Zero rate of change."
                    << endl;

                // Behaviour as per CLAMP
                return 0;
                break;
            }
            case bounds::repeatableBounding::CLAMP:
            {
                return 0;
                break;
            }
            case bounds::repeatableBounding::REPEAT:
            {
                // Adjust lookupValue to >= minLimit
                scalar span = maxLimit-minLimit;
                lookupValue = fmod(lookupValue-minLimit, span) + minLimit;
                break;
            }
        }
    }
    else if (lookupValue >= maxLimit)
    {
        switch (bounding_)
        {
            case bounds::repeatableBounding::ERROR:
            {
                FatalErrorInFunction
                    << "value (" << lookupValue << ") overflow" << nl
                    << exit(FatalError);
                break;
            }
            case bounds::repeatableBounding::WARN:
            {
                WarningInFunction
                    << "value (" << lookupValue << ") overflow" << nl
                    << "    Zero rate of change."
                    << endl;

                // Behaviour as per CLAMP
                return 0;
                break;
            }
            case bounds::repeatableBounding::CLAMP:
            {
                return 0;
                break;
            }
            case bounds::repeatableBounding::REPEAT:
            {
                // Adjust lookupValue <= maxLimit
                scalar span = maxLimit-minLimit;
                lookupValue = fmod(lookupValue-minLimit, span) + minLimit;
                break;
            }
        }
    }

    label lo = 0;
    label hi = 0;

    // look for the correct range
    for (label i = 0; i < n; ++i)
    {
        if (lookupValue >= List<Tuple2<scalar, Type>>::operator[](i).first())
        {
            lo = hi = i;
        }
        else
        {
            hi = i;
            break;
        }
    }

    if (lo == hi)
    {
        // we are at the end of the table - or there is only a single entry
        return 0;
    }
    else if (hi == 0)
    {
        // this treatment should only occur under these conditions:
        //  -> the 'REPEAT' treatment
        //  -> (0 <= value <= minLimit)
        //  -> minLimit > 0
        // Use the value at maxLimit as the value for value=0
        lo = n - 1;

        return
        (
            (
                List<Tuple2<scalar, Type>>::operator[](hi).second()
              - List<Tuple2<scalar, Type>>::operator[](lo).second()
            )
           /(
               List<Tuple2<scalar, Type>>::operator[](hi).first()
             + minLimit
             - List<Tuple2<scalar, Type>>::operator[](lo).first()
            )
        );
    }
    else
    {
        // normal rate of change
        return
        (
            (
                List<Tuple2<scalar, Type>>::operator[](hi).second()
              - List<Tuple2<scalar, Type>>::operator[](lo).second()
            )
           /(
                List<Tuple2<scalar, Type>>::operator[](hi).first()
              - List<Tuple2<scalar, Type>>::operator[](lo).first()
            )
        );
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
const Foam::Tuple2<Foam::scalar, Type>&
Foam::interpolationTable<Type>::operator[](const label i) const
{
    label ii = i;
    label n  = this->size();

    if (n <= 1)
    {
        ii = 0;
    }
    else if (ii < 0)
    {
        switch (bounding_)
        {
            case bounds::repeatableBounding::ERROR:
            {
                FatalErrorInFunction
                    << "index (" << ii << ") underflow" << nl
                    << exit(FatalError);
                break;
            }
            case bounds::repeatableBounding::WARN:
            {
                WarningInFunction
                    << "index (" << ii << ") underflow" << nl
                    << "    Continuing with the first entry"
                    << endl;

                // Behaviour as per 'CLAMP'
                ii = 0;
                break;
            }
            case bounds::repeatableBounding::CLAMP:
            {
                ii = 0;
                break;
            }
            case bounds::repeatableBounding::REPEAT:
            {
                while (ii < 0)
                {
                    ii += n;
                }
                break;
            }
        }
    }
    else if (ii >= n)
    {
        switch (bounding_)
        {
            case bounds::repeatableBounding::ERROR:
            {
                FatalErrorInFunction
                    << "index (" << ii << ") overflow" << nl
                    << exit(FatalError);
                break;
            }
            case bounds::repeatableBounding::WARN:
            {
                WarningInFunction
                    << "index (" << ii << ") overflow" << nl
                    << "    Continuing with the last entry"
                    << endl;

                // Behaviour as per 'CLAMP'
                ii = n - 1;
                break;
            }
            case bounds::repeatableBounding::CLAMP:
            {
                ii = n - 1;
                break;
            }
            case bounds::repeatableBounding::REPEAT:
            {
                while (ii >= n)
                {
                    ii -= n;
                }
                break;
            }
        }
    }

    return List<Tuple2<scalar, Type>>::operator[](ii);
}


template<class Type>
Type Foam::interpolationTable<Type>::operator()(const scalar value) const
{
    label n = this->size();

    if (n <= 1)
    {
        return this->first().second();
    }

    const scalar minLimit = this->first().first();
    const scalar maxLimit = this->last().first();
    scalar lookupValue = value;

    if (lookupValue < minLimit)
    {
        switch (bounding_)
        {
            case bounds::repeatableBounding::ERROR:
            {
                FatalErrorInFunction
                    << "value (" << lookupValue << ") underflow" << nl
                    << exit(FatalError);
                break;
            }
            case bounds::repeatableBounding::WARN:
            {
                WarningInFunction
                    << "value (" << lookupValue << ") underflow" << nl
                    << "    Continuing with the first entry"
                    << endl;

                // Behaviour as per CLAMP
                return this->first().second();
                break;
            }
            case bounds::repeatableBounding::CLAMP:
            {
                return this->first().second();
                break;
            }
            case bounds::repeatableBounding::REPEAT:
            {
                // adjust lookupValue to >= minLimit
                const scalar span = maxLimit-minLimit;
                lookupValue = fmod(lookupValue-minLimit, span) + minLimit;
                break;
            }
        }
    }
    else if (lookupValue >= maxLimit)
    {
        switch (bounding_)
        {
            case bounds::repeatableBounding::ERROR:
            {
                FatalErrorInFunction
                    << "value (" << lookupValue << ") overflow" << nl
                    << exit(FatalError);
                break;
            }
            case bounds::repeatableBounding::WARN:
            {
                WarningInFunction
                    << "value (" << lookupValue << ") overflow" << nl
                    << "    Continuing with the last entry"
                    << endl;

                // Behaviour as per 'CLAMP'
                return this->last().second();
                break;
            }
            case bounds::repeatableBounding::CLAMP:
            {
                return this->last().second();
                break;
            }
            case bounds::repeatableBounding::REPEAT:
            {
                // adjust lookupValue <= maxLimit
                const scalar span = maxLimit-minLimit;
                lookupValue = fmod(lookupValue-minLimit, span) + minLimit;
                break;
            }
        }
    }

    label lo = 0;
    label hi = 0;

    // look for the correct range
    for (label i = 0; i < n; ++i)
    {
        if (lookupValue >= List<Tuple2<scalar, Type>>::operator[](i).first())
        {
            lo = hi = i;
        }
        else
        {
            hi = i;
            break;
        }
    }

    if (lo == hi)
    {
        // we are at the end of the table - or there is only a single entry
        return List<Tuple2<scalar, Type>>::operator[](hi).second();
    }
    else if (hi == 0)
    {
        // this treatment should only occur under these conditions:
        //  -> the 'REPEAT' treatment
        //  -> (0 <= value <= minLimit)
        //  -> minLimit > 0
        // Use the value at maxLimit as the value for value=0
        lo = n - 1;

        return
        (
            List<Tuple2<scalar, Type>>::operator[](lo).second()
          + (
                List<Tuple2<scalar, Type>>::operator[](hi).second()
              - List<Tuple2<scalar, Type>>::operator[](lo).second()
            )
           *(lookupValue / minLimit)
        );
    }
    else
    {
        // normal interpolation
        return
        (
            List<Tuple2<scalar, Type>>::operator[](lo).second()
          + (
                List<Tuple2<scalar, Type>>::operator[](hi).second()
              - List<Tuple2<scalar, Type>>::operator[](lo).second()
            )
           *(
                lookupValue
              - List<Tuple2<scalar, Type>>::operator[](lo).first()
            )
           /(
                List<Tuple2<scalar, Type>>::operator[](hi).first()
              - List<Tuple2<scalar, Type>>::operator[](lo).first()
            )
        );
    }
}


// ************************************************************************* //
