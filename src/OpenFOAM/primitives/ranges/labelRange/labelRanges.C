/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2010-2010, 2017 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2011 OpenFOAM Foundation
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

#include "labelRanges.H"
#include "ListOps.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::labelRanges::insertBefore
(
    const label insert,
    const labelRange& range
)
{
    // insert via copying up
    label nElem = this->size();

    if (labelRange::debug)
    {
        Info<<"before insert "
            << nElem << " elements, insert at " << insert << nl
            << *this << endl;
    }

    StorageContainer::setSize(nElem+1);

    if (labelRange::debug)
    {
        Info<<"copy between " << nElem << " and " << insert << nl;
    }

    for (label i = nElem-1; i >= insert; --i)
    {
        if (labelRange::debug)
        {
            Info<<"copy from " << (i) << " to " << (i+1) << nl;
        }

        StorageContainer::operator[](i+1) = StorageContainer::operator[](i);
    }

    // finally insert the range
    if (labelRange::debug)
    {
        Info<< "finally insert the range at " << insert << nl;
    }
    StorageContainer::operator[](insert) = range;
}


void Foam::labelRanges::purgeEmpty()
{
    // purge empty ranges by copying down
    label nElem = 0;
    forAll(*this, elemI)
    {
        if (!StorageContainer::operator[](elemI).empty())
        {
            if (nElem != elemI)
            {
                StorageContainer::operator[](nElem) =
                    StorageContainer::operator[](elemI);
            }
            ++nElem;
        }
    }

    // truncate
    this->StorageContainer::setSize(nElem);
}


Foam::Ostream& Foam::labelRanges::printRange
(
    Ostream& os,
    const labelRange& range
) const
{
    if (range.empty())
    {
        os  << "empty";
    }
    else
    {
        os  << range << " = " << range.first() << ":" << range.last();
    }
    return os;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::labelRanges::labelRanges(Istream& is)
{
    is  >> *this;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::labelRanges::add(const labelRange& range)
{
    if (range.empty())
    {
        return false;
    }
    else if (this->empty())
    {
        this->append(range);
        return true;
    }

    // find the correct place for insertion
    forAll(*this, elemI)
    {
        labelRange& currRange = StorageContainer::operator[](elemI);

        if (currRange.overlaps(range, true))
        {
            // absorb into the existing (adjacent/overlapping) range
            currRange.join(range);

            // might connect with the next following range(s)
            for (; elemI < this->size()-1; ++elemI)
            {
                labelRange& nextRange = StorageContainer::operator[](elemI+1);
                if (currRange.overlaps(nextRange, true))
                {
                    currRange.join(nextRange);
                    nextRange.clear();
                }
                else
                {
                    break;
                }
            }

            // done - remove any empty ranges that might have been created
            purgeEmpty();
            return true;
            break;
        }
        else if (range < currRange)
        {
            insertBefore(elemI, range);
            return true;
            break;
        }
    }


    // not found: simply append
    this->append(range);

    return true;
}


bool Foam::labelRanges::remove(const labelRange& range)
{
    bool status = false;
    if (range.empty() || this->empty())
    {
        return status;
    }

    forAll(*this, elemI)
    {
        labelRange& currRange = StorageContainer::operator[](elemI);

        if (range.first() > currRange.first())
        {
            if (range.last() < currRange.last())
            {
                // removal of range fragments of currRange

                if (labelRange::debug)
                {
                    Info<<"Fragment removal ";
                    printRange(Info, range) << " from ";
                    printRange(Info, currRange) << endl;
                }

                // left-hand-side fragment: insert before current range
                label lower = currRange.first();
                label upper = range.first() - 1;

                labelRange fragment(lower, upper - lower + 1);

                // right-hand-side fragment
                lower = range.last() + 1;
                upper = currRange.last();

                currRange = labelRange(lower, upper - lower + 1);
                status = true;
                insertBefore(elemI, fragment);

                if (labelRange::debug)
                {
                    Info<<"fragment ";
                    printRange(Info, fragment) << endl;
                    Info<<"yields ";
                    printRange(Info, currRange) << endl;
                }

                // fragmentation can only affect a single range
                // thus we are done
                break;
            }
            else if (range.first() <= currRange.last())
            {
                // keep left-hand-side, remove right-hand-side

                if (labelRange::debug)
                {
                    Info<<"RHS removal ";
                    printRange(Info, range) << " from ";
                    printRange(Info, currRange) << endl;
                }

                const label lower = currRange.first();
                const label upper = range.first() - 1;

                currRange = labelRange(lower, upper - lower + 1);
                status = true;

                if (labelRange::debug)
                {
                    Info<<"yields ";
                    printRange(Info, currRange) << endl;
                }
            }
        }
        else if (range.first() <= currRange.first())
        {
            if (range.last() >= currRange.first())
            {
                // remove left-hand-side, keep right-hand-side

                if (labelRange::debug)
                {
                    Info<<"LHS removal ";
                    printRange(Info, range) << " from ";
                    printRange(Info, currRange) << endl;
                }

                const label lower = range.last() + 1;
                const label upper = currRange.last();

                currRange = labelRange(lower, upper - lower + 1);
                status = true;

                if (labelRange::debug)
                {
                    Info<<"yields ";
                    printRange(Info, currRange) << endl;
                }
            }
        }
    }

    purgeEmpty();

    return status;
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, labelRanges& ranges)
{
    is  >> static_cast<labelRanges::StorageContainer&>(ranges);
    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const labelRanges& ranges)
{
    os  << static_cast<const labelRanges::StorageContainer&>(ranges);
    return os;
}


// ************************************************************************* //
