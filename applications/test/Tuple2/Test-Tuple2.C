/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2010, 2019 OpenCFD Ltd.
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

Application
    Test-Tuple2

Description
    Test construction, comparision etc for Tuple2 and Pair.

\*---------------------------------------------------------------------------*/

#include "labelPair.H"
#include "Tuple2.H"
#include "label.H"
#include "scalar.H"
#include "List.H"
#include "ListOps.H"
#include "ops.H"
#include <functional>

using namespace Foam;

// Test for special comparison operation using compareOp
// Normal sort on label, reverse sort on scalar
struct special1
{
    typedef Tuple2<label, scalar> type;

    bool operator()(const type& a, const type& b) const
    {
        int val = compareOp<label>()(a.first(), b.first());
        return (val == 0) ? (b.second() < a.second()) : (val < 0);
    }
};


// Test for special comparison operation using compareOp
// Normal sort on scalar, reverse sort on label
struct special2
{
    typedef Tuple2<label, scalar> type;

    bool operator()(const type& a, const type& b) const
    {
        scalar val = compareOp<scalar>()(a.second(), b.second());
        return (val == 0) ? (b.first() < a.first()) : (val < 0);
    }
};


// Print info
void printTuple2(const Tuple2<word, word>& t)
{
    Info<< "tuple: " << t << nl;

    Info<< "first  @: " << uintptr_t(t.first().data()) << nl;
    Info<< "second @: " << uintptr_t(t.second().data()) << nl;
}


// Print info
void printTuple2(const Pair<word>& t)
{
    Info<< "tuple: " << t << nl;

    Info<< "first  @: " << uintptr_t(t.first().data()) << nl;
    Info<< "second @: " << uintptr_t(t.second().data()) << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main()
{
    typedef Tuple2<label, scalar> indexedScalar;

    Info<< "null constructed Tuple: " << indexedScalar() << nl;
    Info<< "null constructed Pair: "  << Pair<scalar>() << nl;

    indexedScalar t2(1, 3.2);

    Info<< "Foam::Tuple2: "
        << t2 << " => "
        << t2.first() << ' ' << t2.second() << nl;

    // As list. Generated so that we have duplicate indices
    List<indexedScalar> list1(3*4);
    for (label i = 0; i < 4; ++i)
    {
        const label j = (i+1);
        const label idx = ((i % 2) ? -1 : 1) * (j);

        list1[i]   = indexedScalar(idx, (j*j));
        list1[i+4] = indexedScalar(idx, 2*j);    // duplicate index
        list1[i+8] = indexedScalar(idx+12, 2*j); // duplicate value
    }

    Info<< "Unsorted tuples:" << nl << list1 << nl;

    Foam::sort(list1, std::less<indexedScalar>());

    Info<< "sorted tuples:" << nl << list1 << nl;

    Foam::sort(list1, std::greater<indexedScalar>());

    Info<< "reverse sorted tuples:" << nl << list1 << nl;

    Foam::sort(list1, special1());

    Info<< "special sorted tuples - sort on index, reverse on value:"
        << nl << list1 << nl;

    Foam::sort(list1, special2());

    Info<< "special sorted tuples - sort on value, reverse on index:"
        << nl << list1 << nl;


    {
        Info<< nl << nl << "Foam::Pair" << nl;

        typedef Pair<label> indexedLabel;

        indexedLabel pr(1, 3);

        Info<< "pair: "
            << pr << " => "
            << pr.first() << ' ' << pr.second() << nl;

        List<indexedLabel> list2 = ListOps::create<indexedLabel>
        (
            list1,
            [](const indexedScalar& t2)
            {
                return indexedLabel(t2.first(), t2.second());
            }
        );

        Info<< "Unsorted pairs:" << nl << list2 << nl;
    }


    {
        Info<< nl << nl << "std::pair" << nl;

        typedef std::pair<label, label> indexedLabel;

        indexedLabel pr(1, 3);

        Info<< "pair: "
            << pr << " => "
            << pr.first << ' ' << pr.second << nl;

        List<indexedLabel> list2 = ListOps::create<indexedLabel>
        (
            list1,
            [](const indexedScalar& t2)
            {
                return indexedLabel(t2.first(), t2.second());
            }
        );

        Info<< "Unsorted pairs:" << nl << list2 << nl;
    }


    {
        word word1("hello");
        word word2("word");

        Info<< "create with " << word1 << " @ " << uintptr_t(word1.data())
            << " " << word2 << " @ " << uintptr_t(word2.data()) << nl;

        Tuple2<word, word> tup(std::move(word2), std::move(word1));

        printTuple2(tup);

        Info<< "input is now " << word1 << " @ " << uintptr_t(word1.data())
            << " " << word2 << " @ " << uintptr_t(word2.data()) << nl;
    }

    {
        word word1("hello");
        word word2("word");

        Info<< "create with " << word1 << " @ " << uintptr_t(word1.data())
            << " " << word2 << " @ " << uintptr_t(word2.data()) << nl;

        Pair<word> tup(std::move(word2), std::move(word1));

        printTuple2(tup);

        Info<< "input is now " << word1 << " @ " << uintptr_t(word1.data())
            << " " << word2 << " @ " << uintptr_t(word2.data()) << nl;
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
