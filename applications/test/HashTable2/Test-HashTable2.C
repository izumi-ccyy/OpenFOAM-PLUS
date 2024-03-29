/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2011, 2019 OpenCFD Ltd.
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

Description
    Miscellaneous tests for HashTable

\*---------------------------------------------------------------------------*/

#include "HashTable.H"
#include "HashPtrTable.H"
#include "Map.H"
#include "FlatOutput.H"

using namespace Foam;


// Simple wrapper for testing purposes
class Label
{
    label data_;

public:

    Label()
    :
        data_(0)
    {}

    Label(label val)
    :
        data_(val)
    {}

    ~Label()
    {
        Info<<"delete label: " << data_ << endl;
    }

    // Some arbitrary non-const method (for testing)
    label increment()
    {
        return ++data_;
    }

    // Some arbitrary method (for testing)
    std::string info() const
    {
        return "boxed label=" + std::to_string(data_);
    }

    friend Ostream& operator<<(Ostream& os, const Label& val)
    {
        os  << val.data_;
        return os;
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    HashTable<label, Foam::string> table1
    {
        {"kjhk", 10},
        {"kjhk2", 12}
    };

    Info<< "table1: " << table1 << nl
        << "toc: " << flatOutput(table1.toc()) << nl;

    HashTable<label, label, Hash<label>> table2
    {
        {3, 10},
        {5, 12},
        {7, 16}
    };

    Info<< "table2: " << table2 << nl
        << "toc: " << flatOutput(table2.toc()) << nl;

    Map<label> table3(1);
    table3.transfer(table2);

    Info<< "table2: " << table2 << nl
        << "toc: " << flatOutput(table2.toc()) << nl;

    Info<< "table3: " << table3 << nl
        << "toc: " << flatOutput(table3.toc()) << nl;

    Map<label> table4(std::move(table3));

    Info<< "table3: " << table3 << nl
        << "toc: " << flatOutput(table3.toc()) << nl;

    Info<< "table4: " << table4 << nl
        << "toc: " << flatOutput(table4.toc()) << nl;

    {
        HashTable<Label, Foam::string> table1(0);
        table1.insert("abc", 5);
        table1.insert("def", 10);
        table1.insert("ghi", 15);
        table1.insert("jkl", 20);

        Info<< nl << "Table toc: " << flatOutput(table1.toc()) << nl;

        for (const word k : { "abc" })
        {
            const auto iter = table1.cfind(k);

            if (iter.good())
            {
                Info<< "have " << k << nl
                    << "    info: " << (*iter).info() << nl
                    // Good: does not compile
                    // << "    info: " << iter->info() << nl
                    ;
            }

            auto iter2 = table1.find(k);

            if (iter2.good())
            {
                Info<< "have " << k << nl
                    << "    incr: " << (*iter2).increment() << nl
                    // Good: does not compile
                    // << "    incr: " << iter2->increment() << nl
                    ;
            }
        }
    }

    {
        HashPtrTable<Label> ptable1(0);
        ptable1.insert("abc", autoPtr<Label>::New(5));
        ptable1.insert("def", autoPtr<Label>::New(10));
        ptable1.insert("ghi", autoPtr<Label>::New(15));
        ptable1.insert("jkl", autoPtr<Label>::New(20));

        Info<< nl << "PtrTable toc: " << flatOutput(ptable1.toc()) << nl;

        for (const word k : { "abc" })
        {
            const auto iter = ptable1.cfind(k);

            // Note: increment() changes contents of pointers,
            // not the pointers themselves.
            if (iter.good())
            {
                Info<< "have " << k << nl
                    << "    addr: " << uintptr_t(*iter) << nl
                    << "    info: " << (*iter)->info() << nl
                    << "    info: " << iter->info() << nl
                    << "    incr: " << iter->increment() << nl
                    ;
            }

            auto iter2 = ptable1.find(k);

            if (iter2.good())
            {
                Info<< "have " << k << nl
                    << "    incr: " << (*iter2)->increment() << nl
                    << "    incr: " << iter2->increment() << nl
                    ;
            }
        }


        // Attempt the same again

        HashTable<Label*> tableView1;
        HashTable<const Label*> tableView2;

        forAllConstIters(ptable1, iter)
        {
            tableView1.insert(iter.key(), iter.val());
            tableView2.insert(iter.key(), iter.val());
        }

        Info<< nl << "Table<pointer> toc: "
            << flatOutput(tableView1.toc()) << nl;

        for (const word k : { "abc" })
        {
            const auto iter1 = tableView1.cfind(k);

            // Note that increment changes contents of the pointers
            // not the table
            if (iter1.good())
            {
                Info<< "have " << k << nl
                    << "    addr: " << uintptr_t(*iter1) << nl
                    << "    info: " << (*iter1)->info() << nl
                    << "    info: " << iter1->info() << nl
                    << "    incr: " << iter1->increment() << nl
                    ;
            }

            auto iter2 = tableView2.cfind(k);
            if (iter2.good())
            {
                Info<< "have " << k << nl
                    << "    addr: " << uintptr_t(*iter2) << nl
                    << "    info: " << (*iter2)->info() << nl
                    << "    info: " << iter2->info() << nl
                    // Good: does not compile
                    // << "    incr: " << iter2->increment() << nl
                    ;
            }

            auto iter3 = tableView2.find(k);
            if (iter3.good())
            {
                Info<< "have " << k << nl
                    << "    addr: " << uintptr_t(*iter3) << nl
                    << "    info: " << (*iter3)->info() << nl
                    << "    info: " << iter3->info() << nl
                    // Good: does not compile
                    // << "    incr: " << iter3->increment() << nl
                    ;
            }
        }

        Info<< nl << "Ending scope" << nl;
    }

    {
        Info<< nl << "Table<labelList> copy/move/emplace insertion" << nl;

        HashTable<labelList> ltable1(0);
        ltable1.insert("abc", identity(2));
        ltable1.insert("def", identity(3));
        ltable1.insert("ghi", identity(4));
        ltable1.emplace("jkl", 10, -35);
        ltable1.emplace("mno");

        labelList list1(identity(4, -4));

        Info<<"move insert " << list1 << nl;

        ltable1.insert("pqr", std::move(list1));

        Info<<"after insert " << list1 << nl;

        Info<< nl << "HashTable<labelList>: "
            << ltable1 << nl;
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
