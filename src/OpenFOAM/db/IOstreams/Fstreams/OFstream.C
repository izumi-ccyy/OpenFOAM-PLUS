/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2010, 2017-2019 OpenCFD Ltd.
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

#include "OFstream.H"
#include "OSspecific.H"
#include "gzstream.h"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(OFstream, 0);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::Detail::OFstreamAllocator::OFstreamAllocator
(
    const fileName& pathname,
    IOstream::compressionType compression,
    const bool append
)
:
    allocatedPtr_(nullptr)
{
    if (pathname.empty())
    {
        if (OFstream::debug)
        {
            InfoInFunction << "Cannot open null file " << endl;
        }
    }

    std::ios_base::openmode mode(std::ios_base::out|std::ios_base::binary);
    if (append)
    {
        mode |= std::ios_base::app;
    }

    if (compression == IOstream::COMPRESSED)
    {
        // Get identically named uncompressed version out of the way
        fileName::Type pathType = Foam::type(pathname, false);
        if (pathType == fileName::FILE || pathType == fileName::LINK)
        {
            rm(pathname);
        }
        fileName gzPathName(pathname + ".gz");

        if (!append && Foam::type(gzPathName) == fileName::LINK)
        {
            // Disallow writing into softlink to avoid any problems with
            // e.g. softlinked initial fields
            rm(gzPathName);
        }

        allocatedPtr_ = new ogzstream(gzPathName.c_str(), mode);
    }
    else
    {
        // get identically named compressed version out of the way
        fileName gzPathName(pathname + ".gz");
        fileName::Type gzType = Foam::type(gzPathName, false);
        if (gzType == fileName::FILE || gzType == fileName::LINK)
        {
            rm(gzPathName);
        }
        if (!append && Foam::type(pathname, false) == fileName::LINK)
        {
            // Disallow writing into softlink to avoid any problems with
            // e.g. softlinked initial fields
            rm(pathname);
        }

        allocatedPtr_ = new std::ofstream(pathname, mode);
    }
}


Foam::Detail::OFstreamAllocator::~OFstreamAllocator()
{
    deallocate();
}


void Foam::Detail::OFstreamAllocator::deallocate()
{
    if (allocatedPtr_)
    {
        delete allocatedPtr_;
        allocatedPtr_ = nullptr;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::OFstream::OFstream
(
    const fileName& pathname,
    streamFormat format,
    versionNumber version,
    compressionType compression,
    const bool append
)
:
    Detail::OFstreamAllocator(pathname, compression, append),
    OSstream
    (
        *allocatedPtr_,
        pathname,
        format,
        version,
        compression
    )
{
    setClosed();
    setState(allocatedPtr_->rdstate());

    if (!good())
    {
        if (debug)
        {
            InfoInFunction
                << "Could not open file " << pathname
                << " for output" << nl << info() << Foam::endl;
        }

        setBad();
    }
    else
    {
        setOpened();
    }

    lineNumber_ = 1;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

std::ostream& Foam::OFstream::stdStream()
{
    if (!allocatedPtr_)
    {
        FatalErrorInFunction
            << "No stream allocated." << abort(FatalError);
    }
    return *allocatedPtr_;
}


const std::ostream& Foam::OFstream::stdStream() const
{
    if (!allocatedPtr_)
    {
        FatalErrorInFunction
            << "No stream allocated." << abort(FatalError);
    }
    return *allocatedPtr_;
}


void Foam::OFstream::print(Ostream& os) const
{
    os  << "OFstream: ";
    OSstream::print(os);
}


// ************************************************************************* //
