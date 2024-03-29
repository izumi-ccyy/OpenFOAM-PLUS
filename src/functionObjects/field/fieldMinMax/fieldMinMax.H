/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2010, 2015-2017 OpenCFD Ltd.
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

Class
    Foam::functionObjects::fieldMinMax

Group
    grpFieldFunctionObjects

Description
    Calculates the value and location of scalar minimum and maximum for a list
    of user-specified fields.

    For variables with a rank greater than zero, either the min/max of a
    component value or the magnitude is reported.  When operating in parallel,
    the processor owning the value is also given.

Usage
    Example of function object specification:
    \verbatim
    fieldMinMax1
    {
        type        fieldMinMax;
        libs        ("libfieldFunctionObjects.so");
        ...
        writeToFile yes;
        log         yes;
        location    yes;
        mode        magnitude;
        fields      (U p);
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property     | Description             | Required    | Default value
        type         | type name: fieldMinMax  | yes         |
        writeToFile  | write min/max data to file | no       | yes
        log          | write min/max data to standard output | no | yes
        location     | write location of the min/max value | no | yes
        mode         | calculation mode: magnitude or component | no | magnitude
        fields       | list of fields to process | yes       |
    \endtable

    Output data is written to the file \<timeDir\>/fieldMinMax.dat

See also
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::writeFile

SourceFiles
    fieldMinMax.C
    fieldMinMaxTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_fieldMinMax_H
#define functionObjects_fieldMinMax_H

#include "Switch.H"
#include "Enum.H"
#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "vector.H"
#include "volFieldSelection.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                         Class fieldMinMax Declaration
\*---------------------------------------------------------------------------*/

class fieldMinMax
:
    public fvMeshFunctionObject,
    public writeFile
{
public:

    // Public enumerations

        enum modeType
        {
            mdMag,      //!< magnitude
            mdCmpt      //!< component
        };


protected:

    // Protected data

        //- Mode type names
        static const Enum<modeType> modeTypeNames_;

        //- Write location of min/max values?
        bool location_;

        //- Mode for min/max - only applicable for ranks > 0
        modeType mode_;

        //- Fields to assess min/max
        volFieldSelection fieldSet_;


    // Protected Member Functions

        //- Helper function to write the output
        template<class Type>
        void output
        (
            const word& fieldName,
            const word& outputName,
            const label minCell,
            const label maxCell,
            const vector& minC,
            const vector& maxC,
            const label minProci,
            const label maxProci,
            const Type& minValue,
            const Type& maxValue
        );


        //- Output file header information
        virtual void writeFileHeader(Ostream& os);

        //- No copy construct
        fieldMinMax(const fieldMinMax&) = delete;

        //- No copy assignment
        void operator=(const fieldMinMax&) = delete;

        //- Calculate the field min/max for a given field type
        template<class Type>
        void calcMinMaxFieldType
        (
            const GeometricField<Type, fvPatchField, volMesh>& field,
            const word& outputFieldName
        );

        //- Calculate the field min/max
        template<class Type>
        void calcMinMaxFields
        (
            const word& fieldName,
            const modeType& mode
        );


public:

    //- Runtime type information
    TypeName("fieldMinMax");


    // Constructors

        //- Construct from Time and dictionary
        fieldMinMax
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~fieldMinMax() = default;


    // Member Functions

        //- Read the field min/max data
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Write the fieldMinMax
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fieldMinMaxTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
