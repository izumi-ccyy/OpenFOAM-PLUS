/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2016-2017 Wikki Ltd
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
    Foam::zeroGradientFaPatchField

Description

Author
    Zeljko Tukovic, FMENA
    Hrvoje Jasak, Wikki Ltd.

SourceFiles
    zeroGradientFaPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef zeroGradientFaPatchField_H
#define zeroGradientFaPatchField_H

#include "faPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class zeroGradientFaPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class zeroGradientFaPatchField
:
    public faPatchField<Type>
{

public:

    //- Runtime type information
    TypeName("zeroGradient");


    // Constructors

        //- Construct from patch and internal field
        zeroGradientFaPatchField
        (
            const faPatch&,
            const DimensionedField<Type, areaMesh>&
        );

        //- Construct from patch, internal field and dictionary
        zeroGradientFaPatchField
        (
            const faPatch&,
            const DimensionedField<Type, areaMesh>&,
            const dictionary&
        );

        //- Construct by mapping the given zeroGradientFaPatchField<Type>
        //  onto a new patch
        zeroGradientFaPatchField
        (
            const zeroGradientFaPatchField<Type>&,
            const faPatch&,
            const DimensionedField<Type, areaMesh>&,
            const faPatchFieldMapper&
        );

        //- Construct as copy
        zeroGradientFaPatchField
        (
            const zeroGradientFaPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<faPatchField<Type>> clone() const
        {
            return tmp<faPatchField<Type>>
            (
                new zeroGradientFaPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        zeroGradientFaPatchField
        (
            const zeroGradientFaPatchField<Type>&,
            const DimensionedField<Type, areaMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<faPatchField<Type>> clone
        (
            const DimensionedField<Type, areaMesh>& iF
        ) const
        {
            return tmp<faPatchField<Type>>
            (
                new zeroGradientFaPatchField<Type>(*this, iF)
            );
        }


    //- Destructor
    virtual ~zeroGradientFaPatchField<Type>()
    {}


    // Member functions

        // Evaluation functions

            //- Return gradient at boundary
            virtual tmp<Field<Type>> snGrad() const
            {
                return tmp<Field<Type>>
                (
                    new Field<Type>(this->size(), Zero)
                );
            }

            //- Evaluate the patch field
            virtual void evaluate
            (
                const Pstream::commsTypes commsType = Pstream::commsTypes::blocking
            );

            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<Field<Type>> valueInternalCoeffs
            (
                const tmp<scalarField>&
            ) const;

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<Field<Type>> valueBoundaryCoeffs
            (
                const tmp<scalarField>&
            ) const;

            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<Field<Type>> gradientInternalCoeffs() const;

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<Field<Type>> gradientBoundaryCoeffs() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "zeroGradientFaPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
