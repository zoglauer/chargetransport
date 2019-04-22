/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) YEAR OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "fixedValueFvPatchFieldTemplate.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "unitConversion.H"
//{{{ begin codeInclude

//}}} end codeInclude


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode

//}}} end localCode


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

extern "C"
{
    // dynamicCode:
    // SHA1 = 43bf39ff0df38d6ac1dcb9452712b670694f53ec
    //
    // unique function name that can be checked if the correct library version
    // has been loaded
    void frontBC_43bf39ff0df38d6ac1dcb9452712b670694f53ec(bool load)
    {
        if (load)
        {
            // code that can be explicitly executed after loading
        }
        else
        {
            // code that can be explicitly executed before unloading
        }
    }
}

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

makeRemovablePatchTypeField
(
    fvPatchScalarField,
    frontBCFixedValueFvPatchScalarField
);


const char* const frontBCFixedValueFvPatchScalarField::SHA1sum =
    "43bf39ff0df38d6ac1dcb9452712b670694f53ec";


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

frontBCFixedValueFvPatchScalarField::
frontBCFixedValueFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(p, iF)
{
    if (false)
    {
        Info<<"construct frontBC sha1: 43bf39ff0df38d6ac1dcb9452712b670694f53ec"
            " from patch/DimensionedField\n";
    }
}


frontBCFixedValueFvPatchScalarField::
frontBCFixedValueFvPatchScalarField
(
    const frontBCFixedValueFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<scalar>(ptf, p, iF, mapper)
{
    if (false)
    {
        Info<<"construct frontBC sha1: 43bf39ff0df38d6ac1dcb9452712b670694f53ec"
            " from patch/DimensionedField/mapper\n";
    }
}


frontBCFixedValueFvPatchScalarField::
frontBCFixedValueFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<scalar>(p, iF, dict)
{
    if (false)
    {
        Info<<"construct frontBC sha1: 43bf39ff0df38d6ac1dcb9452712b670694f53ec"
            " from patch/dictionary\n";
    }
}


frontBCFixedValueFvPatchScalarField::
frontBCFixedValueFvPatchScalarField
(
    const frontBCFixedValueFvPatchScalarField& ptf
)
:
    fixedValueFvPatchField<scalar>(ptf)
{
    if (false)
    {
        Info<<"construct frontBC sha1: 43bf39ff0df38d6ac1dcb9452712b670694f53ec"
            " as copy\n";
    }
}


frontBCFixedValueFvPatchScalarField::
frontBCFixedValueFvPatchScalarField
(
    const frontBCFixedValueFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(ptf, iF)
{
    if (false)
    {
        Info<<"construct frontBC sha1: 43bf39ff0df38d6ac1dcb9452712b670694f53ec "
            "as copy/DimensionedField\n";
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

frontBCFixedValueFvPatchScalarField::
~frontBCFixedValueFvPatchScalarField()
{
    if (false)
    {
        Info<<"destroy frontBC sha1: 43bf39ff0df38d6ac1dcb9452712b670694f53ec\n";
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void frontBCFixedValueFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (false)
    {
        Info<<"updateCoeffs frontBC sha1: 43bf39ff0df38d6ac1dcb9452712b670694f53ec\n";
    }

//{{{ begin code
    #line 44 "/home/jgarg20/chargetransport/CorysCode/planDet/0/phi.boundaryField.front"
const fvPatch& boundaryPatch = patch();
             const vectorField& Cf = boundaryPatch.Cf();

             scalarField& field = *this;
             field = patchInternalField();

             //scalar min = 0.00325;
             //scalar max = 0.005;
             scalar xmin = 0.00325;
             scalar xmax = 0.005;
             scalar ymin = 0.00325;
             scalar ymax = 0.005;
             /*
             scalar xmin = 0.00325;
             scalar xmax = 0.005;
             scalar ymin = 0.001;
             scalar ymax = 0.002;
             /**/

             forAll(Cf, faceI)
             {
                if (
                     (Cf[faceI].x() > xmin) &&
                     (Cf[faceI].x() < xmax) &&
                     (Cf[faceI].y() > ymin) &&
                     (Cf[faceI].y() < ymax)
                   )
                 {
                     field[faceI] = 1;
                 }
             }
//}}} end code

    this->fixedValueFvPatchField<scalar>::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

