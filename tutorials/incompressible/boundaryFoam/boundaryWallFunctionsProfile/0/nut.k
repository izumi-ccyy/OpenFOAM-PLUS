/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v1906                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    location    "0";
    object      nut.k;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 2 -1 0 0 0 0];

internalField   uniform 0;

boundaryField
{
    lowerWall
    {
        type            nutkWallFunction;
        value           $internalField;
    }
    upperWall
    {
        type            nutkWallFunction;
        value           $internalField;
    }
    front
    {
        type            cyclic;
    }
    back
    {
        type            cyclic;
    }
    defaultFaces
    {
        type            empty;
    }
}


// ************************************************************************* //
