Info<< "Reading velocity field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Initialise the velocity internal field to zero
U = dimensionedVector(U.dimensions(), Zero);

surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fvc::flux(U)
);

if (args.found("initialiseUBCs"))
{
    U.correctBoundaryConditions();
    phi = fvc::flux(U);
}


// Construct a pressure field
// If it is available read it otherwise construct from the velocity BCs
// converting fixed-value BCs to zero-gradient and vice versa.

// Allow override from command-line -pName option
const word pName = args.opt<word>("pName", "p");

// Infer the pressure BCs from the velocity
wordList pBCTypes
(
    U.boundaryField().size(),
    fixedValueFvPatchScalarField::typeName
);

forAll(U.boundaryField(), patchi)
{
    if (U.boundaryField()[patchi].fixesValue())
    {
        pBCTypes[patchi] = zeroGradientFvPatchScalarField::typeName;
    }
}

// Note that registerObject is false for the pressure field. The pressure
// field in this solver doesn't have a physical value during the solution.
// It shouldn't be looked up and used by sub models or boundary conditions.
Info<< "Constructing pressure field " << pName << nl << endl;
volScalarField p
(
    IOobject
    (
        pName,
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE,
        false
    ),
    mesh,
    dimensionedScalar(sqr(dimVelocity), Zero),
    pBCTypes
);

// Infer the velocity potential BCs from the pressure
wordList PhiBCTypes
(
    p.boundaryField().size(),
    zeroGradientFvPatchScalarField::typeName
);

forAll(p.boundaryField(), patchi)
{
    if (p.boundaryField()[patchi].fixesValue())
    {
        PhiBCTypes[patchi] = fixedValueFvPatchScalarField::typeName;
    }
}

Info<< "Constructing velocity potential field Phi\n" << endl;
volScalarField Phi
(
    IOobject
    (
        "Phi",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(dimLength*dimVelocity, Zero),
    PhiBCTypes
);

label PhiRefCell = 0;
scalar PhiRefValue = 0;
setRefCell
(
    Phi,
    potentialFlow.dict(),
    PhiRefCell,
    PhiRefValue
);
mesh.setFluxRequired(Phi.name());

#include "createMRF.H"

// Add solver-specific interpolations
{
    wordHashSet& nonInt =
        const_cast<wordHashSet&>(Stencil::New(mesh).nonInterpolatedFields());

    nonInt.insert("cellMask");
    nonInt.insert("interpolatedCells");

}

// Mask field for zeroing out contributions on hole cells
#include "createCellMask.H"

// Create bool field with interpolated cells
#include "createInterpolatedCells.H"
