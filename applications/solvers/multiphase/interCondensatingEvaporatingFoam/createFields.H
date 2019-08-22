Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh // pressure
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U // velocity
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

#include "createPhi.H" // phi

// Creating e based thermo
autoPtr<twoPhaseMixtureEThermo> thermo // thermo, defined by twoPhaseMixtureEThermo
(
    new twoPhaseMixtureEThermo(U, phi)
);

// Create mixture and
Info<< "Creating temperaturePhaseChangeTwoPhaseMixture\n" << endl;
autoPtr<temperaturePhaseChangeTwoPhaseMixture> mixture =
    temperaturePhaseChangeTwoPhaseMixture::New(thermo(), mesh); // mixture, defined by left


volScalarField& alpha1(thermo->alpha1()); // alpha is a member of thermo
volScalarField& alpha2(thermo->alpha2());

const dimensionedScalar& rho1 = thermo->rho1();
const dimensionedScalar& rho2 = thermo->rho2();

// Need to store rho for ddt(rho, U)
volScalarField rho // rho calculated by alpha and rho for every phase
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    alpha1*rho1 + alpha2*rho2
);
rho.oldTime();

// Construct interface from alpha1 distribution
interfaceProperties interface // interface
(
    alpha1,
    U,
    thermo->transportPropertiesDict()
);

// Construct incompressible turbulence model
autoPtr<incompressible::turbulenceModel> turbulence // incompressible turbulence model
(
    incompressible::turbulenceModel::New(U, phi, thermo())
);

#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"

volScalarField& p = thermo->p(); // define pressure from thermo

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
    p_rgh = p - rho*gh;
}

mesh.setFluxRequired(p_rgh.name());
mesh.setFluxRequired(alpha1.name());

// Turbulent Prandtl number
dimensionedScalar Prt("Prt", dimless, thermo->transportPropertiesDict()); // define Prt from thermo

volScalarField kappaEff // define kappaEff by thermo
(
    IOobject
    (
        "kappaEff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    thermo->kappa()
);

// Need to store rho for ddt(rhoCp, U)
volScalarField rhoCp // define rho * cp by thermo
(
    IOobject
    (
        "rhoCp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho*thermo->Cp()
);

rhoCp.oldTime();
