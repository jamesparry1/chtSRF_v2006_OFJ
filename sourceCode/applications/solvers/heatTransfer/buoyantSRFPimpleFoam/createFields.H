Info<< "Reading thermophysical properties\n" << endl;

autoPtr<rhoThermo> pThermo(rhoThermo::New(mesh));
rhoThermo& thermo = pThermo();
thermo.validate(args.executable(), "h", "e");

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    thermo.rho()
);

volScalarField& p = thermo.p();

Info<< "Reading field Urel\n" << endl;
volVectorField Urel
(
    IOobject
    (
        "Urel",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Creating SRF model\n" << endl;
autoPtr<SRF::SRFModel> SRF
(
      SRF::SRFModel::New(Urel)
);

// Create the absolute velocity
Info<< "Creating field U\n" << endl;
volVectorField U
(
     IOobject
     (
	"U",
        runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::AUTO_WRITE
     ),
     Urel + SRF->U()
);

#include "compressibleCreatePhiSRF.H"

pressureControl pressureControl(p, rho, pimple.dict(), false);

Info<< "Creating turbulence model\n" << endl;
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New
    (
        rho,
        Urel,
        phi,
        thermo
    )
);


#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"


Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
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

// Force p_rgh to be consistent with p
p_rgh = p - rho*gh;

mesh.setFluxRequired(p_rgh.name());

label pRefCell = 0;
scalar pRefValue = 0.0;

if (p_rgh.needReference())
{
    setRefCell
    (
        p,
        p_rgh,
        pimple.dict(),
        pRefCell,
        pRefValue
    );

    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
}

dimensionedScalar initialMass("initialMass", fvc::domainIntegrate(rho));

#include "createDpdt.H"

#include "createK.H"

//#include "createMRF.H"
#include "createRadiationModel.H"
#include "createFvOptions.H"


Info<< "Creating field r\n" << endl;
volScalarField r
(
    IOobject
    (
	"r",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::AUTO_WRITE
    ),
    sqrt(magSqr(mesh.C().component(1))+magSqr(mesh.C().component(2)))
);

Info<< "Creating field phiAngle\n" << endl;
volScalarField phiAngle
(
    IOobject
    (
	"phiAngle",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::AUTO_WRITE
    ),
    Foam::atan2(mesh.C().component(1),mesh.C().component(2))
);

Info<< "Creating field Urad\n" << endl;
volScalarField Urad
(
    IOobject
    (
	"Urad",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::AUTO_WRITE
    ),
    (U.component(2)*cos(phiAngle)+U.component(1)*sin(phiAngle))
);

Info<< "Creating field Utan\n" << endl;
volScalarField Utan
(
    IOobject
    (
	"Utan",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::AUTO_WRITE
    ),
    (U.component(2)*sin(phiAngle)-U.component(1)*cos(phiAngle))
);

const dimensionedScalar rhoMax("rhoMax", dimDensity, GREAT, pimple.dict());
const dimensionedScalar rhoMin("rhoMin", dimDensity, Zero, pimple.dict());
