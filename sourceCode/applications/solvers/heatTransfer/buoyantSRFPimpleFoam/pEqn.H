//Changes to code to use Urel instead of U for solving in the rotating frame of reference

dimensionedScalar compressibility = fvc::domainIntegrate(psi);
bool compressible = (compressibility.value() > SMALL);

rho = thermo.rho();

// Thermodynamic density needs to be updated by psi*d(p) after the
// pressure solution
const volScalarField psip0(psi*p);

volScalarField rAUrel(1.0/UrelEqn.A());
surfaceScalarField rhorAUrelf("rhorAUrelf", fvc::interpolate(rho*rAUrel));
volVectorField HbyA(constrainHbyA(rAUrel*UrelEqn.H(), Urel, p_rgh));

surfaceScalarField phig(-rhorAUrelf*ghf*fvc::snGrad(rho)*mesh.magSf());

surfaceScalarField phiHbyA
(
    "phiHbyA",
    (
        (fvc::interpolate(rho*HbyA) & mesh.Sf())
      + rhorAUrelf*fvc::ddtCorr(rho, Urel, phi)
      //+ MRF.zeroFilter(rhorAUf*fvc::ddtCorr(rho, U, phi))
    )
  + phig
);

//fvOptions.makeRelative(fvc::interpolate(rho), phiHbyA);

// Update the pressure BCs to ensure flux consistency
constrainPressure(p_rgh, rho, Urel, phiHbyA, rhorAUrelf);

fvScalarMatrix p_rghDDtEqn
(
    fvc::ddt(rho) + psi*correction(fvm::ddt(p_rgh))
  + fvc::div(phiHbyA)
  ==
    fvOptions(psi, p_rgh, rho.name())
);

while (pimple.correctNonOrthogonal())
{
    fvScalarMatrix p_rghEqn
    (
        p_rghDDtEqn
      - fvm::laplacian(rhorAUrelf, p_rgh)
    );

    p_rghEqn.setReference
    (
        pRefCell,
        compressible ? getRefCellValue(p_rgh, pRefCell) : pRefValue
    );

    p_rghEqn.solve(mesh.solver(p_rgh.select(pimple.finalInnerIter())));

    if (pimple.finalNonOrthogonalIter())
    {
        // Calculate the conservative fluxes
        phi = phiHbyA + p_rghEqn.flux();

        // Explicitly relax pressure for momentum corrector
        p_rgh.relax();

        // Correct the momentum source with the pressure gradient flux
        // calculated from the relaxed pressure
        Urel = HbyA + rAUrel*fvc::reconstruct((phig + p_rghEqn.flux())/rhorAUrelf);
        Urel.correctBoundaryConditions();
        fvOptions.correct(Urel);
        K = 0.5*magSqr(Urel);
    }
}

p = p_rgh + rho*gh;

pressureControl.limit(p);

if (p_rgh.needReference())
{
    if (!compressible)
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
    }
    else
    {
        p += (initialMass - fvc::domainIntegrate(psi*p))
            /compressibility;
        thermo.correctRho(psi*p - psip0, rhoMin, rhoMax);
        rho = thermo.rho();
        p_rgh = p - rho*gh;
        p_rgh.correctBoundaryConditions();
    }
}
else
{
    thermo.correctRho(psi*p - psip0, rhoMin, rhoMax);
}

// Solve density continuity
#include "rhoEqn.H"
// Update continuity errors
#include "compressibleContinuityErrs.H"

rho = thermo.rho();

// Update pressure time derivative if needed
if (thermo.dpdt())
{
    dpdt = fvc::ddt(p);
}
