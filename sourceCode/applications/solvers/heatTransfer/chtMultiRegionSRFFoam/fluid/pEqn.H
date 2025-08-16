//Changes to code to use Urel instead of U for solving in the rotating frame of reference

bool closedVolume = p_rgh.needReference();
dimensionedScalar compressibility = fvc::domainIntegrate(psi);
bool compressible = (compressibility.value() > SMALL);

rho = thermo.rho();

// Thermodynamic density needs to be updated by psi*d(p) after the
// pressure solution
const volScalarField psip0(psi*p);

volScalarField rAUrel("rAUrel", 1.0/UrelEqn.A());
surfaceScalarField rhorAUrelf("rhorAUrelf", fvc::interpolate(rho*rAUrel));
volVectorField HbyA(constrainHbyA(rAUrel*UrelEqn.H(), Urel, p_rgh));

surfaceScalarField phig(-rhorAUrelf*ghf*fvc::snGrad(rho)*mesh.magSf());

surfaceScalarField phiHbyA
(
    "phiHbyA",
    (
	(fvc::interpolate(rho*HbyA) & mesh.Sf())
      + rhorAUrelf*fvc::ddtCorr(rho, Urel, phi)
//        fvc::flux(rho*HbyA)
//      + MRF.zeroFilter(rhorAUf*fvc::ddtCorr(rho, U, phi))
    )
  + phig
);

//fvOptions.makeRelative(fvc::interpolate(rho), phiHbyA);
//MRF.makeRelative(fvc::interpolate(rho), phiHbyA);

// Update the pressure BCs to ensure flux consistency
constrainPressure(p_rgh, rho, Urel, phiHbyA, rhorAUrelf, MRF);

{
    fvScalarMatrix p_rghDDtEqn
    (
        fvc::ddt(rho) + psi*correction(fvm::ddt(p_rgh))
      + fvc::div(phiHbyA)
    );

    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
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

        p_rghEqn.solve
        (
            mesh.solver
            (
                p_rgh.select
                (
                    (
                        oCorr == nOuterCorr-1
                     && corr == nCorr-1
                     && nonOrth == nNonOrthCorr
                    )
                )
            )
        );

        if (nonOrth == nNonOrthCorr)
        {
            phi = phiHbyA + p_rghEqn.flux();

            p_rgh.relax();

            Urel = HbyA
              + rAUrel*fvc::reconstruct((phig + p_rghEqn.flux())/rhorAUrelf);
            Urel.correctBoundaryConditions();
            fvOptions.correct(Urel);
            K = 0.5*magSqr(Urel);
        }
    }

    p = p_rgh + rho*gh;

}

pressureControl.limit(p);

// For closed-volume cases adjust the pressure and density levels
// to obey overall mass continuity
if (closedVolume)
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
    thermo.correctRho(psi*p - psip0,  rhoMin, rhoMax);
}

#include "rhoEqn.H"
#include "compressibleContinuityErrors.H"

rho = thermo.rho();

// Update pressure time derivative if needed
if (thermo.dpdt())
{
    dpdt = fvc::ddt(p);
}
