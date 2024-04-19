    // Solve the Momentum equation

//    MRF.correctBoundaryVelocity(Urel);

    tmp<fvVectorMatrix> tUrelEqn
    (
        fvm::ddt(rho, Urel) 
      + fvm::div(phi, Urel)
      + (SRF.Su())*rho
//      + MRF.DDt(rho, Urel)
      + turbulence.divDevRhoReff(Urel)
     ==
        fvOptions(rho, Urel)
    );
    fvVectorMatrix& UrelEqn = tUrelEqn.ref();

    UrelEqn.relax();

    fvOptions.constrain(UrelEqn);

    if (momentumPredictor)
    {
        solve
        (
            UrelEqn
          ==
            fvc::reconstruct
            (
                (
                  - ghf*fvc::snGrad(rho)
                  - fvc::snGrad(p_rgh)
                )*mesh.magSf()
            ),
            mesh.solver(Urel.select(finalIter))
        );

        fvOptions.correct(Urel);
        K = 0.5*magSqr(Urel);
    }
