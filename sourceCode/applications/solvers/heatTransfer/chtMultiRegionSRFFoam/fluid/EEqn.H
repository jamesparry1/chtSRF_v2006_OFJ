{
    volScalarField& he = thermo.he();

    Info<< "Solving EEqn" << nl << endl;

    fvScalarMatrix EEqn
    (
        fvm::ddt(rho, he) + fvm::div(phi, he)
      + fvc::ddt(rho, K) + fvc::div(phi, K)
      + (
            he.name() == "e"
          ? fvc::div
            (
                fvc::absolute(phi/fvc::interpolate(rho), Urel),
                p,
                "div(phiv,p)"
            )
          : -dpdt
        )
      - fvm::laplacian(turbulence.alphaEff(), he)
      ==
      //  rho*(U&g)
        rho*r()*magSqr(SRF.omega())*Urad()
      + rad.Sh(thermo, he)
      + Qdot
      + fvOptions(rho, he)
    );

    EEqn.relax();

    fvOptions.constrain(EEqn);

    EEqn.solve(mesh.solver(he.select(finalIter)));

    fvOptions.correct(he);

    thermo.correct();
    rad.correct();

    Info<< "Min/max T:" << min(thermo.T()).value() << ' '
        << max(thermo.T()).value() << endl;
    Info<< "Min/max h:" << min(thermo.he()).value() << ' '
        << max(thermo.he()).value() << endl;
    Info<< "Min/max K:" << min(K()).value() << ' '
        << max(K()).value() << endl;
    Info<< "Omega: " << SRF.omega() << endl;
    Info<< "magSq: " << magSqr(SRF.omega()) << endl;
}
