{
    volScalarField& he = thermo.he();

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
      - fvm::laplacian(turbulence->alphaEff(), he)
     ==
        //rho*(U&g)
	rho*r()*magSqr(SRF->omega())*Urad()
      + radiation->Sh(thermo, he)
      + fvOptions(rho, he)
    );

    EEqn.relax();

    fvOptions.constrain(EEqn);

    EEqn.solve();

    fvOptions.correct(he);

    thermo.correct();
    radiation->correct();

    // Output temperature information
    Info<< "Min/max T:" << min(thermo.T()).value() << ' ' 
	<< max(thermo.T()).value() << endl;
    Info<< "Min/max h:" << min(thermo.he()).value() << ' '
	<< max(thermo.he()).value() << endl;
    Info<< "Min/max K:" << min(K()).value() << ' '
	<< max(K()).value() << endl;
    Info<< "Omega: " << SRF->omega() << endl;
    Info<< "magSq: " << magSqr(SRF->omega()) << endl;
}
