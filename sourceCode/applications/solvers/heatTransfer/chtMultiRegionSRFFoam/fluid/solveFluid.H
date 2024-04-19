if (finalIter)
{
    mesh.data::add("finalIteration", true);
}

if (frozenFlow)
{
    #include "EEqn.H"
}
else
{
    if (oCorr == 0)
    {
        #include "rhoEqn.H"
    }

    #include "UrelEqn.H"

    //Update absolute velcoity for calculating Urad
    U = Urel + SRF.U();
    //Update Urad for EEqn.H
    Urad = (U.component(2)*cos(phiAngle)+(U.component(1)*sin(phiAngle)));

    #include "YEqn.H"
    #include "EEqn.H"

    // --- PISO loop
    for (int corr=0; corr<nCorr; corr++)
    {
        #include "pEqn.H"
    }

    //Update pEqn corrected absolute velocity for 
    //calculating U and Urad
    U = Urel + SRF.U();
    Urad = (U.component(2)*cos(phiAngle)+(U.component(1)*sin(phiAngle)));
    Utan = (U.component(2)*sin(phiAngle)-(U.component(1)*cos(phiAngle)));

    turbulence.correct();

    rho = thermo.rho();
}

if (finalIter)
{
    mesh.data::remove("finalIteration");
}
