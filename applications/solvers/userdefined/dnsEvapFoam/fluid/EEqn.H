{
    volScalarField& he = thermo.he();

    volScalarField Yk(0.0*Y[0]);
    
    // reset enthalpy diffusion term
    enthDiff *= 0.0;

    // enthalpy diffusion term
    forAll(Y, k)
    {
        if (composition.active(k))
        {
            volScalarField& Yk = Y[k];
            volScalarField enthK(0.0*he);
            forAll(mesh.C(), celli)
            {
                scalar enthKcell = composition.Hs(k, p[celli], thermo.T()[celli]);
                enthK[celli] += enthKcell;
            }
            
            enthDiff += turbulence->muEff()*invrhoInf/Sc*fvc::grad(Yk)*enthK;
        }
    }

    // sensible enthalpy equation
    fvScalarMatrix EEqn
    (
        fvm::ddt(he) 
      + mvConvection->fvmDiv(phi, he)
      - fvm::laplacian(turbulence->alphaEff()*invrhoInf, he)
     ==
        invrhoInf*parcels.Sh(he)
      - fvc::div(enthDiff)
      + fvOptions(he)
    );

    EEqn.relax();

    fvOptions.constrain(EEqn);

    EEqn.solve();

    fvOptions.correct(he);

    thermo.correct();
}
