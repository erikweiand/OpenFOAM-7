tmp<fv::convectionScheme<scalar>> mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.divScheme("div(phi,Yk_h)")
    )
);


{
    volScalarField Yt(0.0*Y[0]);

    forAll(Y, k)
    {
        if (k != inertIndex && composition.active(k))
        {
            volScalarField& Yk = Y[k];

            fvScalarMatrix YEqn
            (
                fvm::ddt(Yk)
              + mvConvection->fvmDiv(phi, Yk)
              - fvm::laplacian(turbulence->muEff()/Sc*invrhoInf, Yk)
              ==
                invrhoInf*parcels.SYi(k, Yk)
              + fvOptions(Yk)
            );

            YEqn.relax();

            fvOptions.constrain(YEqn);

            YEqn.solve("Yk");

            fvOptions.correct(Yk);

            Yk.max(0.0);
            Yt += Yk;
        }
    }

    Y[inertIndex] = scalar(1) - Yt;
    Y[inertIndex].max(0.0);
}
