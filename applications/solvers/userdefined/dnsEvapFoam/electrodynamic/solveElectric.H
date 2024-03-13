{
    fvScalarMatrix phiQEqn
    (
        fvm::laplacian(phiQ) + rhoQ()/constant::electromagnetic::epsilon0
    );

    phiQEqn.solve();

    ElectricFieldStrength == -fvc::grad(phiQ);
}
