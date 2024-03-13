// particle-laden channel flow
fvVectorMatrix UEqn
(
    fvm::ddt(U)
  + fvm::div(phi, U)
  - fvm::laplacian(nu, U)
 ==
    invrhoInf*parcels.SU(U)
    + fvOptions(U)
);

fvOptions.constrain(UEqn);

solve(UEqn == -fvc::grad(p));

fvOptions.correct(U);
