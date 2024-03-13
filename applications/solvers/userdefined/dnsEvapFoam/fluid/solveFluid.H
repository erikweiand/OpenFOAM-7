#include "globalProperties.H"
#include "UEqn.H"

// --- PISO loop
while (pisos.correct())
{
    #include "pEqn.H"
}

#include "YEqn.H"
#include "EEqn.H"

// update TKE and TDR
TKE = 1.0/2.0 * magSqr(U);
epsilon = 0.5 * nuGlobal * magSqr(twoSymm(fvc::grad(U)));

Info<< "\nEvolving " << parcels.name() << endl;
parcels.evolve();

// update charge density
rhoQFluid = parcels.rhoQ();
