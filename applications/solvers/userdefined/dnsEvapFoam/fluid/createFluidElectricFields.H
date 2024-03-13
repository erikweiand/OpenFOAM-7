Info<< "    Adding to phiQ\n" << endl;
volScalarField phiQFluid
(
    IOobject
    (
        "phiQ",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("phiQ", dimensionSet(1,2,-3,0,0,-1,0), scalar(0.0))
);

// electric field strength
volVectorField ElectricFieldStrengthFluid
(
    IOobject
    (
        "ElectricFieldStrength",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("ElectricFieldStrength", dimensionSet(1,1,-3,0,0,-1,0), Foam::vector(0,0,0))
);
