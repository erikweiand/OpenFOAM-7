    Info<< "Reading transportProperties\n" << endl;

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    dimensionedScalar nu
    (
        "nu",
        dimViscosity,
        transportProperties
    );

    dimensionedScalar rhoInfValue
    (
	"rhoInf",
	dimDensity,
	transportProperties
    );
    dimensionedScalar invrhoInf("invrhoInf",(1.0/rhoInfValue));
