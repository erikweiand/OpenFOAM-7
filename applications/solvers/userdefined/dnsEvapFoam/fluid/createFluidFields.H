Info<< "*** Reading fluid mesh properties for region "
    << mesh.name() << nl << endl;

// transport properties
Info<< "    Adding to transportProperties" << endl;
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

dimensionedScalar nuValue
(
    "nu",
    dimViscosity,
    transportProperties
);

dimensionedScalar Sc
(
    "Sc",
    dimensionSet(0, 0, 0, 0, 0, 0, 0),
    transportProperties
);

IOdictionary ScDict
(
    IOobject
    (
        "ScDict",
        runTime.constant(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    )
);

// move to runtime if change desired
ScDict.set("ScVal", Sc);

dimensionedScalar Pr
(
    "Pr",
    dimensionSet(0, 0, 0, 0, 0, 0, 0),
    transportProperties
);

dimensionedScalar rhoInfValue
(
    "rhoInf",
    dimDensity,
    transportProperties
);
dimensionedScalar invrhoInf("invrhoInf",(1.0/rhoInfValue));

// thermophysical properties
Info<< "    Adding to thermophysical properties\n" << endl;
autoPtr<rhoReactionThermo> pThermo(rhoReactionThermo::New(mesh));
rhoReactionThermo& thermo = pThermo();
thermo.validate(args.executable(), "h", "e");

SLGThermo slgThermo(mesh, thermo);
basicSpecieMixture& composition = thermo.composition();
PtrList<volScalarField>& Y = composition.Y();

const word inertSpecie(thermo.lookup("inertSpecie"));
if (!composition.species().found(inertSpecie))
{
    FatalIOErrorIn(args.executable().c_str(), thermo)
    	<< "Inert specie " << inertSpecie << " not found in available species "
    	<< composition.species()
    	<< exit(FatalIOError);
}

// turbulence properties
Info<< "    Adding to turbulenceProperties" << endl;
IOdictionary turbulenceProperties
(
    IOobject
    (
    	"turbulenceProperties",
    	runTime.constant(),
    	mesh,
    	IOobject::MUST_READ_IF_MODIFIED,
    	IOobject::NO_WRITE
    )
);

// pressure
Info<< "    Adding to p" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

// velocity
Info<< "    Adding to U" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

    // density
Info<< "    Adding to rho" << endl;
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    rhoInfValue
);

// surface flux
Info<< "    Adding to phi" << endl;
surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    fvc::flux(U)
);

// turbulence model placeholder
Info<< "    Adding to turbulence model" << endl;
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);

// gravitational constant
Info<< "    Adding to g" << endl;
uniformDimensionedVectorField g
(
    IOobject
    (
        "g",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

// dynamic viscosity
Info<< "    Adding to mu" << endl;
volScalarField mu
(
    IOobject
    (
        "mu",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    nuValue*rhoInfValue
);

// kinematic viscosity
Info<< "    Adding to nu" << endl;
volScalarField nu
(
    IOobject
    (
        "nu",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    nuValue
);
   
// nuGlobal
Info<< "    Adding to nuGlobal" << endl;
dimensionedScalar nuGlobal
(
    "nuGlobal", nuValue
);

volVectorField enthDiff
(
    IOobject
    (
        "enthDiff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("enthDiff", dimensionSet(0, 3, -3, 0, 0, 0, 0), vector::zero)
);

// kinematic cloud
const word reactingMultiphaseCloudName
(
    args.optionLookupOrDefault<word>("cloudName","reactingMultiphaseCloud")
);

// reacting multiphase cloud
Info<< "    Adding to reactingMultiphaseCloud" << endl;
basicReactingMultiphaseCloud parcels
( 
    reactingMultiphaseCloudName,
    rho,
    U,
    g,
    slgThermo
);

// charge density
Info<< "    Adding to rhoQ" << endl;
volScalarField rhoQFluid
(
    IOobject
    (
        "rhoQ",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("rhoQ", dimensionSet(0, -3, 1, 0, 0, 1, 0), scalar(0.0))
);

// TKE
Info<< "    Adding to k" << endl;
volScalarField TKE
(
    IOobject
    (
        "TKE",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("TKE", dimensionSet(0, 2, -2, 0, 0, 0, 0), scalar(0.0))
);

// TDR
Info<< "    Adding to epsilon" << endl;
volScalarField epsilon
(
    IOobject
    (
        "epsilon",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("epsilon", dimensionSet(0, 2, -3, 0, 0, 0, 0), scalar(0.0))
);

// create new interpolation object for meshToMesh
meshToMesh interpolObjFluid(electricRegion, mesh, meshToMesh::imDirect, false);

mesh.setFluxRequired(p.name());

multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;
forAll(Y, k)
{
    fields.add(Y[k]);
}
fields.add(thermo.he());

#include "createFvOptions.H"
