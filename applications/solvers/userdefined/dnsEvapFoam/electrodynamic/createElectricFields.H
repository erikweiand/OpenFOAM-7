Info<< "*** Reading electrodynamic mesh physical properties for region "
    << electricRegion.name() << endl;

Info<< "    Adding electric potential psiQ" << endl;
volScalarField phiQEl
(
    IOobject
    (
    	"phiQ",
    	runTime.timeName(),
   	electricRegion,
    	IOobject::MUST_READ,
    	IOobject::AUTO_WRITE
    ),
    electricRegion
);

Info<< "    Adding charge density rhoQ" << endl;
volScalarField rhoQEl
(
    IOobject
    (
    	"rhoQEl",
    	runTime.timeName(),
   	electricRegion,
    	IOobject::NO_READ,
    	IOobject::AUTO_WRITE
    ),
    electricRegion,
    dimensionedScalar("rhoQEl", dimensionSet(0, -3, 1, 0, 0, 1, 0), scalar(0.0))		// improve by reading initial charge density from fluid mesh
);

Info<< "    Adding electric field E" << endl;
volVectorField ElectricFieldStrengthEl
(
    IOobject
    (
    	"ElectricFieldStrengthEl",
    	runTime.timeName(),
    	electricRegion,
    	IOobject::NO_READ,
    	IOobject::AUTO_WRITE
    ),
    electricRegion,
    dimensionedVector("ElectricFieldStrengthEl", dimensionSet(1,1,-3,0,0,-1,0), Foam::vector(0,0,0))
);

// create new interpolation object for meshToMesh
meshToMesh interpolObjEl(mesh, electricRegion, meshToMesh::imDirect, false);
