    const wordList electricsNames(rp["electric"]);

    Info<< "Create electric mesh for region " << electricsNames[0]
        << " for time = " << runTime.timeName() << nl << endl;

    PtrList<fvMesh> electricRegions(electricsNames.size());

    electricRegions.set
    (
	0,
	new fvMesh
	(
 	    IOobject
	    (
            	electricsNames[0],
            	runTime.timeName(),
            	runTime,
            	IOobject::MUST_READ
            )
	)
    );

    fvMesh& electricRegion = electricRegions[0];
