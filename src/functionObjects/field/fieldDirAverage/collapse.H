    fileName path
    (
        Pstream::parRun()
          ? runTime.path()/".."/"postProcessing"/"graphs"/mesh_.time().timeName()
          : runTime.path()/"postProcessing"/"graphs"/mesh_.time().timeName()
    );
    mkDir(path);

    const scalarField& y = channelIndexing_.y();

    forAll(faItems_, i)
    {
	const word& fieldName = faItems_[i].fieldName();

	const word meanFieldName = fieldName + "DirMean";

    	if (obr_.foundObject<volScalarField>(fieldName))
    	{
            const volScalarField baseField = obr_.lookupObject<volScalarField>(fieldName);

            scalarField meanValues(channelIndexing_.collapse(baseField));

	    if (Pstream::master())
            {
                makeGraph(y, meanValues, meanFieldName, path, gFormat_);
	    }
    	}
    }
