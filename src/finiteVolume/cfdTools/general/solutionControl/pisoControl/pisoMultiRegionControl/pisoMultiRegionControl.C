/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "pisoMultiRegionControl.H"
#include "pisoControl.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pisoMultiRegionControl, 0);
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

const Foam::Time& Foam::pisoMultiRegionControl::time
(
    const PtrList<fvMesh>& pisoMeshes,
    const PtrList<fvMesh>& solidMeshes
)
{
    if (pisoMeshes.empty() && solidMeshes.empty())
    {
        FatalErrorInFunction
            << "There needs to be at least one region"
            << exit(FatalError);
    }

    if (!pisoMeshes.empty())
    {
        return pisoMeshes[0].time();
    }

    return solidMeshes[0].time();
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pisoMultiRegionControl::pisoMultiRegionControl
(
    PtrList<fvMesh>& pisoMeshes,
    PtrList<fvMesh>& solidMeshes,
    const word& algorithmName
)
:
    multiRegionSolutionControl(time(pisoMeshes, solidMeshes), algorithmName),
    pisoControls_(),
    solidControls_(),
    nCorrPiso_(-1),
    corrPiso_()
{
    forAll(pisoMeshes, i)
    {
        pisoControls_.append
        (
            new fluidSolutionControl(pisoMeshes[i], algorithmName)
        );

	corrPiso_.append(0);
    }

    forAll(solidMeshes, i)
    {
        solidControls_.append
        (
            new solidNoLoopControl(solidMeshes[i], algorithmName, *this)
        );
    }

    read();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pisoMultiRegionControl::~pisoMultiRegionControl()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::pisoMultiRegionControl::read()
{
    forAll(pisoControls_, i)
    {
        if (!pisoControls_[i].read())
        {
            return false;
        }
    }
    forAll(solidControls_, i)
    {
        if (!solidControls_[i].read())
        {
            return false;
        }
    }

    const dictionary& solutionDict = dict();

    nCorrPiso_ = solutionDict.lookupOrDefault<label>("nCorrectors", 1);

    return true;
}


bool Foam::pisoMultiRegionControl::isFinal() const
{
    bool orth = true;
    forAll(pisoControls_, i)
    {
	orth = orth &&
		((!pisoControls_[i].anyNonOrthogonalIter() && finalPisoIter(i))
     		|| (pisoControls_[i].finalNonOrthogonalIter() && finalPisoIter(i))
     		|| (pisoControls_[i].finalNonOrthogonalIter() && !anyPisoIter(i)));
	
	// shorten loop for false answers	
	if (!orth)
	{
	    return false;
	}
    }

    return orth;
}


bool Foam::pisoMultiRegionControl::correct()
{
    read();

    // check for final iter in all PISO controls
    bool isFin = true;
    forAll(pisoControls_, i)
    {
	isFin = isFin && finalPisoIter(i);
    }

    if (isFin)
    {
	forAll(pisoControls_, i)
	{
	    corrPiso_[i] = 0;

	    pisoControls_[i].updateFinal();
	}
	forAll(solidControls_, i)
	{
	    solidControls_[i].updateFinal();
	}

        return false;
    }

    forAll(pisoControls_, i)
    {
	++ corrPiso_[i];
	pisoControls_[i].updateFinal();
    }
    forAll(solidControls_, i)
    {
	solidControls_[i].updateFinal();
    }

    return true;
}


bool Foam::pisoMultiRegionControl::run(Time& time)
{
    return time.run();
}


bool Foam::pisoMultiRegionControl::loop(Time& time)
{
    return time.loop();
}


// ************************************************************************* //
