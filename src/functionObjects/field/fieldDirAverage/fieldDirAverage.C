/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "fieldDirAverage.H"
#include "volFields.H"
#include "fieldDirAverageItem.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(fieldDirAverage, 0);
    addToRunTimeSelectionTable(functionObject, fieldDirAverage, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

// empty

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::fieldDirAverage::fieldDirAverage
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    faItems_(),
    channelIndexing_(mesh_ ,dict)
{
    read(dict);

    gFormat_ = runTime.graphFormat();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::fieldDirAverage::~fieldDirAverage()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::fieldDirAverage::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    Log << type() << " " << name() << ":" << nl;

    dict.lookup("fields") >> faItems_;

    Log << endl;

    return true;
}


bool Foam::functionObjects::fieldDirAverage::execute()
{
    return true;
}


bool Foam::functionObjects::fieldDirAverage::write()
{
    const Time& runTime = obr_.time();

//    if (Pstream::master())
//    {
    Info<< "Collapsing fields for time " << runTime.timeName() << endl;

    // Average fields over channel down to a line
    #include "collapse.H"
//    }

    return true;
}


// ************************************************************************* //
