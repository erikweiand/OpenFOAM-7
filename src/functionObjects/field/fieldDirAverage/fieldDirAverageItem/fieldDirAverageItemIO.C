/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

#include "fieldDirAverageItem.H"
#include "IOstreams.H"
#include "dictionaryEntry.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::fieldDirAverageItem::fieldDirAverageItem(Istream& is)
:
    fieldName_("unknown")
{
    is.check
    (
        "Foam::functionObjects::fieldDirAverageItem::fieldDirAverageItem"
        "(Foam::Istream&)"
    );

    const dictionaryEntry entry(dictionary::null, is);

    fieldName_ = entry.keyword();
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::functionObjects::operator>>
(
    Istream& is,
    fieldDirAverageItem& faItem
)
{
    is.check
    (
        "Foam::Istream& Foam::operator>>"
        "(Foam::Istream&, Foam::functionObjects::fieldDirAverageItem&)"
    );

    const dictionaryEntry entry(dictionary::null, is);

    faItem.fieldName_ = entry.keyword();

    return is;
}


Foam::Ostream& Foam::functionObjects::operator<<
(
    Ostream& os,
    const fieldDirAverageItem& faItem
)
{
    os.check
    (
        "Foam::Ostream& Foam::operator<<"
        "(Foam::Ostream&, const Foam::functionObjects::fieldDirAverageItem&)"
    );

    os  << faItem.fieldName_ << nl << token::BEGIN_BLOCK << nl;

    os  << token::END_BLOCK << nl;

    os.check
    (
        "Foam::Ostream& Foam::operator<<"
        "(Foam::Ostream&, const Foam::functionObjects::fieldDirAverageItem&)"
    );

    return os;
}


// ************************************************************************* //
