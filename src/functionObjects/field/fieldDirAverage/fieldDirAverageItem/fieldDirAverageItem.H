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

Class
    Foam::functionObjects::fieldDirAverageItem

Description
    Helper class to describe what form of averaging to apply.  A set will be
    applied to each base field in Foam::fieldDirAverage, of the form:

    \verbatim
    {
        mean            on;
        prime2Mean      on;
        base            time; // iteration
        window          200;  // optional averaging window
        windowName      w1;   // optional window name (default = "")
    }
    \endverbatim

    The averaging window corresponds to the averaging interval (iters or time)
    If not specified, the averaging is over 'all iters/time'

SourceFiles
    fieldDirAverageItem.C
    fieldDirAverageItemIO.C

\*---------------------------------------------------------------------------*/

#ifndef fieldDirAverageItem_H
#define fieldDirAverageItem_H

#include "NamedEnum.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Istream;
class Ostream;

namespace functionObjects
{

// Forward declaration of friend functions and operators
class fieldDirAverageItem;
Istream& operator>>(Istream&, fieldDirAverageItem&);
Ostream& operator<<(Ostream&, const fieldDirAverageItem&);

/*---------------------------------------------------------------------------*\
                      Class fieldDirAverageItem Declaration
\*---------------------------------------------------------------------------*/

class fieldDirAverageItem
{
public:

    // Public data

        // File and field name extensions

            //- Mean average
            static const word EXT_MEAN;

        //- Enumeration defining the averaging base type
        enum class baseType
        {
            iter,
            time
        };


private:

    // Private Data

        //- Field name
        word fieldName_;

public:

    // Constructors

        //- Construct null
        fieldDirAverageItem();

        //- Construct from Istream
        fieldDirAverageItem(Istream&);

        //- Copy constructor
        fieldDirAverageItem(const fieldDirAverageItem&);


    //- Destructor
    ~fieldDirAverageItem();


    // Member Functions

        // Access

            //- Return const access to the field name
            const word& fieldName() const
            {
                return fieldName_;
            }

    // Member Operators

        void operator=(const fieldDirAverageItem&);


    // Friend Operators

        friend bool operator==
        (
            const fieldDirAverageItem& a,
            const fieldDirAverageItem& b
        )
        {
            return
                a.fieldName_ == b.fieldName_;
        }

        friend bool operator!=
        (
            const fieldDirAverageItem& a,
            const fieldDirAverageItem& b
        )
        {
            return !(a == b);
        }


    // IOstream Operators

        friend Istream& operator>>(Istream&, fieldDirAverageItem&);
        friend Ostream& operator<<(Ostream&, const fieldDirAverageItem&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
