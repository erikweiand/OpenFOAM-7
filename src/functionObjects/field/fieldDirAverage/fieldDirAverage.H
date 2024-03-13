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
    Foam::functionObjects::fieldDirAverage

Description
    Calculates average in specified carthesian direction.

    Example of function object specification:
    \verbatim
    fieldDirAverage1
    {
        type                fieldDirAverage;
        libs                ("libfieldFunctionObjects.so");

        writeControl        writeTime;

	symmetric	    no;
	component	    y;
	patches		    wallA;

        fields
        (
            U
            p
        );
    }
    \endverbatim

Usage
    \table
        Property          | Description               | Required | Default
        type              | type name: fieldDirAverage              | yes |
	symmetric	  | domain symmetry			    | yes |
	component	  | averaging component			    | yes |
	patches		  | seed patches			    | yes |
        fields            | list of fields			    | yes |
    \endtable


Note
    To employ the \c prime2Mean option, the \c mean option must be selected.

See also
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObject

SourceFiles
    fieldDirAverage.C
    fieldDirAverageTemplates.C
    fieldDirAverageItem.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_fieldDirAverage_H
#define functionObjects_fieldDirAverage_H

#include "fvMeshFunctionObject.H"
#include "channelIndex.H"
#include "makeGraph.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

// Forward declaration of classes
class fieldDirAverageItem;

/*---------------------------------------------------------------------------*\
                        Class fieldDirAverage Declaration
\*---------------------------------------------------------------------------*/

class fieldDirAverage
:
    public fvMeshFunctionObject
{
protected:

    // Protected data

        //- Graph format
        word gFormat_;

        //- List of field average items, describing what averages to be
        //  calculated and output
        List<fieldDirAverageItem> faItems_;

	//- Channel indexing
	channelIndex channelIndexing_;

        // Counters

	    // *empty*


    // Protected Member Functions

        // Initialisation routines
	    
	    // *empty*

        // Calculation functions

	    // *empty*

        // I-O

	    // *empty*


public:

    //- Runtime type information
    TypeName("fieldDirAverage");


    // Constructors

        //- Construct from Time and dictionary
        fieldDirAverage
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        fieldDirAverage(const fieldDirAverage&) = delete;


    //- Destructor
    virtual ~fieldDirAverage();


    // Member Functions

        //- Read the field average data
        virtual bool read(const dictionary&);

        //- Calculate the field averages
        virtual bool execute();

        //- Write the field averages
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const fieldDirAverage&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
