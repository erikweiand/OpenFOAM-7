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

Class
    Foam::interpolationPatchCellConstrained

Description
    Uses the cell center and vertex values for any point in the cell apart 
    from a boundary face where it uses the cell value directly.
    Note: will not work on an empty patch.

\*---------------------------------------------------------------------------*/

#ifndef interpolationPatchCellConstrained_H
#define interpolationPatchCellConstrained_H

#include "interpolation.H"
#include "cellPointWeight.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;

/*---------------------------------------------------------------------------*\
                           Class interpolationPatchCellConstrained Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class interpolationPatchCellConstrained
:
    public interpolation<Type>
{
protected:

    // Protected data

        //- Interpolated volfield
        const GeometricField<Type, pointPatchField, pointMesh> psip_;


public:

    //- Runtime type information
    TypeName("patchCellConstrained");


    // Constructors

        //- Construct from components
        interpolationPatchCellConstrained
        (
            const GeometricField<Type, fvPatchField, volMesh>& psi
        );

        //- Construct from components
        interpolationPatchCellConstrained
        (
            const GeometricField<Type, fvPatchField, volMesh>& psi,
            tmp<GeometricField<Type, pointPatchField, pointMesh>> psip
        );


    // Member Functions

	//- Interpolate field to the given coordinates in the tetrahedron
        //  defined by the given indices.
        inline Type interpolate
        (
            const vector& position,
    	    const label celli,
            const label facei = -1
        ) const;

        //- Interpolate field for the given cellPointWeight
        inline Type interpolateCP(const cellPointWeight& cpw) const;

        //- Interpolate field to the given point in the given cell
        inline Type interpolateCP
        (
            const vector& position,
            const label celli,
            const label facei = -1
        ) const;

        //- Interpolate field to the given coordinates in the tetrahedron
        //  defined by the given indices.
        inline Type interpolateCP
        (
            const barycentric& coordinates,
            const tetIndices& tetIs,
            const label facei = -1
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "interpolationPatchCellConstrainedI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "interpolationPatchCellConstrained.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
