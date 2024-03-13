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

Class
    Foam::pisoMultiRegionControl

Description
    Piso multi-region control class. As Foam::pisoControl, but for a multi-
    region simulation comprising piso and solid regions. More region types
    could be added.

SourceFiles
    pisoMultiRegionControlI.H
    pisoMultiRegionControl.C

\*---------------------------------------------------------------------------*/

#ifndef pisoMultiRegionControl_H
#define pisoMultiRegionControl_H

#include "fluidSolutionControl.H"
#include "solidNoLoopControl.H"
#include "multiRegionSolutionControl.H"

#define PISOMULTIREGION_CONTROL

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class pisoMultiRegionControl Declaration
\*---------------------------------------------------------------------------*/

class pisoMultiRegionControl
:
    public multiRegionSolutionControl
{
protected:

    // Protected data

	//- The per-region piso controls
        PtrList<fluidSolutionControl> pisoControls_;

        //- The per-region solid controls
        PtrList<solidNoLoopControl> solidControls_;

    // Protected Static Functions

        //- Get the time from the region meshes
        static const Time& time
        (
            const PtrList<fvMesh>& pimpleMeshes,
            const PtrList<fvMesh>& solidMeshes
        );

        //- Maximum number of piso correctors
        label nCorrPiso_;

        //- Current piso corrector
        labelList corrPiso_;


public:

    // Static Data Members

        //- Run-time type information
        TypeName("pisoMultiRegionControl");


    // Constructors

        //- Construct from a meshes and the name of the algorithm
        pisoMultiRegionControl
	(
	    PtrList<fvMesh>& pisoMeshes,
	    PtrList<fvMesh>& solidMeshes,
	    const word& algorithmName="PISO"
	);


    //- Destructor
    virtual ~pisoMultiRegionControl();


    // Member Functions

        // IO

            //- Read controls
            virtual bool read();


        // Access

	    //- Access the i-th pimple control
            inline fluidSolutionControl& piso(const label i);

            //- Access the i-th solid control
            inline solidNoLoopControl& solid(const label i);

            //- Maximum number of piso correctors
            inline label nCorrPiso() const;

            //- Flag to indicate any piso iteration
            inline bool anyPisoIter(const label i) const;

            //- Flag to indicate the first piso iteration
            inline bool firstPisoIter(const label i) const;

            //- Flag to indicate the last piso iteration
            inline bool finalPisoIter(const label i) const;

            //- Flag to indicate the last inner iteration (last piso and last
            //  non-orthogonal)
            inline bool finalInnerIter(const label i) const;

            //- Flag to indicate whether in final state
            virtual bool isFinal() const;


        // Evolution

            //- Piso loop
            bool correct();

            //- Time run loop
            bool run(Time& time);

            //- Time loop loop
            bool loop(Time& time);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "pisoMultiRegionControlI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
