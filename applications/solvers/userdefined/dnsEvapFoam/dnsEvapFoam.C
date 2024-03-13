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

Application
    dnsEvapFoam

Description
    Solver for steady or transient fluid flow and solid heat conduction, with
    conjugate heat transfer between regions, buoyancy effects, turbulence,
    reactions and radiation modelling.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fixedGradientFvPatchFields.H"
#include "regionProperties.H"

#include "graph.H"
#include "writeFile.H"
#include "meshToMesh.H"
#include "basicReactingMultiphaseCloud.H"
#include "rhoReactionThermo.H"
#include "SLGThermo.H"
#include "fvOptions.H"
#include "localEulerDdtScheme.H"
#include "pisoMultiRegionControl.H"
#include "electromagneticConstants.H"
#include "turbulentFluidThermoModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #define NO_CONTROL
    #define CREATE_MESH createMeshesPostProcess.H
    #include "postProcess.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMeshes.H"
    #include "createFields.H"

    #include "initContinuityErrs.H"

    const label inertIndex(composition.species()[inertSpecie]);

    pisoMultiRegionControl pisos(fluidRegions, electricRegions);

    #include "createControl.H"

    // time step control
    #include "createTimeControls.H"

    // set pressure reference cell
    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p, pisos.dict(), pRefCell, pRefValue);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< nl << "Starting time loop" << endl;

    while (pisos.run(runTime))
    {
        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

	// courant time step control
	#include "CourantNo.H"
	#include "setDeltaT.H" 

        Info<< "\nSolving for fluid region "
            << mesh.name() << endl;
        #include "setRegionFluidFields.H"
        #include "solveFluid.H"

        Info<< "\nSolving for electrodynamic region "
            << electricRegion.name() << endl;
        #include "setRegionElectricFields.H"
        #include "solveElectric.H"

	Info << "pmax: " << max(p.internalField()) << nl << "pmin: " << min(p.internalField()) << endl;
	Info << "Emax: " << max(ElectricFieldStrength.internalField()) << nl << "Emin: " << min(ElectricFieldStrength.internalField()) << endl;
	Info << "Ymax: " << max(Y[1].internalField()) << nl << "Ymin: " << min(Y[1].internalField()) << endl;
	Info << "phiQmax: " << max(phiQ.internalField()) << nl << "phiQmin: " << min(phiQ.internalField()) << endl;
        Info << "Tmax: " << max(thermo.T()).value() << nl << "Tmin: " << min(thermo.T()).value() << endl;

	if(runTime.write())
	{
      	    p.write();
	}

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
