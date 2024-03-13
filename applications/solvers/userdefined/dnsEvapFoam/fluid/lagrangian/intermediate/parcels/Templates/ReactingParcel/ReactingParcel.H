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
    Foam::ReactingParcel

Description
    Reacting parcel class with one/two-way coupling with the continuous
    phase.

SourceFiles
    ReactingParcelI.H
    ReactingParcel.C
    ReactingParcelIO.C

\*---------------------------------------------------------------------------*/

#ifndef ReactingParcel_H
#define ReactingParcel_H

#include "particle.H"
#include "SLGThermo.H"
#include "demandDrivenEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class ParcelType>
class ReactingParcel;

template<class ParcelType>
Ostream& operator<<
(
    Ostream&,
    const ReactingParcel<ParcelType>&
);


/*---------------------------------------------------------------------------*\
                        Class ReactingParcel Declaration
\*---------------------------------------------------------------------------*/

template<class ParcelType>
class ReactingParcel
:
    public ParcelType
{
    // Private Data

        //- Size in bytes of the fields
        static const std::size_t sizeofFields_;


public:

    //- Class to hold reacting parcel constant properties
    class constantProperties
    :
        public ParcelType::constantProperties
    {
        // Private Data

            //- Minimum pressure [Pa]
            demandDrivenEntry<scalar> pMin_;

            //- Constant volume flag - e.g. during mass transfer
            demandDrivenEntry<bool> constantVolume_;


    public:

        // Constructors

            //- Null constructor
            constantProperties();

            //- Copy constructor
            constantProperties(const constantProperties& cp);

            //- Construct from dictionary
            constantProperties(const dictionary& parentDict);


        // Access

            //- Return const access to the minimum pressure
            inline scalar pMin() const;

            //- Return const access to the constant volume flag
            inline bool constantVolume() const;
    };


    class trackingData
    :
        public ParcelType::trackingData
    {
    private:

        // Private Data

            // Interpolators for continuous phase fields

                //- Interpolator for continuous phase pressure field
                autoPtr<interpolation<scalar>> pInterp_;


            // Cached continuous phase properties

                //- Pressure [Pa]
                scalar pc_;

    public:

        typedef typename ParcelType::trackingData::trackPart trackPart;

        // Constructors

            //- Construct from components
            template<class TrackCloudType>
            inline trackingData
            (
                const TrackCloudType& cloud,
                trackPart part = ParcelType::trackingData::tpLinearTrack
            );


        // Member Functions

            //- Return const access to the interpolator for continuous phase
            //  pressure field
            inline const interpolation<scalar>& pInterp() const;

            //- Return the continuous phase pressure
            inline scalar pc() const;

            //- Access the continuous phase pressure
            inline scalar& pc();
    };


protected:

    // Protected data

        // Parcel properties

            //- Initial mass [kg]
            scalar mass0_;

            //- Mass fractions of mixture []
            scalarField Y_;

	    // E. Weiand - 2020/02/21
	    //- Schmidt number []
	    autoPtr<scalar> Sc_;

	    // E. Weiand - 2020/03/04 
	    // Initial diameter [m^2/s]
	    scalar d0_;

	    // E. Weiand - 2020/03/04 
	    // Evaporation rate [m^2/s]
	    scalar Kevap_;
	
	    // E. Weiand - 2020/03/04
	    // Particle diffusivity [m^2/s]
	    scalar Dp_;

	    // E. Weiand - 2020/03/04
	    // Droplet surface blockage [-]
	    scalar beta_;

	    // E. Weiand - 2020/09/29
	    bool isDenselyPacked_;

	    // E. Weiand - 2020/09/29
	    scalar densMult_;

	    //- Adams-Bashforth variables
    	    scalar f0_;
	    scalar f0Q_;
 	    scalar dtPrev_;
 	    scalar dtPrevQ_;
	    bool isFirstStep_;
	    bool isFirstStepQ_;


    // Protected Member Functions

	// M.P.Sitte calcSurfaceValuesAS in .H
        //- Calculate surface thermo properties
              template<class TrackCloudType>
              void calcSurfaceValuesAS
              (
		  TrackCloudType& cloud,
                  trackingData& td,
                  const label cellI,
                  const scalar T,
                  const label idPhase,             // M.P.Sitte
                  const scalarField& YComponents,  // M.P.Sitte
                  scalar& Ts,
                  scalarField& YYs,
                  scalar& Ws,
                  scalar& rhos,
                  scalar& mus,
                  scalar& Pr,
                  scalar& kappas,
                  scalar& BM
              ) const;

	//- Calculate Phase change, modified.
        // M.P. Sitte: <calcPhaseChangeAS>
        template<class TrackCloudType>
        void calcPhaseChangeAS
        (
	    TrackCloudType& cloud,
            trackingData& td,
            const scalar dt,           // timestep
            const label cellI,         // owner cell
            const scalar Re,           // Reynolds number
            const scalar Pr,           // Prandtl number
            const scalar Ts,           // Surface temperature
            const scalar rhos,         // Surface density   // M.P. Sitte, 09/10/2015
            const scalar nus,          // Surface kinematic viscosity
            const scalar d,            // diameter
            const scalar T,            // temperature
            const scalar mass,         // mass
            const label idPhase,       // id of phase involved in phase change
            const scalar YPhase,       // total mass fraction
            const scalarField& YComponents, // component mass fractions
	    scalar& BM,		// Spalding mass transfer number
            scalarField& dMassPC,      // mass transfer - local to particle
            scalar& fmixS,             // ....... M.P. Sitte, 12/10/2015
            scalar& grFmixS,      // ....... M.P. Sitte, 09/10/2015
            scalar& DS,           // ....... M.P. Sitte, 12/10/2015
            scalar& Sh,                // explicit particle enthalpy source
            scalar& N,                 // flux of species emitted from particle
            scalar& NCpW,              // sum of N*Cp*W of emission species
            scalarField& Cs            // carrier conc. of emission species
        );

        //- Calculate new particle temperature
        // calcHeatTransferAS in .H // M.P.Sitte
        template<class TrackCloudType>
        scalar calcHeatTransferAS
        (
	    TrackCloudType& cloud,
            trackingData& td,
            const scalar dt,           // timestep
            const label cellI,         // owner cell
            const label idPhase,
            const scalar Re,           // Reynolds number
            const scalar Pr,           // Prandtl number - surface
            const scalar Ts,           // ref. temperature
            const scalar rhos,         // density @ ref. cond              // Sitte
            const scalar nus,           // viscosity @ ref. cond
            const scalar Cps,           // gas Cp @ ref. cond.             // Sitte
            const scalar kappas,        // Thermal conductivity - surface
            const scalar BM,           // mass tranfer number             // Sitte
            const scalarField& dMass,  // change in mass, due to evap.    // Sitte
            const scalar Sh,           // explicit particle enthalpy source
            scalar& dhsTrans,          // sensible enthalpy transfer to carrier
            scalar& Sph                // linearised heat transfer coefficient
        );

        //- Calculate Phase change
        template<class TrackCloudType>
        void calcPhaseChange
        (
            TrackCloudType& cloud,
            trackingData& td,
            const scalar dt,           // timestep
            const scalar Re,           // Reynolds number
            const scalar Pr,           // Prandtl number
            const scalar Ts,           // Surface temperature
            const scalar nus,          // Surface kinematic viscosity
            const scalar d,            // diameter
            const scalar T,            // temperature
            const scalar mass,         // mass
            const label idPhase,       // id of phase involved in phase change
            const scalar YPhase,       // total mass fraction
            const scalarField& YComponents, // component mass fractions
            scalarField& dMassPC,      // mass transfer - local to parcel
            scalar& Sh,                // explicit parcel enthalpy source
            scalar& N,                 // flux of species emitted from parcel
            scalar& NCpW,              // sum of N*Cp*W of emission species
            scalarField& Cs            // carrier conc. of emission species
        );

        //- Update mass fraction
        scalar updateMassFraction
        (
            const scalar mass0,
            const scalarField& dMass,
            scalarField& Y
        ) const;


public:

    // Static Data Members

        //- Runtime type information
        TypeName("ReactingParcel");

        //- String representation of properties
        AddToPropertyList
        (
            ParcelType,
            " mass0"
          + " nPhases(Y1..YN)"
	  + " d0"
	  + " Kevap"
	  + " Dp"
	  + " beta"
	  + " isDenselyPacked"
	  + " densMult"
        );


    // Constructors

        //- Construct from mesh, coordinates and topology
        //  Other properties initialised as null
        inline ReactingParcel
        (
            const polyMesh& mesh,
            const barycentric& coordinates,
            const label celli,
            const label tetFacei,
            const label tetPti
        );

        //- Construct from a position and a cell, searching for the rest of the
        //  required topology. Other properties are initialised as null.
        inline ReactingParcel
        (
            const polyMesh& mesh,
            const vector& position,
            const label celli
        );

        //- Construct from components
        inline ReactingParcel
        (
            const polyMesh& mesh,
            const barycentric& coordinates,
            const label celli,
            const label tetFacei,
            const label tetPti,
            const label typeId,
            const scalar nParticle0,
            const scalar d0,
            const scalar dTarget0,
            const vector& U0,
            const vector& f0,
            const vector& angularMomentum0,
            const vector& torque0,
            const scalarField& Y0,
	    const scalar Qdens0,	// E. Weiand - 27/02/2020 - electric particle charge density
	    const vector& E0,		// E. Weiand - 30/10/2019 - electric field strength
	    const vector& FLorentz0,	// E. Weiand - 30/10/2019 - Lorentz force
            const constantProperties& constProps
        );

        //- Construct from Istream
        ReactingParcel
        (
            const polyMesh& mesh,
            Istream& is,
            bool readFields = true
        );

        //- Construct as a copy
        ReactingParcel
        (
            const ReactingParcel& p,
            const polyMesh& mesh
        );

        //- Construct as a copy
        ReactingParcel(const ReactingParcel& p);

        //- Construct and return a (basic particle) clone
        virtual autoPtr<particle> clone() const
        {
            return autoPtr<particle>(new ReactingParcel<ParcelType>(*this));
        }

        //- Construct and return a (basic particle) clone
        virtual autoPtr<particle> clone(const polyMesh& mesh) const
        {
            return autoPtr<particle>
            (
                new ReactingParcel<ParcelType>(*this, mesh)
            );
        }

        //- Factory class to read-construct particles used for
        //  parallel transfer
        class iNew
        {
            const polyMesh& mesh_;

        public:

            iNew(const polyMesh& mesh)
            :
                mesh_(mesh)
            {}

            autoPtr<ReactingParcel<ParcelType>> operator()(Istream& is) const
            {
                return autoPtr<ReactingParcel<ParcelType>>
                (
                    new ReactingParcel<ParcelType>(mesh_, is, true)
                );
            }
        };


    // Member Functions

        // Access

            //- Return const access to initial mass [kg]
            inline scalar mass0() const;

            //- Return const access to mass fractions of mixture []
            inline const scalarField& Y() const;
	
	    // E. Weiand - 2020/03/04
	    //- Return const access to initial diameter [m]
            inline scalar d0() const;

	    // E. Weiand - 2020/03/04
	    //- Return const access to evaporation rate [m^2/s]
            inline scalar Kevap() const;

	    // E. Weiand - 2020/03/04
	    //- Return const access to particle diffusivity [m^2/s]
            inline scalar Dp() const;

	    // E. Weiand - 2020/03/04
            //- Return const access to droplet surface blockage [-]
	    inline scalar beta() const;

	    // E. Weiand - 2020/09/29
            inline bool isDenselyPacked() const;
	    
	    // E. Weiand - 2020/09/29
            inline scalar densMult() const;

	    inline scalar f0() const;
	    inline scalar dtPrev() const;
	    inline bool isFirstStep() const;


        // Edit

            //- Return access to initial mass [kg]
            inline scalar& mass0();

            //- Return access to mass fractions of mixture []
            inline scalarField& Y();

	    // E. Weiand - 2020/03/04
            //- Return access to initial diameter [m]
            inline scalar& d0();

	    // E. Weiand - 2020/03/04
            //- Return access to evaporation rate [m^2/s]
            inline scalar& Kevap();

	    // E. Weiand - 2020/03/04
            //- Return access to particle diffusivity [m^2/s]
            inline scalar& Dp();
	
	    // E. Weiand - 2020/03/04
            //- Return access to droplet surface blockage [-]
	    inline scalar& beta();

	    // E. Weiand - 2020/09/29
	    inline bool& isDenselyPacked();
	    
	    // E. Weiand - 2020/09/29
	    inline scalar& densMult();
	    
	    inline scalar& f0();
	    inline scalar& dtPrev();
	    inline bool& isFirstStep();


        // Main calculation loop

            //- Set cell values
            template<class TrackCloudType>
            void setCellValues(TrackCloudType& cloud, trackingData& td);

            //- Correct cell values using latest transfer information
            template<class TrackCloudType>
            void cellValueSourceCorrection
            (
                TrackCloudType& cloud,
                trackingData& td,
                const scalar dt
            );

            //- Correct surface values due to emitted species
            template<class TrackCloudType>
            void correctSurfaceValues
            (
                TrackCloudType& cloud,
                trackingData& td,
                const scalar T,
                const scalarField& Cs,
                scalar& rhos,
                scalar& mus,
                scalar& Prs,
                scalar& kappas
            );

	    // M.P. Sitte, 05/04/2016: correctSurfaceValuesAS
            //- Correct surface values due to emitted species
            template<class TrackCloudType>
            void correctSurfaceValuesAS
            (
    	        TrackCloudType& cloud,
                trackingData& td,
                const scalar T,           // film temperature, ref. cond.
                const scalarField& Ys,   // film composition, ref. cond.    // M.P.Sitte
                const scalar Ws,          // mean molecular mass, ref cond.  // M.P.Sitte
                scalar& rhos,
                scalar& mus,
                scalar& Prs,
                scalar& kappas,
                scalar& Cps
            );

            //- Update parcel properties over the time interval
            template<class TrackCloudType>
            void calc
            (
                TrackCloudType& cloud,
                trackingData& td,
                const scalar dt
            );


        // I-O

            //- Read
            template<class CloudType, class CompositionType>
            static void readFields
            (
                CloudType& c,
                const CompositionType& compModel
            );

            //- Read - no composition
            template<class CloudType>
            static void readFields(CloudType& c);

            //- Write
            template<class CloudType, class CompositionType>
            static void writeFields
            (
                const CloudType& c,
                const CompositionType& compModel
            );

            //- Write - composition supplied
            template<class CloudType>
            static void writeFields(const CloudType& c);


    // Ostream Operator

        friend Ostream& operator<< <ParcelType>
        (
            Ostream&,
            const ReactingParcel<ParcelType>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ReactingParcelI.H"
#include "ReactingParcelTrackingDataI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ReactingParcel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
