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
    Foam::NanofluidEvaporation

Description
    Nanofluid evaporation model
    - As used in Wei et al. (2016).
    - Implementation by Erik Weiand, December 2019.

\*---------------------------------------------------------------------------*/

#ifndef NanofluidEvaporation_H
#define NanofluidEvaporation_H

#include "PhaseChangeModel.H"
#include "liquidMixtureProperties.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                     Class NanofluidEvaporation Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class NanofluidEvaporation
:
    public PhaseChangeModel<CloudType>
{
protected:

    // Protected data

        //- Global liquid properties data
        const liquidMixtureProperties& liquids_;

        //- List of active liquid names
        List<word> activeLiquids_;
	
        //- Mapping between liquid and carrier species
        List<label> liqToCarrierMap_;

        //- Mapping between local and global liquid species
        List<label> liqToLiqMap_;

	//- Sine of contact angle squared
	scalar alpha_;

	//- nanoparticle radius
	scalar rP_;

	//- initial nanoparticle mass fraction per droplet
	scalar Yvo_;

	//- liquid dynamic viscosity
	scalar etal_;
	

    // Protected Member Functions

        //- Sherwood number as a function of Reynolds and Schmidt numbers
        scalar Sh(const scalar Re, const scalar Sc) const;

        //- Calculate the carrier phase component volume fractions at celli
        tmp<scalarField> calcXc(const label celli) const;


public:

    //- Runtime type information
    TypeName("nanofluidEvaporation");


    // Constructors

        //- Construct from dictionary
        NanofluidEvaporation(const dictionary& dict, CloudType& cloud);

        //- Construct copy
        NanofluidEvaporation(const NanofluidEvaporation<CloudType>& pcm);

        //- Construct and return a clone
        virtual autoPtr<PhaseChangeModel<CloudType>> clone() const
        {
            return autoPtr<PhaseChangeModel<CloudType>>
            (
                new NanofluidEvaporation<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~NanofluidEvaporation();


    // Member Functions

        //- Update model
        virtual void calculate
        (
            const scalar dt,
            const label celli,
            const scalar Re,
            const scalar Pr,
            const scalar d,
            const scalar nu,
            const scalar T,
            const scalar Ts,
            const scalar pc,
            const scalar Tc,
            const scalarField& X,
            scalarField& dMassPC
        ) const {}

        //- Update model, Abramzon & Sirignano  //  M.P. Sitte, 09/10/2015
        virtual void calculateAS
        (
            const scalar dt,
            const label cellI,
            const scalar Re,
            const scalar Pr,
            const scalar d,
            const scalar rhos, // M.P.S.
            const scalar nu,
            const scalar T,
            const scalar Ts,
            const scalar pc,
            const scalar Tc,
            const scalarField& Yl,
	    const scalar d0,	  // E. Weiand - 2020/03/04
	    scalar& BM,	  // E. Weiand - 2020/03/30
            scalarField& dMassPC, // M.P.S.
            scalar& fmixS,        // M.P.S.
            scalar& grFmixS,      // M.P.S.
            scalar& DS,           // M.P.S.
	    scalar& Kevap,	  // E. Weiand - 2020/03/04
	    scalar& Dp,		  // E. Weiand - 2020/03/04
	    scalar& beta,	  // E. Weiand - 2020/03/04
	    scalar& tauD,	  // E. Weiand - 2020/03/06
	    bool& isDenselyPacked,	// E. Weiand - 2020/09/21
	    scalar& f0,
	    scalar& dtPrev,
	    bool& isFirstStep
        ) const;

        //- Return the enthalpy per unit mass
        virtual scalar dh
        (
            const label idc,
            const label idl,
            const scalar p,
            const scalar T
        ) const;

        //- Return vapourisation temperature
        virtual scalar Tvap(const scalarField& X) const;

        //- Return maximum/limiting temperature
        virtual scalar TMax(const scalar p, const scalarField& X) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "NanofluidEvaporation.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
