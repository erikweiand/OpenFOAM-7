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
    Foam::LorentzForce

Description
    Calculates particle Lorentz (electric field) force

SourceFiles
    LorentzForceI.H
    LorentzForce.C

\*---------------------------------------------------------------------------*/

#ifndef LorentzForce_H
#define LorentzForce_H

#include "ParticleForce.H"
#include "interpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;

/*---------------------------------------------------------------------------*\
                       Class LorentzForce Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class LorentzForce
:
    public ParticleForce<CloudType>
{
    // Private Data

        //- Name of electric field strength field - default = "ElectricFieldStrength"
        const word ElectricFieldStrengthName_;

        //- ElectricFieldStrength interpolator
        const interpolation<vector>* ElectricFieldStrengthInterpPtr_;


public:

    //- Runtime type information
    TypeName("lorentz");


    // Constructors

        //- Construct from mesh
        LorentzForce
        (
            CloudType& owner,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Construct copy
        LorentzForce(const LorentzForce& gf);

        //- Construct and return a clone
        virtual autoPtr<ParticleForce<CloudType>> clone() const
        {
            return autoPtr<ParticleForce<CloudType>>
            (
                new LorentzForce<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~LorentzForce();


    // Member Functions

        // Access

            //- Return the name of electric field strength field
            const word& ElectricFieldStrengthName() const;


        // Evaluation

            //- Cache fields
            virtual void cacheFields(const bool store);

            //- Calculate the non-coupled force
            virtual forceSuSp calcNonCoupled
            (
                const typename CloudType::parcelType& p,
                const typename CloudType::parcelType::trackingData& td,
                const scalar dt,
                const scalar mass,
                const scalar Re,
                const scalar muc
            ) const;

	    // Erik Weiand - 26/11/2019
	    //- Return the added electric field
            virtual vector elFieldAdd
            (
                const typename CloudType::parcelType& p,
                const typename CloudType::parcelType::trackingData& td
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "LorentzForceI.H"

#ifdef NoRepository
    #include "LorentzForce.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
