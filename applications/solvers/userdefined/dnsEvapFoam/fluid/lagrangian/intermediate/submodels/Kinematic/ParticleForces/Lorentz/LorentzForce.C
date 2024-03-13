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

#include "LorentzForce.H"
#include "demandDrivenData.H"
#include "electromagneticConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::LorentzForce<CloudType>::LorentzForce
(
    CloudType& owner,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    ParticleForce<CloudType>(owner, mesh, dict, typeName, true),
    ElectricFieldStrengthName_
    (
        this->coeffs().template lookupOrDefault<word>("ElectricFieldStrength", "ElectricFieldStrength")
    ),
    ElectricFieldStrengthInterpPtr_(nullptr)
{}


template<class CloudType>
Foam::LorentzForce<CloudType>::LorentzForce
(
    const LorentzForce& pf
)
:
    ParticleForce<CloudType>(pf),
    ElectricFieldStrengthName_(pf.ElectricFieldStrengthName_),
    ElectricFieldStrengthInterpPtr_(pf.ElectricFieldStrengthInterpPtr_)
{}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::LorentzForce<CloudType>::~LorentzForce()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::LorentzForce<CloudType>::cacheFields(const bool store)
{
    if (store)
    {
        const volVectorField& ElectricFieldStrength =
            this->mesh().template lookupObject<volVectorField>(ElectricFieldStrengthName_);

        ElectricFieldStrengthInterpPtr_ = interpolation<vector>::New
        (
            this->owner().solution().interpolationSchemes(),
            ElectricFieldStrength
        ).ptr();
    }
    else
    {
        deleteDemandDrivenData(ElectricFieldStrengthInterpPtr_);
    }
}


template<class CloudType>
Foam::forceSuSp Foam::LorentzForce<CloudType>::calcNonCoupled
(
    const typename CloudType::parcelType& p,
    const typename CloudType::parcelType::trackingData& td,
    const scalar dt,
    const scalar mass,
    const scalar Re,
    const scalar muc
) const
{
    forceSuSp value(Zero, 0.0);

    value.Su() = p.FLorentz();

    return value;
}


// Erik Weiand - 26/11/2019
template<class CloudType>
Foam::vector Foam::LorentzForce<CloudType>::elFieldAdd
(
    const typename CloudType::parcelType& p,
    const typename CloudType::parcelType::trackingData& td
) const
{
    const interpolation<vector>& ElectricFieldStrengthInterp = *ElectricFieldStrengthInterpPtr_;
    const vector interp = ElectricFieldStrengthInterp.interpolate(p.position(), p.cell());

    return interp;
}


// ************************************************************************* //
