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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
inline
Foam::ReactingParcel<ParcelType>::constantProperties::constantProperties()
:
    ParcelType::constantProperties(),
    pMin_(this->dict_, 0.0),
    constantVolume_(this->dict_, false)
{}


template<class ParcelType>
inline Foam::ReactingParcel<ParcelType>::constantProperties::constantProperties
(
    const constantProperties& cp
)
:
    ParcelType::constantProperties(cp),
    pMin_(cp.pMin_),
    constantVolume_(cp.constantVolume_)
{}


template<class ParcelType>
inline Foam::ReactingParcel<ParcelType>::constantProperties::constantProperties
(
    const dictionary& parentDict
)
:
    ParcelType::constantProperties(parentDict),
    pMin_(this->dict_, "pMin", 1000.0),
    constantVolume_(this->dict_, word("constantVolume"))
{}


template<class ParcelType>
inline Foam::ReactingParcel<ParcelType>::ReactingParcel
(
    const polyMesh& mesh,
    const barycentric& coordinates,
    const label celli,
    const label tetFacei,
    const label tetPti
)
:
    ParcelType(mesh, coordinates, celli, tetFacei, tetPti),
    mass0_(0.0),
    Y_(0),
    d0_(0.0),		// E. Weiand - 2020/03/04
    Kevap_(0.0),	// E. Weiand - 2020/03/04
    Dp_(0.0),		// E. Weiand - 2020/03/04
    beta_(1.0),		// E. Weiand - 2020/03/04
    isDenselyPacked_(false),	// E. Weiand - 2020/09/29
    densMult_(1.0),		// E. Weiand - 2020/09/29
    f0_(0.0),
    f0Q_(0.0),
    dtPrev_(0.0),
    dtPrevQ_(0.0),
    isFirstStep_(true),
    isFirstStepQ_(true)
{
    d0_ = this->d();	// E. Weiand - 2020/03/04
}


template<class ParcelType>
inline Foam::ReactingParcel<ParcelType>::ReactingParcel
(
    const polyMesh& mesh,
    const vector& position,
    const label celli
)
:
    ParcelType(mesh, position, celli),
    mass0_(0.0),
    Y_(0),
    d0_(0.0),		// E. Weiand - 2020/03/04
    Kevap_(0.0),	// E. Weiand - 2020/03/04
    Dp_(0.0),		// E. Weiand - 2020/03/04
    beta_(1.0),		// E. Weiand - 2020/03/04
    isDenselyPacked_(false),	// E. Weiand - 2020/09/29
    densMult_(1.0),		// E. Weiand - 2020/09/29
    f0_(0.0),
    f0Q_(0.0),
    dtPrev_(0.0),
    dtPrevQ_(0.0),
    isFirstStep_(true),
    isFirstStepQ_(true)
{
    d0_ = this->d();	// E. Weiand - 2020/03/04
}


template<class ParcelType>
inline Foam::ReactingParcel<ParcelType>::ReactingParcel
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
    const scalar Qdens0,		// E. Weiand - 27/02/2019 - electric particle charge density
    const vector& E0,			// E. Weiand - 30/10/2019 - electric field strength
    const vector& FLorentz0,		// E. Weiand - 30/10/2019 - Lorentz force
    const constantProperties& constProps
)
:
    ParcelType
    (
        mesh,
        coordinates,
        celli,
        tetFacei,
        tetPti,
        typeId,
        nParticle0,
        d0,
        dTarget0,
        U0,
        f0,
        angularMomentum0,
        torque0,
	Qdens0,			// E. Weiand - 27/02/2020 - electric particle charge density
	E0,			// E. Weiand - 30/10/2019 - electric field strength
	FLorentz0,		// E. Weiand - 30/10/2019 - Lorentz force
        constProps
    ),
    mass0_(0.0),
    Y_(Y0),
    d0_(0.0),			// E. Weiand - 2020/03/04 - initial diameter
    Kevap_(0.0),		// E. Weiand - 2020/03/04 - evaporation rate
    Dp_(0.0),			// E. Weiand - 2020/03/04 - particle diffusivity
    beta_(1.0),			// E. Weiand - 2020/03/04 - surface blockage
    isDenselyPacked_(false),	// E. Weiand - 2020/09/29 - surface blockage flag
    densMult_(1.0),		// E. Weiand - 2020/09/29 - density multiplier
    f0_(0.0),
    f0Q_(0.0),
    dtPrev_(0.0),
    dtPrevQ_(0.0),
    isFirstStep_(true),
    isFirstStepQ_(true)
{
    // Set initial parcel mass
    mass0_ = this->mass();

    d0_ = this->d();		// E. Weiand - 2020/03/04 - initial particle diameter
}


// * * * * * * * * * constantProperties Member Functions * * * * * * * * * * //

template<class ParcelType>
inline Foam::scalar
Foam::ReactingParcel<ParcelType>::constantProperties::pMin() const
{
    return pMin_.value();
}


template<class ParcelType>
inline bool
Foam::ReactingParcel<ParcelType>::constantProperties::constantVolume() const
{
    return constantVolume_.value();
}


// * * * * * * * * * * ThermoParcel Member Functions * * * * * * * * * * * * //

template<class ParcelType>
inline Foam::scalar Foam::ReactingParcel<ParcelType>::mass0() const
{
    return mass0_;
}


// E. Weiand - 2020/03/04
template<class ParcelType>
inline Foam::scalar Foam::ReactingParcel<ParcelType>::d0() const
{
    return d0_;
}


// E. Weiand - 2020/03/04
template<class ParcelType>
inline Foam::scalar Foam::ReactingParcel<ParcelType>::Kevap() const
{
    return Kevap_;
}


// E. Weiand - 2020/03/04
template<class ParcelType>
inline Foam::scalar Foam::ReactingParcel<ParcelType>::Dp() const
{
    return Dp_;
}


// E. Weiand - 2020/03/04
template<class ParcelType>
inline Foam::scalar Foam::ReactingParcel<ParcelType>::beta() const
{
    return beta_;
}

// E. Weiand - 2020/09/29
template<class ParcelType>
inline bool Foam::ReactingParcel<ParcelType>::isDenselyPacked() const
{
    return isDenselyPacked_;
}

// E. Weiand - 2020/09/29
template<class ParcelType>
inline Foam::scalar Foam::ReactingParcel<ParcelType>::densMult() const
{
    return densMult_;
}

template<class ParcelType>
inline Foam::scalar Foam::ReactingParcel<ParcelType>::f0() const
{
    return f0_;
}

template<class ParcelType>
inline Foam::scalar Foam::ReactingParcel<ParcelType>::dtPrev() const
{
    return dtPrev_;
}

template<class ParcelType>
inline bool Foam::ReactingParcel<ParcelType>::isFirstStep() const
{
    return isFirstStep_;
}


template<class ParcelType>
inline const Foam::scalarField& Foam::ReactingParcel<ParcelType>::Y() const
{
    return Y_;
}


template<class ParcelType>
inline Foam::scalar& Foam::ReactingParcel<ParcelType>::mass0()
{
    return mass0_;
}


// E. Weiand - 2020/03/04
template<class ParcelType>
inline Foam::scalar& Foam::ReactingParcel<ParcelType>::d0()
{
    return d0_;
}


// E. Weiand - 2020/03/04
template<class ParcelType>
inline Foam::scalar& Foam::ReactingParcel<ParcelType>::Kevap()
{
    return Kevap_;
}


// E. Weiand - 2020/03/04
template<class ParcelType>
inline Foam::scalar& Foam::ReactingParcel<ParcelType>::Dp()
{
    return Dp_;
}


// E. Weiand - 2020/03/04
template<class ParcelType>
inline Foam::scalar& Foam::ReactingParcel<ParcelType>::beta()
{
    return beta_;
}

// E. Weiand - 2020/03/04
template<class ParcelType>
inline bool& Foam::ReactingParcel<ParcelType>::isDenselyPacked()
{
    return isDenselyPacked_;
}

// E. Weiand - 2020/03/04
template<class ParcelType>
inline Foam::scalar& Foam::ReactingParcel<ParcelType>::densMult()
{
    return densMult_;
}

template<class ParcelType>
inline Foam::scalar& Foam::ReactingParcel<ParcelType>::f0()
{
    return f0_;
}

template<class ParcelType>
inline Foam::scalar& Foam::ReactingParcel<ParcelType>::dtPrev()
{
    return dtPrev_;
}


template<class ParcelType>
inline bool& Foam::ReactingParcel<ParcelType>::isFirstStep()
{
    return isFirstStep_;
}


template<class ParcelType>
inline Foam::scalarField& Foam::ReactingParcel<ParcelType>::Y()
{
    return Y_;
}


// ************************************************************************* //
