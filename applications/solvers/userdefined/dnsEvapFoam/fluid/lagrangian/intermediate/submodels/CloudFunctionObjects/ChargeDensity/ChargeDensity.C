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

#include "ChargeDensity.H"

// * * * * * * * * * * * * * Protectd Member Functions * * * * * * * * * * * //

template<class CloudType>
void Foam::ChargeDensity<CloudType>::write()
{
    if (rhoQPtr_.valid())
    {
        rhoQPtr_->write();
    }
    else
    {
        FatalErrorInFunction
            << "rhoQPtr not valid" << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ChargeDensity<CloudType>::ChargeDensity
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName
)
:
    CloudFunctionObject<CloudType>(dict, owner, modelName, typeName),
    rhoQPtr_(nullptr)
{}


template<class CloudType>
Foam::ChargeDensity<CloudType>::ChargeDensity
(
    const ChargeDensity<CloudType>& rq
)
:
    CloudFunctionObject<CloudType>(rq),
    rhoQPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ChargeDensity<CloudType>::~ChargeDensity()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::ChargeDensity<CloudType>::preEvolve()
{
    if (rhoQPtr_.valid())
    {
        rhoQPtr_->primitiveFieldRef() = 0.0;
    }
    else
    {
        const fvMesh& mesh = this->owner().mesh();

        rhoQPtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    this->owner().name() + "RhoQ",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar(dimensionSet(0, -3, 1, 0, 0, 1, 0), 0)
            )
        );
    }
}


template<class CloudType>
void Foam::ChargeDensity<CloudType>::postEvolve()
{
    volScalarField& rhoQ = rhoQPtr_();

    const fvMesh& mesh = this->owner().mesh();

    rhoQ.primitiveFieldRef() /= mesh.time().deltaTValue()*mesh.V();

    CloudFunctionObject<CloudType>::postEvolve();
}


template<class CloudType>
void Foam::ChargeDensity<CloudType>::postMove
(
    parcelType& p,
    const scalar dt,
    const point&,
    bool&
)
{
    volScalarField& rhoQ = rhoQPtr_();

    rhoQ[p.cell()] += dt*p.nParticle()*p.Qdens()*p.volume();
}


// ************************************************************************* //
