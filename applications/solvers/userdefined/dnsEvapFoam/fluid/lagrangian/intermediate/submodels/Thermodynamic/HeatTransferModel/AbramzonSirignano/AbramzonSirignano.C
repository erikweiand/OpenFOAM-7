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

#include "AbramzonSirignano.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::AbramzonSirignano<CloudType>::AbramzonSirignano
(
    const dictionary& dict,
    CloudType& cloud
)
:
    HeatTransferModel<CloudType>(dict, cloud, typeName)
{}


template<class CloudType>
Foam::AbramzonSirignano<CloudType>::AbramzonSirignano(const AbramzonSirignano<CloudType>& htm)
:
    HeatTransferModel<CloudType>(htm)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::AbramzonSirignano<CloudType>::~AbramzonSirignano()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::AbramzonSirignano<CloudType>::Nu
(
    const scalar Re,
    const scalar Pr
) const
{
    scalar fc = 1.0;
    if (Re > 1.0 && Re <= 400.0)
 	fc = pow(Re, 0.077);
    else if (Re > 400.0)
	FatalErrorInFunction << "Re_p > 400 for Nusselt approx." << abort(FatalError);

    const scalar Nu0 = 1.0 + pow(1.0 + Re * Pr, 1.0/3.0) * fc;

    // to be implemented:
    const scalar Bt = 1.0;

    const scalar Ft = pow(1.0 + Bt, 0.7) * log(1.0 + Bt) / Bt;
    const scalar NuPrime = 2.0 + (Nu0 - 2.0) / Ft;
    return log(1.0 + Bt) / Bt * NuPrime;
}


// ************************************************************************* //
