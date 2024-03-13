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

\*---------------------------------------------------------------------------*/

#include "NanofluidEvaporation.H"
#include "specie.H"
#include "mathematicalConstants.H"
#include "fundamentalConstants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class CloudType>
Foam::tmp<Foam::scalarField> Foam::NanofluidEvaporation<CloudType>::calcXc
(
    const label celli
) const
{
    scalarField Xc(this->owner().thermo().carrier().Y().size());

    forAll(Xc, i)
    {
        Xc[i] =
            this->owner().thermo().carrier().Y()[i][celli]
           /this->owner().thermo().carrier().Wi(i);
    }

    return Xc/sum(Xc);
}

template<class CloudType>
Foam::scalar Foam::NanofluidEvaporation<CloudType>::Sh
(
    const scalar Re,
    const scalar Sc
) const
{
    // E. Weiand - 2020/03/06: changed to Clift correlation
    //return 2.0 + 0.6*Foam::sqrt(Re)*cbrt(Sc);
    const scalar fRe = (Re < 1.0) ? 1.0 : pow(Re, 0.077);
    return 1.0 + pow(1.0 + Re*Sc, 1.0/3.0) * fRe;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::NanofluidEvaporation<CloudType>::NanofluidEvaporation
(
    const dictionary& dict,
    CloudType& owner
)
:
    PhaseChangeModel<CloudType>(dict, owner, typeName),
    liquids_(owner.thermo().liquids()),
    activeLiquids_(this->coeffDict().lookup("activeLiquids")),
    liqToCarrierMap_(activeLiquids_.size(), -1),
    liqToLiqMap_(activeLiquids_.size(), -1),
    alpha_(readScalar(this->coeffDict().lookup("alpha"))),
    rP_(readScalar(this->coeffDict().lookup("rP"))),
    Yvo_(readScalar(this->coeffDict().lookup("Yvo"))),
    etal_(readScalar(this->coeffDict().lookup("etaLiquid")))
{
    if (activeLiquids_.size() == 0)
    {
        WarningInFunction
            << "Evaporation model selected, but no active liquids defined"
            << nl << endl;
    }
    else
    {
        Info<< "Participating liquid species:" << endl;

        // Determine mapping between liquid and carrier phase species
        forAll(activeLiquids_, i)
        {
            Info<< "    " << activeLiquids_[i] << endl;
            liqToCarrierMap_[i] =
                owner.composition().carrierId(activeLiquids_[i]);
        }

        // Determine mapping between model active liquids and global liquids
        const label idLiquid = owner.composition().idLiquid();
        forAll(activeLiquids_, i)
        {
            liqToLiqMap_[i] =
                owner.composition().localId(idLiquid, activeLiquids_[i]);
        }
    }
}


template<class CloudType>
Foam::NanofluidEvaporation<CloudType>::NanofluidEvaporation
(
    const NanofluidEvaporation<CloudType>& pcm
)
:
    PhaseChangeModel<CloudType>(pcm),
    liquids_(pcm.owner().thermo().liquids()),
    activeLiquids_(pcm.activeLiquids_),
    liqToCarrierMap_(pcm.liqToCarrierMap_),
    liqToLiqMap_(pcm.liqToLiqMap_),
    alpha_(pcm.alpha_),
    rP_(pcm.rP_),
    Yvo_(pcm.Yvo_),
    etal_(pcm.etal_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::NanofluidEvaporation<CloudType>::~NanofluidEvaporation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// M.P. Sitte, 03/12/2015
// E. Weiand, 05/01/2020
template<class CloudType>
void Foam::NanofluidEvaporation<CloudType>::calculateAS
(
    const scalar dt,
    const label cellI,
    const scalar Re,   // Reynolds @ ref. cond.
    const scalar Pr,   // Prandtl  @ ref. cond.
    const scalar d,    // droplet diameter
    const scalar rhos, // density at ref. cond.   // M.P.S.
    const scalar nu,   // nu @ ref. cond.
    const scalar T,    // droplet temp. 
    const scalar Ts,   // temperature @ ref. cond.
    const scalar pc,   // pressure in cell
    const scalar Tc,   // temperatrue in cell
    const scalarField& Yl,
    const scalar d0,   	// E. Weiand - 2020/03/04
    scalar& BM,	// E. Weiand - 2020/03/30 - global BM test
    scalarField& dMassPC,
    scalar& fmixS,    	// M.P.S.
    scalar& grFmixS,  	// M.P.S.
    scalar& DS,       	// M.P.S.
    scalar& Kevap,    	// E. Weiand - 2020/03/04
    scalar& Dp,		// E. Weiand - 2020/03/04
    scalar& beta,	// E. Weiand - 2020/03/04
    scalar& tauV,	// E. Weiand - 2020/03/06
    bool& isDenselyPacked,	// E. Weiand - 2020/09/21
    scalar& f0,
    scalar& dtPrev,
    bool& isFirstStep
) const
{

    // liquid volume fraction
    const scalarField X(liquids_.X(Yl));


    fmixS= 0.0;
    grFmixS= 0.0;
    DS= 0.0;
 

    // immediately evaporate mass that has reached critical condition
    // liquids_.Tc(X) ... critical temperature
    if ((liquids_.Tc(X) - T) < SMALL)
    {
        std::cerr<<  "WARNING ! T > Tcritical ! " <<  "\n";
    }


    // construct carrier phase species volume fractions for cell, cellI
    const scalarField XcMix(calcXc(cellI));

    // species mass fractions in cell
    scalarField YcMix(XcMix);
    forAll(YcMix,i)
    {
        YcMix[i]= this->owner().thermo().carrier().Y()[i][cellI];
    }


    // carrier mean molecular weight, excluding vapour
    scalar Wcarr= 0.0;
    scalarField Xcarr(XcMix);
    forAll(activeLiquids_, i)
    {
        const label gid = liqToCarrierMap_[i];
        Xcarr[gid]= 0.0;
    }
    Xcarr /= sum(Xcarr);
    forAll(Xcarr, i)
    {
        Wcarr += Xcarr[i] * this->owner().thermo().carrier().Wi(i);
    }


    // calculate mass transfer of each specie in liquid
    forAll(activeLiquids_, i)
    {
        const label gid = liqToCarrierMap_[i];
        const label lid = liqToLiqMap_[i];

        // boiling temperature at cell pressure for liquid species lid [K]
        const scalar TBoil = liquids_.properties()[lid].pvInvert(pc);

        // limit droplet temperature to boiling/critical temperature
        const scalar Td = min(T, 0.9999*TBoil);
  
        // saturation pressure for liquid species lid [Pa]
        //const scalar pSat = liquids_.properties()[lid].pv(pc, Td);

	// E. Weiand - 2020/03/06
	// saturation pressure from Antoine equation for ethanol
	const scalar Apv = 5.37229;
	const scalar Bpv = 1670.409;
	const scalar Cpv = -40.191;
	scalar pSat = pow(10, Apv - Bpv/(Td+Cpv)) * 1e5;
	pSat = max(pSat, 0.0);

        // surface equilibrium mole fraction
        const scalar Xs= X[lid] * pSat / pc;

        // surface equilibrium mass fraction
        const scalar Ys= Xs / (Xs + (1.0-Xs)*(Wcarr/liquids_.properties()[lid].W()));
 
        // carrier phase vapour mass fraction
        const scalar Yc = YcMix[gid];

      // vapour diffusivity [m2/s]
//	const scalar Dab = this->owner().constProps().D0();		// E. Weiand - 2020/03/04: constant diffusivity
	const scalar Dfactor = this->owner().constProps().D0();
        const scalar Dab = Dfactor*(8.967e-11*pow(Ts, 2.0) + 2.539e-08*Ts - 3.108e-06);
//        const scalar Dab = liquids_.properties()[lid].D(pc, Ts, Wcarr);  // M.P.Sitte. corrected: pc instead of ps

        // Schmidt number
        const scalar Sc = nu/Dab;

        // Sherwood number
        const scalar Sh = this->Sh(Re, Sc);

        // evaporation
        // using notation as in Miller et al. 1998. // M.P. Sitte
        BM = max((Ys - Yc), 0.0) / max(1.0 - Ys, SMALL);
        BM = max(BM, 1.0e-10);

        // Abramzoan & Sirignano: correcttion for Stefan flow
        const scalar FM = pow(1.0 + BM,0.7) * log(1.0 + BM) / BM;
        const scalar ShStar = 2.0 + (Sh-2.0)/FM;

	/*Info << "nu = " << nu << endl;
	Info << "FM = " << FM << endl;
	Info << "ShStar = " << ShStar << endl;
	Info << "diam = " << d << endl;
	Info << "rhos = " << rhos << endl;	
	Info << "Yc = " << Yc << endl;
	Pout << "Ys = " << Ys << endl;*/

	///////////////////////////////////////////////////////////////////////////////////////
	// Nanofuel correction ////////////////////////////////////////////////////////////////
	// E. Weiand //////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////

	// nanoparticle radius [m]
	const scalar rp = rP_;

	// sine of contact angle squared
	const scalar alpha = alpha_;

	// initial nanoparticle mass fraction per droplet [-]
	const scalar Yvo = Yvo_;

	// liquid density
	const scalar rhol = this->owner().composition().liquids().properties()[lid].rho(pc, Td);

	// pi
	const scalar pi = constant::mathematical::pi;

	// Boltzmann constant
	const scalar kb = physicoChemical::k.value();

	// liquid dynamic viscosity [Pa s]
	const scalar eta = etal_;

	// particle diffusivity in carrier phase
	Dp = kb * T / (6.0 * pi * eta * rp);

	// evaporation rate constant [m^2/s]
	const scalar K = 4.0*ShStar*Dab*rhos/rhol*log(1.0 + BM);

	// droplet Péclet number [-]
	const scalar Pe = K / (8.0 * Dp);

	// surface enrichment [-]
	/* piece-wise polynomials for low and high Péclet regimes
		evaluated from E = e^{Pe/2} / (3 int_0^1 R^2 exp(Pe R^2 / 2)dR; */
	scalar E;
	if (Pe > 10.0)
	    E = 0.3964 + 0.3316 * Pe + 2.102e-5 * pow(Pe, 2.0) - 1.392e-7 * pow(Pe, 3.0) +
		2.147e-10 * pow(Pe, 4.0);
	else
	    E = 0.9987 + 0.2008 * Pe + 0.01185 * pow(Pe, 2.0) - 0.0005069 * pow(Pe, 3.0) +
		8.382e-6 * pow(Pe, 4.0);

	// ratio of effective to actual evaporation rate constant [-]
	const scalar betaPrev = beta;
	beta = max(1.0 - 1.5 * alpha * Yvo * E * pow(d0 / d, 3.0), 0.0);

	// surface volume fraction
	const scalar Yvs = Yvo * pow(d0 / d, 3.0) * E;

	/*Info << "Yvs = " << Yvs << endl;
	Info << "d = " << d << endl;
	Info << "E = " << E << endl;
	Info << "beta = " << beta << endl;
	Info << "Pe = " << Pe << endl;
	Info << "K = " << K << endl;
	Info << "T = " << T << endl;*/

	// end of first stage evaporation => supress further reduction due to surface blockage
	if (Yvs > 0.6)
	{
	    /* Assumptions for constant beta:
		- no pore diffusion/Knudsen effect
		- no influence of nanoparticles on droplet heat transfer (very rough assumption)
		- film evaporation unaffected by shell-formation
	    */
	    beta = betaPrev;

	    // maximum particle packing reached
	    isDenselyPacked = true;

	    // Future refined modeling here. Consider heat transfer, pore diffusion etc.
	}

	///////////////////////////////////////////////////////////////////////////////////////


        // mass transfer [kg]
        dMassPC[lid] += beta*pi*rhos*d*Dab*ShStar*log(1.0 + BM)*dt;

/*        Info << "dtPrev = " << dtPrev << endl;
	Info << "f0 = " << f0 << endl;
	scalar f1 = beta*pi*rhos*d*Dab*ShStar*log(1.0 + BM);
	scalar dy2 = 0.0;
	if (isFirstStep == true)
	{
	    // Euler method
	    dy2 = f1*dt;
	    isFirstStep = false;
	}
	else
	{
	    // 2nd order Adams-Bashforth method
	    dy2 = 0.5*(dt + 2*dtPrev)*f1 - 0.5*dt*f0;
	}
        f0 = f1;
	dtPrev = dt;
	dMassPC[lid] += dy2;*/

	// evaporation time scale
	tauV = rhol * pow(d, 2.0) / (6 * Dab * ShStar * rhos * log(1.0 + BM));

        // Surface mixture fraction // M.P. Sitte, 11/05/2016
        fmixS += max(Yc, Ys);

        grFmixS += grFmixS + abs( -ShStar * log(1.0 + BM) * (1.0-Ys)/d  );

        // Surface diffusivity // M.P. Sitte, 04/10/2017
        DS += Dab;

	// E. Weiand - 2020/03/04
	// Evaporation rate
	Kevap = K;
    }
}

template<class CloudType>
Foam::scalar Foam::NanofluidEvaporation<CloudType>::dh
(
    const label idc,
    const label idl,
    const scalar p,
    const scalar T
) const
{
    scalar dh = 0;

    typedef PhaseChangeModel<CloudType> parent;
    switch (parent::enthalpyTransfer_)
    {
        case (parent::etLatentHeat):
        {
            dh = liquids_.properties()[idl].hl(p, T);
            break;
        }
        case (parent::etEnthalpyDifference):
        {
            scalar hc = this->owner().composition().carrier().Ha(idc, p, T);
            scalar hp = liquids_.properties()[idl].h(p, T);

            dh = hc - hp;
            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Unknown enthalpyTransfer type" << abort(FatalError);
        }
    }

    return dh;
}


template<class CloudType>
Foam::scalar Foam::NanofluidEvaporation<CloudType>::Tvap
(
    const scalarField& X
) const
{
    return liquids_.Tpt(X);
}


template<class CloudType>
Foam::scalar Foam::NanofluidEvaporation<CloudType>::TMax
(
    const scalar p,
    const scalarField& X
) const
{
    return liquids_.pvInvert(p, X);
}


// ************************************************************************* //
