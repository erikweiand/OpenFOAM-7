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

#include "ReactingMultiphaseParcel.H"
#include "mathematicalConstants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class ParcelType>
const Foam::label Foam::ReactingMultiphaseParcel<ParcelType>::GAS(0);

template<class ParcelType>
const Foam::label Foam::ReactingMultiphaseParcel<ParcelType>::LIQ(1);

template<class ParcelType>
const Foam::label Foam::ReactingMultiphaseParcel<ParcelType>::SLD(2);


// * * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * //

template<class ParcelType>
template<class TrackCloudType>
Foam::scalar Foam::ReactingMultiphaseParcel<ParcelType>::CpEff
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar p,
    const scalar T,
    const label idG,
    const label idL,
    const label idS
) const
{
    return
        this->Y_[GAS]*cloud.composition().Cp(idG, YGas_, p, T)
      + this->Y_[LIQ]*cloud.composition().Cp(idL, YLiquid_, p, T)
      + this->Y_[SLD]*cloud.composition().Cp(idS, YSolid_, p, T);
}


template<class ParcelType>
template<class TrackCloudType>
Foam::scalar Foam::ReactingMultiphaseParcel<ParcelType>::HsEff
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar p,
    const scalar T,
    const label idG,
    const label idL,
    const label idS
) const
{
    return
        this->Y_[GAS]*cloud.composition().Hs(idG, YGas_, p, T)
      + this->Y_[LIQ]*cloud.composition().Hs(idL, YLiquid_, p, T)
      + this->Y_[SLD]*cloud.composition().Hs(idS, YSolid_, p, T);
}


template<class ParcelType>
template<class TrackCloudType>
Foam::scalar Foam::ReactingMultiphaseParcel<ParcelType>::LEff
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar p,
    const scalar T,
    const label idG,
    const label idL,
    const label idS
) const
{
    return
        this->Y_[GAS]*cloud.composition().L(idG, YGas_, p, T)
      + this->Y_[LIQ]*cloud.composition().L(idL, YLiquid_, p, T)
      + this->Y_[SLD]*cloud.composition().L(idS, YSolid_, p, T);
}


template<class ParcelType>
Foam::scalar Foam::ReactingMultiphaseParcel<ParcelType>::updateMassFractions
(
    const scalar mass0,
    const scalarField& dMassGas,
    const scalarField& dMassLiquid,
    const scalarField& dMassSolid
)
{
    scalarField& YMix = this->Y_;

    scalar massGas =
        this->updateMassFraction(mass0*YMix[GAS], dMassGas, YGas_);
    scalar massLiquid =
        this->updateMassFraction(mass0*YMix[LIQ], dMassLiquid, YLiquid_);
    scalar massSolid =
        this->updateMassFraction(mass0*YMix[SLD], dMassSolid, YSolid_);

    scalar massNew = max(massGas + massLiquid + massSolid, rootVSmall);

    YMix[GAS] = massGas/massNew;
    YMix[LIQ] = massLiquid/massNew;
    YMix[SLD] = 1.0 - YMix[GAS] - YMix[LIQ];

    return massNew;
}


// * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * * //

template<class ParcelType>
template<class TrackCloudType>
void Foam::ReactingMultiphaseParcel<ParcelType>::setCellValues
(
    TrackCloudType& cloud,
    trackingData& td
)
{
    ParcelType::setCellValues(cloud, td);
}


template<class ParcelType>
template<class TrackCloudType>
void Foam::ReactingMultiphaseParcel<ParcelType>::cellValueSourceCorrection
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt
)
{
    // Re-use correction from reacting parcel
    ParcelType::cellValueSourceCorrection(cloud, td, dt);
}


// <calc> modified by M.P.Sitte
template<class ParcelType>
template<class TrackCloudType>
void Foam::ReactingMultiphaseParcel<ParcelType>::calc
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt
)
{
    typedef typename TrackCloudType::reactingCloudType reactingCloudType;
    const CompositionModel<reactingCloudType>& composition =
        cloud.composition();

    // Define local properties at beginning of timestep
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


    const label idG = composition.idGas();
    const label idL = composition.idLiquid();
    const label idS = composition.idSolid();

    // E. Weiand - 06/01/2020
    const label cellI = this->cell();

    // fuel index in liquids. Assumes single component fuel  // M.P.Sitte
    label lid= 100;
    forAll(YLiquid_, i)
    {
	// E. Weiand - 06/02/2020
        //lid= composition.globalIds(idL)[i];
	lid = i;
    }

    // ag813 - Correct for negative pressure
    scalar pGas = td.pc();
    if (pGas < 0.0)
    {
        pGas  = 101325;
    }


    // Compute TBoil        // M.P.Sitte, 2-17-06-09
    scalar Tboil = composition.liquids().properties()[lid].pvInvert(pGas);

    const scalar massLoc = this->mass();

    // Limit temperature    // M.P.Sitte, 2017-06-09
    if (this->T_ > Tboil) 
    {
        std::cerr << "WARNING 1:: T>TBoil.   T=" << this->T_ << ".  Limiting T" << "\n";
        this->T_= min(this->T_, 0.9999*Tboil);
    }
    else if (this->T_ < 200.0)
    {
        std::cerr << "WARNING 1:: T<200 K.   T=" << this->T_ << ".  Limiting T" << "\n";
        this->T_= max(this->T_, 200.0);
    }

    const bool isDenselyPacked0 = this->isDenselyPacked();
    scalar& densMult = this->densMult();
    
    // E. Weiand - 2020/03/02: constant droplet density
    this->rho_ = max(cloud.constProps().rho0()*densMult,1e-9);
    
    // update d after limiting T   // M.P.Sitte, 2017-06-09
    // - NOTE: important to use massLoc instead of this->mass() - otherwise d=cst, mass decreases
    this->d_ = (isDenselyPacked0) ? this->dPacked() : cbrt(massLoc/this->rho_*6.0/pi);

   
    const scalar np0 = this->nParticle_;
    const scalar d0 = this->d_;
    const vector& U0 = this->U_;
    const scalar T0 = this->T_;
    const scalar mass0 = this->mass();

    scalar dCurr = (isDenselyPacked0) ? this->dPacked() : d0;
    
    scalar pc = td.pc();
    //ag813 - correction for negative pressure
    if (pc < 0.0)
    {
        pc  = 101325;
    }

    const scalarField& YMix = this->Y_;

    // Calc surface values
    scalar Ts, rhos, mus, Prs, kappas, Cps;

    // Mean molecular mass at ref. cond.
    scalar Ws;

    // Composition vector at ref. cond. YYs // M.P.Sitte
    scalarField YYs(composition.carrier().species().size(), 0.0);

    // Spalding transfer number for mass
    scalar BM;

    // M.P.Sitte 15/03/2016: call <calcSurfaceValuesAS>
    this->calcSurfaceValuesAS
    (   cloud,
	td, 
        cellI, 
        T0, 
        idL, 
        YLiquid_, 
        Ts, 
        YYs,
        Ws,
        rhos, 
        mus, 
        Prs, 
        kappas,
        BM
    );

    // M.P. Sitte: call <correctSurfaceValuesAS>
    this->correctSurfaceValuesAS
    (   cloud,
	td,  
        Ts, 
        YYs, 
        Ws, 
        rhos, 
        mus, 
        Prs, 
        kappas,
        Cps
    );

    //Pout <<  "parcel " << this->origId() << ": isDenselyPacked = " << isDenselyPacked0 << endl;

    // Reynolds number @ ref. cond.
    scalar Res = this->Re(td.rhoc(), U0, td.Uc(), dCurr, mus);  // rhoc according to Abramzon & Sazhin


    // Sources
    //~~~~~~~~

    // Explicit momentum source for particle
    vector Su = Zero;

    // Linearised momentum source coefficient
    scalar Spu = 0.0;

    // Momentum transfer from the particle to the carrier phase
    vector dUTrans = Zero;

    // Explicit enthalpy source for particle
    scalar Sh = 0.0;

    // Linearised enthalpy source coefficient
    scalar Sph = 0.0;

    // Sensible enthalpy transfer from the particle to the carrier phase
    scalar dhsTrans = 0.0;


    // 1. Compute models that contribute to mass transfer - U, T held constant
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Phase change in liquid phase
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Mass transfer due to phase change
    scalarField dMassPC(YLiquid_.size(), 0.0);

    // Phase change leading to source of species variance ... // M.P. Sitte, 12/10/2015
    scalar fmixS = 0.0;

    // Phase change leading to source of species variance ... // M.P. Sitte, 12/10/2015
    scalar grFmixS = 0.0;

    // Phase change leading to source of species variance ... // M.P. Sitte, 12/10/2015
    scalar DS = 0.0;


    // Molar flux of species emitted from the particle (kmol/m^2/s)
    scalar Ne = 0.0;

    // Sum Ni*Cpi*Wi of emission species
    scalar NCpW = 0.0;

    // Surface concentrations of emitted species
    scalarField Cs(composition.carrier().species().size(), 0.0);

    // M.P. Sitte, 11/10/2015: use <calcPhaseChangeAS>
    // Calc mass and enthalpy transfer due to phase change
    // Abramzon and Sirignano model
 
    if ((YMix[LIQ] > SMALL) && (np0*mass0*YMix[LIQ] > cloud.constProps().minParcelMass()))
    {
        this->calcPhaseChangeAS
        (
	    cloud,
            td,
            dt,
            cellI,
            Res,
            Prs,
            Ts,
            rhos, // density at ref. cond. // M.P.Sitte, 11/10/2015
            mus/rhos,
            dCurr,
            T0,
            mass0,
            idL,
            YMix[LIQ],
            YLiquid_,
	    BM,
            dMassPC,
            fmixS,       // added ... M.P. Sitte, 12/10/2015
            grFmixS,     // added ... M.P. Sitte, 11/10/2015
            DS,          // added ... M.P. Sitte, 12/10/2015
            Sh,  // enthalpy change of droplet in [kg * J/kg / s]
            Ne,
            NCpW,
            Cs
        );
    }

    const label isDenselyPacked1 = this->isDenselyPacked();

    if ((isDenselyPacked0 == false) && (isDenselyPacked1 == true))
	this->dPacked_ = this->d();

    // Mass transfer due to devolatilisation
    scalarField dMassDV(YGas_.size(), 0.0);

    // Change in carrier phase composition due to surface reactions
    scalarField dMassSRGas(YGas_.size(), 0.0);
    scalarField dMassSRLiquid(YLiquid_.size(), 0.0);
    scalarField dMassSRSolid(YSolid_.size(), 0.0);
    scalarField dMassSRCarrier(composition.carrier().species().size(), 0.0);


    // 2. Update the parcel properties due to change in mass
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    scalarField dMassGas(dMassDV + dMassSRGas);
    scalarField dMassLiquid(dMassPC + dMassSRLiquid);
    scalarField dMassSolid(dMassSRSolid);

    scalar mass1 =
        updateMassFractions(mass0, dMassGas, dMassLiquid, dMassSolid);

    // E. Weiand - 2020/09/29
    // update density multiplier if densely packed
    if (isDenselyPacked1)
	densMult = mass1 / (pi/6.0*cloud.constProps().rho0()*pow(this->dPacked(),3.0));

    this->Cp_ = CpEff(cloud, td, pc, T0, idG, idL, idS);

/*
    // Update particle density or diameter
    if (cloud.constProps().constantVolume())
    {
        this->rho_ = mass1/this->volume();
    }
    else
    {
        this->rho_ = cloud.constProps().rho0();		// E. Weiand - 2020/03/07: constant properties
        this->d_ = cbrt(mass1/this->rho_*6.0/pi);
    }
*/

    // update rho and d    // M.P.Sitte, 2017-06-09
//    this->rho_ = composition.liquids().properties()[lid].rho(pc, min(this->T_, Tboil));

    // E. Weiand - 2020/09/22 (if condition)
    if ((isDenselyPacked1 == false) && (YMix[LIQ] > SMALL) && (np0*mass1*YMix[LIQ] > cloud.constProps().minParcelMass()))	// very crude fix for now
    {
        // Calculate new particle temperature
        this->T_ =
            this->calcHeatTransferAS // M.P.Sitte
            (
	        cloud,
                td,
                dt,
                cellI,
                idL,
                Res,    // Re @ ref. cond.
                Prs,    // Pr @ ref. cond.
                Ts,     // ref. temperature           // new
                rhos,   // density @ ref. cond.     // new
                mus/rhos,// viscosity @ ref. cond.     // new
                Cps,    // gas Cp @ ref. cond.        // new
                kappas, // kappa @ ref. cond.
                BM,     // mass transf. number        // new
                dMassPC,// change in mass             // new
                Sh,     // evaporation enthalpy in [J/s]
                dhsTrans,
                Sph
            );
    }

    // E. Weiand - 2020/03/02: changed to constant density
    this->rho_ = max(cloud.constProps().rho0()*densMult, 1e-9);
    this->d_ = (isDenselyPacked1) ? this->dPacked() : cbrt(mass1/this->rho_*6.0/pi);

    // Remove the particle when mass falls below minimum threshold
    if (np0*mass1*YMix[LIQ] < cloud.constProps().minParcelMass())
    {
        //td.keepParticle = false;	// E. Weiand - 2020/09/29: commented to keep solid particles

        if (cloud.solution().coupled())
        {
            scalar dm = np0*mass0*YMix[LIQ];

            // Absorb parcel into carrier phase
            /*forAll(YGas_, i)
            {
                label gid = composition.localToCarrierId(GAS, i);
                cloud.rhoTrans(gid)[this->cell()] += dm*YMix[GAS]*YGas_[i];
            }*/
            forAll(YLiquid_, i)
            {
                label gid = composition.localToCarrierId(LIQ, i);
                cloud.rhoTrans(gid)[this->cell()] += dm*YMix[LIQ]*YLiquid_[i];
            }

            // Note: Here no code was added ... M.P. Sitte, 11/10/2015
            //   i.e. mixture fraction variance source do to absorbed particles is neglected

/*
            // No mapping between solid components and carrier phase
            forAll(YSolid_, i)
            {
                label gid = composition.localToCarrierId(SLD, i);
                cloud.rhoTrans(gid)[this->cell()] += dm*YMix[SLD]*YSolid_[i];
            }
*/

            cloud.UTrans()[this->cell()] += dm*U0;

            cloud.hsTrans()[this->cell()] +=
                dm*HsEff(cloud, td, pc, T0, idG, idL, idS);

	    // M.P.Sitte, 23/11/2016
            scalar hL =  composition.liquids().properties()[lid].h(pc,  min(this->T_, Tboil) );
            cloud.hTrans()[cellI] += dm * hL;

            cloud.phaseChange().addToPhaseChangeMass(np0*mass1*YMix[LIQ]);
        }

        //return;	// E. Weiand - 2020/10/01: removed to ensure continuous tracking
    }


    // Correct surface values due to emitted species
    //this->correctSurfaceValues(td, cellI, Ts, Cs, rhos, mus, Prs, kappas);
    // M.P. Sitte: call <correctSurfaceValuesAS>
    this->correctSurfaceValuesAS(cloud, td, Ts, YYs, Ws, rhos, mus, Prs, kappas, Cps);

    Res = this->Re(td.rhoc(), U0, td.Uc(), dCurr, mus);  // rhoc according to Abramzon & Sazhin
    //Res = this->Re(U0, this->d_, rhos, mus);



    // 3. Compute heat- and momentum transfers
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Heat transfer
    // ~~~~~~~~~~~~~

    this->Cp_ = CpEff(cloud, td, pc, this->T_, idG, idL, idS);


    // Motion
    // ~~~~~~

    // E. Weiand - 26/11/2019
    // Calculate new particle Reynolds number (output var)
    const scalar Re = this->Re(td);
    this->particleRe_ = Re;

    // E. Weiand - 26/11/2019
    // Update electric field & Lorentz force
    this->E_ = this->calcElField(cloud, td);
    const scalar Q = this->Qdens_ * pi/6.0 * pow(this->d0_, 3.0);	// charge contained in nanoparticles
    this->FLorentz_ = this->E_ * Q;


    Pout << "Res = " << Res << "; mus = " << mus << endl;

    // Calculate new particle velocity
    this->U_ =
        this->calcVelocity(cloud, td, dt, Res, mus, mass1, Su, dUTrans, Spu);

    // update rho and d     // M.P.Sitte, 2017-06-09
//    this->rho_ = composition.liquids().properties()[lid].rho(pc, min(this->T_, Tboil)); 

    // E. Weiand - 2020/03/02: constant droplet density
    /*this->rho_ = cloud.constProps().rho0();
    this->d_ = cbrt(mass1/this->rho_*6.0/pi);*/


    // 4. Accumulate carrier phase source terms
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (cloud.solution().coupled())
    {
        // Transfer mass lost to carrier mass, momentum and enthalpy sources
        forAll(YGas_, i)
        {
            scalar dm = np0*dMassGas[i];
            label gid = composition.localToCarrierId(GAS, i);
            scalar hs = composition.carrier().Hs(gid, pc, T0);	// E. Weiand - 28/02/2020: T0 correct? 
            cloud.rhoTrans(gid)[this->cell()] += dm;
            cloud.UTrans()[this->cell()] += dm*U0;
            cloud.hsTrans()[this->cell()] += dm*hs;
        }

        scalar dmtot = 0.0;  // M.P.Sitte

        forAll(YLiquid_, i)
        {
            scalar dm = np0*dMassLiquid[i];	// dm > 0 for evap
            label gid = composition.localToCarrierId(LIQ, i);
            scalar hs = composition.carrier().Hs(gid, pc, T0);	// E. Weiand - 28/02/2020: T0 correct? 
            cloud.rhoTrans(gid)[this->cell()] += dm;
            cloud.UTrans()[this->cell()] += dm*U0;
            

	    // E. Weiand, 28/02/2020: required when making use of hsTrans
/*           // M.P. Sitte, 22/11/2016 ... commented, replaced, see below		*/
            cloud.hsTrans()[cellI] += dm*hs;

            // M.P.Sitte, 23/11/2016
            scalar hL =  composition.liquids().properties()[lid].h(pc,  min(this->T_, Tboil));
            //scalar Lv =  composition.liquids().properties()[lid].hl(pc, min(this->T_, Tboil*0.999));

	    cloud.hTrans()[cellI] += dm * hL;

            dmtot += dm; // M.P.Sitte
        }


        // M.P.Sitte, 23/11/2016: hTrans Heat Transfer source
        //   dhsTrans= dt*(-md_dot)*CpV/BT *(Tres.average() - Tc) 


        // M.P. Sitte, 15/10/2015
        // Averages of mixture fraction, its gradient and diffusivity @ droplet surface
        // Average calculated by weighting
        // Weighting based on evaporated mass
        cloud.FsNp()()[cellI]  += fmixS   * dmtot;  //  np0*(this->areaP());
        cloud.GrFsNp()()[cellI]+= grFmixS * dmtot;  //  np0*(this->areaP()); 
        cloud.DsNp()()[cellI]  += DS      * dmtot;  //  np0*(this->areaP());        
        cloud.rdNp()()[cellI]  += (d0/2.0)* dmtot;
        cloud.ndNp()()[cellI]  += np0     * dmtot;

        cloud.Np()()[cellI]    += dmtot;            // np0*(this->areaP());


        cloud.VolNpDt()()[cellI] += np0 * pow3(d0)*pi/6.0 * dmtot * dt; // use avg diameter
        cloud.NpDt()()[cellI]    += dmtot *dt;

        cloud.massd()()[cellI]   = np0 * mass1;


        forAll(dMassSRCarrier, i)
        {
            scalar dm = np0*dMassSRCarrier[i];
            scalar hs = composition.carrier().Hs(i, pc, T0);	// E. Weiand - 28/02/2020: T0 correct? 
            cloud.rhoTrans(i)[this->cell()] += dm;
            cloud.UTrans()[this->cell()] += dm*U0;
            cloud.hsTrans()[this->cell()] += dm*hs;
        }

        // Update momentum transfer
        cloud.UTrans()[this->cell()] += np0*dUTrans;
        cloud.UCoeff()[this->cell()] += np0*Spu;

        // Update sensible enthalpy transfer
        cloud.hsTrans()[this->cell()] += np0*dhsTrans;
        cloud.hsCoeff()[this->cell()] += np0*Sph;

        // Update radiation fields
        if (cloud.radiation())
        {
            const scalar ap = this->areaP();
            const scalar T4 = pow4(T0);
            cloud.radAreaP()[this->cell()] += dt*np0*ap;
            cloud.radT4()[this->cell()] += dt*np0*T4;
            cloud.radAreaPT4()[this->cell()] += dt*np0*ap*T4;
        }
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class ParcelType>
template<class TrackCloudType>
void Foam::ReactingMultiphaseParcel<ParcelType>::calcDevolatilisation
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt,
    const scalar age,
    const scalar Ts,
    const scalar d,
    const scalar T,
    const scalar mass,
    const scalar mass0,
    const scalarField& YGasEff,
    const scalarField& YLiquidEff,
    const scalarField& YSolidEff,
    label& canCombust,
    scalarField& dMassDV,
    scalar& Sh,
    scalar& N,
    scalar& NCpW,
    scalarField& Cs
) const
{
    // Check that model is active
    if (!cloud.devolatilisation().active())
    {
        if (canCombust != -1)
        {
            canCombust = 1;
        }
        return;
    }

    // Initialise demand-driven constants
    (void)cloud.constProps().TDevol();
    (void)cloud.constProps().LDevol();

    // Check that the parcel temperature is within necessary limits for
    // devolatilisation to occur
    if (T < cloud.constProps().TDevol() || canCombust == -1)
    {
        return;
    }

    typedef typename TrackCloudType::reactingCloudType reactingCloudType;
    const CompositionModel<reactingCloudType>& composition =
        cloud.composition();


    // Total mass of volatiles evolved
    cloud.devolatilisation().calculate
    (
        dt,
        age,
        mass0,
        mass,
        T,
        YGasEff,
        YLiquidEff,
        YSolidEff,
        canCombust,
        dMassDV
    );

    scalar dMassTot = sum(dMassDV);

    cloud.devolatilisation().addToDevolatilisationMass
    (
        this->nParticle_*dMassTot
    );

    Sh -= dMassTot*cloud.constProps().LDevol()/dt;

    // Update molar emissions
    if (cloud.heatTransfer().BirdCorrection())
    {
        // Molar average molecular weight of carrier mix
        const scalar Wc = max(small, td.rhoc()*RR*td.Tc()/td.pc());

        // Note: hardcoded gaseous diffusivities for now
        // TODO: add to carrier thermo
        const scalar beta = sqr(cbrt(15.0) + cbrt(15.0));

        forAll(dMassDV, i)
        {
            const label id = composition.localToCarrierId(GAS, i);
            const scalar Cp = composition.carrier().Cp(id, td.pc(), Ts);
            const scalar W = composition.carrier().Wi(id);
            const scalar Ni = dMassDV[i]/(this->areaS(d)*dt*W);

            // Dab calc'd using API vapour mass diffusivity function
            const scalar Dab =
                3.6059e-3*(pow(1.8*Ts, 1.75))
               *sqrt(1.0/W + 1.0/Wc)
               /(td.pc()*beta);

            N += Ni;
            NCpW += Ni*Cp*W;
            Cs[id] += Ni*d/(2.0*Dab);
        }
    }
}


template<class ParcelType>
template<class TrackCloudType>
void Foam::ReactingMultiphaseParcel<ParcelType>::calcSurfaceReactions
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt,
    const scalar d,
    const scalar T,
    const scalar mass,
    const label canCombust,
    const scalar N,
    const scalarField& YMix,
    const scalarField& YGas,
    const scalarField& YLiquid,
    const scalarField& YSolid,
    scalarField& dMassSRGas,
    scalarField& dMassSRLiquid,
    scalarField& dMassSRSolid,
    scalarField& dMassSRCarrier,
    scalar& Sh,
    scalar& dhsTrans
) const
{
    // Check that model is active
    if (!cloud.surfaceReaction().active())
    {
        return;
    }

    // Initialise demand-driven constants
    (void)cloud.constProps().hRetentionCoeff();
    (void)cloud.constProps().TMax();

    // Check that model is active
    if (canCombust != 1)
    {
        return;
    }


    // Update surface reactions
    const scalar hReaction = cloud.surfaceReaction().calculate
    (
        dt,
        this->cell(),
        d,
        T,
        td.Tc(),
        td.pc(),
        td.rhoc(),
        mass,
        YGas,
        YLiquid,
        YSolid,
        YMix,
        N,
        dMassSRGas,
        dMassSRLiquid,
        dMassSRSolid,
        dMassSRCarrier
    );

    cloud.surfaceReaction().addToSurfaceReactionMass
    (
        this->nParticle_
       *(sum(dMassSRGas) + sum(dMassSRLiquid) + sum(dMassSRSolid))
    );

    const scalar xsi = min(T/cloud.constProps().TMax(), 1.0);
    const scalar coeff =
        (1.0 - xsi*xsi)*cloud.constProps().hRetentionCoeff();

    Sh += coeff*hReaction/dt;

    dhsTrans += (1.0 - coeff)*hReaction;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::ReactingMultiphaseParcel<ParcelType>::ReactingMultiphaseParcel
(
    const ReactingMultiphaseParcel<ParcelType>& p
)
:
    ParcelType(p),
    YGas_(p.YGas_),
    YLiquid_(p.YLiquid_),
    YSolid_(p.YSolid_),
    canCombust_(p.canCombust_)
{}


template<class ParcelType>
Foam::ReactingMultiphaseParcel<ParcelType>::ReactingMultiphaseParcel
(
    const ReactingMultiphaseParcel<ParcelType>& p,
    const polyMesh& mesh
)
:
    ParcelType(p, mesh),
    YGas_(p.YGas_),
    YLiquid_(p.YLiquid_),
    YSolid_(p.YSolid_),
    canCombust_(p.canCombust_)
{}


// * * * * * * * * * * * * * * IOStream operators  * * * * * * * * * * * * * //

#include "ReactingMultiphaseParcelIO.C"

// ************************************************************************* //
