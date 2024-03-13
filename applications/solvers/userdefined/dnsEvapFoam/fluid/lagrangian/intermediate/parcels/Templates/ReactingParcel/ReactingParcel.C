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

#include "ReactingParcel.H"
#include "specie.H"
#include "CompositionModel.H"
#include "PhaseChangeModel.H"
#include "mathematicalConstants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * * //

// M.P. Sitte, 15/03/2016: calcSurfaceValuesAS
template<class ParcelType>
template<class TrackCloudType>
void Foam::ReactingParcel<ParcelType>::calcSurfaceValuesAS
(
    TrackCloudType& cloud,
    trackingData& td,
    const label cellI,
    const scalar T,
    const label idPhase,  // M.P.Sitte
    const scalarField& YComponents, // M.P.Sitte
    scalar& Ts,
    scalarField& YYs,
    scalar& Ws,
    scalar& rhos,
    scalar& mus,
    scalar& Pr,
    scalar& kappas,
    scalar& BM
) const
{

    typedef typename TrackCloudType::reactingCloudType reactingCloudType;

    const SLGThermo& thermo = cloud.thermo();

    const CompositionModel<reactingCloudType>& composition = cloud.composition();


    // Molar mass of carrier
    //const scalar Wc = td.rhoc()*specie::RR*td.Tc()/td.pc();

    // Index of fuel in gas phase composition vector
    label idc = 100;
    label idl = 100;
    forAll(YComponents, i)
    {
        idc = composition.localToCarrierId(idPhase, i);

	// E. Weiand - 06/02/2020
//        idl = composition.globalIds(idPhase)[i]; 
	idl = i;
    }

    // Boiling temperature
    const scalar TBoil = composition.liquids().properties()[idl].pvInvert(td.pc());

    // limit droplet temperature
    scalar Td= max(T, cloud.constProps().TMin());
    Td= min(Td, 0.999*TBoil);

    // Reference condition: Ts and Ys
    // Film temperature using two thirds rule
    Ts = (2.0*Td + td.Tc())/3.0;

    // Far field carrier composition
    scalarField YYinf(thermo.carrier().species().size());
    forAll(YYinf, i)
    {
        YYinf[i]= thermo.carrier().Y(i)[cellI];
    }

    // Far field vapour mass fraction
    scalar Yc = YYinf[idc];

    // Molar mass of fuel - only correct for a one-component fuel
    scalar Wfuel = thermo.carrier().Wi(idc);

    // Carrier composition, excluding vapour
    scalarField YYcarr(YYinf);
    YYcarr[idc]= 0.0;
    YYcarr /= sum(YYcarr);

    // Carrier mean molecular mass, excluding vapour
    scalar Wcarr= 0.0;
    forAll(YYcarr, i)
    {
        Wcarr += YYcarr[i]/thermo.carrier().Wi(i);
    }
    Wcarr= 1.0/Wcarr; 

    // Surface vapour fraction
    //scalar Xsurf =  composition.liquids().properties()[idl].pv(td.pc(), Td) / td.pc();

    // E. Weiand - 2020/03/06
    // saturation pressure from Antoine equation for ethanol
    const scalar Apv = 5.37229;
    const scalar Bpv = 1670.409;
    const scalar Cpv = -40.191;
    scalar Xsurf = pow(10, Apv - Bpv/(Td+Cpv))/1.01325;		// normalized by pboil (in bar)

       // scalar Ysurf= Xsurf / ( Xsurf + (1.0-Xsurf)*(Wc/Wfuel) );  // wrong, carrier Wc should not contain fuel
    Xsurf= max(Xsurf, 0.0);

    scalar Ysurf= Xsurf*Wfuel / (Xsurf*Wfuel + (1.0-Xsurf)*Wcarr);

    scalarField YYsurf(thermo.carrier().species().size());
    forAll(YYsurf, i)
    {
        YYsurf[i]= YYinf[i] * (1.0-max(Ysurf,Yc)) / (1.0-Yc+SMALL);
    }
    YYsurf[idc]= 0.0;
    YYsurf[idc]= 1.0 - sum(YYsurf);


    // Film vapour fraction according to 2/3 rule
//    scalar Ys= (2.0*Ysurf + Yc)/3.0;	// E. Weiand - 06/02/2020 - unused variable

    Ws= 0.0;
    forAll(YYs, i)
    {
        YYs[i]= (2.0*YYsurf[i] +  YYinf[i])/3.0;
        Ws += YYs[i] / thermo.carrier().Wi(i);
    }

    Ws= 1.0 / Ws;

    // Spalding transfer number for mass
    BM= max((Ysurf - Yc), 0.0) / max(1.0 - Ysurf, SMALL);
    BM= max(BM, 1.0e-10);

    // Only a first guess of ref. cond. // M.P.Sitte
    // Assuming thermo props vary linearly with T for small d(T)
    const scalar TRatio = td.Tc()/Ts;

    rhos = td.rhoc()*TRatio;

    //Info << nl << "calcSurfaceValuesAS: mus = " << mus << nl << endl;

    tetIndices tetIs = this->currentTetIndices();
    mus = td.muInterp().interpolate(this->coordinates(), tetIs)/TRatio;
    kappas = td.kappaInterp().interpolate(this->coordinates(), tetIs)/TRatio;

    Pr = td.Cpc()*mus/kappas;
    Pr = max(ROOTVSMALL, Pr);
}




// M.P. Sitte, 05/04/2016: correctSurfValuesAS
template<class ParcelType>
template<class TrackCloudType>
void Foam::ReactingParcel<ParcelType>::correctSurfaceValuesAS
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar T,        // film temperature, ref. cond.
    const scalarField& Ys, // film composition, ref. cond.
    const scalar Ws,       // mean molecular mass, ref cond.
    scalar& rhos,
    scalar& mus,
    scalar& Prs,
    scalar& kappas,
    scalar& Cps
)
{

    const SLGThermo& thermo = cloud.thermo();


    // Far field carrier  molar fractions
    scalarField Xinf(thermo.carrier().species().size());

    // Surface carrier composition (molar fraction)
    scalarField Xs(Xinf.size());

    forAll(Xs, i)
    {
        Xs[i]= Ys[i] *  Ws/thermo.carrier().Wi(i);
    }

    rhos = 0;
    mus = 0;
    kappas = 0;
    Cps = 0;
    scalar sumYiSqrtW = 0;
    scalar sumYiCbrtW = 0;

    forAll(Ys, i)
    {
        const scalar W = thermo.carrier().Wi(i);
        const scalar sqrtW = sqrt(W);
        const scalar cbrtW = cbrt(W);

        rhos += Xs[i]*W;
        mus += Ys[i]*sqrtW*thermo.carrier().mu(i, td.pc(), T);
        kappas += Ys[i]*cbrtW*thermo.carrier().kappa(i, td.pc(), T);
//        Cps += Xs[i]*thermo.carrier().Cp(i, td.pc(), T);
        Cps += Ys[i]*thermo.carrier().Cp(i, td.pc(), T);	// E. Weiand - 2020/05/06: changed molar fraction to mixture fraction

        sumYiSqrtW += Ys[i]*sqrtW;
        sumYiCbrtW += Ys[i]*cbrtW;
    }

    Cps = max(Cps, ROOTVSMALL);

    rhos *= td.pc()/(RR*T);
    rhos = max(rhos, ROOTVSMALL);

    mus /= sumYiSqrtW;
    mus = max(mus, ROOTVSMALL);

    kappas /= sumYiCbrtW;
    kappas = max(kappas, ROOTVSMALL);

    Prs = Cps*mus/kappas;
}

// M.P. Sitte, 10/05/2016: <calcPhaseChangeAS>
// calcPhaseChange for Abramzon & Sirignano model
template<class ParcelType>
template<class TrackCloudType>
void Foam::ReactingParcel<ParcelType>::calcPhaseChangeAS
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt,
    const label cellI,
    const scalar Re,
    const scalar Pr,
    const scalar Ts,
    const scalar rhos,        // density at ref. cond. // M.P. Sitte
    const scalar nus,
    const scalar d,
    const scalar T,
    const scalar mass,
    const label idPhase,
    const scalar YPhase,
    const scalarField& YComponents,
    scalar& BM,		// E. Weiand - 2020/03/30
    scalarField& dMassPC,
    scalar& fmixS,           // M.P. Sitte
    scalar& grFmixS,    // M.P. Sitte
    scalar& DS,         // M.P. Sitte
    scalar& Sh,
    scalar& N,
    scalar& NCpW,
    scalarField& Cs
)
{
    typedef typename TrackCloudType::reactingCloudType reactingCloudType;
    PhaseChangeModel<reactingCloudType>& phaseChange = cloud.phaseChange();

    /*Pout << "T = " << T << endl;
    Pout << "Tvap = " << phaseChange.Tvap(YComponents) << endl;
    Pout << "YPhase = " << YPhase << endl;
    Pout << "active = " << phaseChange.active() << endl;*/

    if (!phaseChange.active())
    {
        return;
    }

    scalar Tvap = phaseChange.Tvap(YComponents);

    if (T < Tvap || YPhase < SMALL)
    {
        return;
    }

    const scalar TMax = phaseChange.TMax(td.pc(), YComponents);
    const scalar Tdash = min(T, TMax);
    //const scalar Tsdash = min(Ts, TMax);

    // E. Weiand - 2020/03/04: nanofuel evaporation properties
    const scalar d0 = this->d0();
    scalar& Kevap = this->Kevap();
    scalar& Dp = this->Dp();
    scalar& beta = this->beta(); 
    bool& isDenselyPacked = this->isDenselyPacked(); 
    scalar& tauD = this->localMaxCo();
    scalar& f0 = this->f0();
    scalar& dtPrev = this->dtPrev();
    bool& isFirstStep = this->isFirstStep();

    // Calculate mass transfer due to phase change
    // M.P. Sitte: <phaseChange.calculateAS>
    phaseChange.calculateAS
    (
        dt,
        cellI,
        Re,
        Pr,
        d,
        rhos,      // rho @ ref. cond.  // M.P.Sitte
        nus,
        T,    //Tdash, // droplet temp. // don't limit, M.P.Sitte
        Ts,   //Tsdash,                 // don't limit, M.P.Sitte
        td.pc(),
        td.Tc(),
        YComponents,
	d0,
	BM,
        dMassPC,
        fmixS,     // M.P. Sitte, 04/10/2017 
        grFmixS,   // M.P. Sitte, 04/10/2017 
        DS,        // M.P. Sitte, 04/10/2017 
	Kevap,	   // E. Weiand - 2020/03/04
	Dp,	   // E. Weiand - 2020/03/04
	beta,	   // E. Weiand - 2020/03/04
	tauD,	   // E. Weiand - 2020/03/06
	isDenselyPacked,	// E. Weiand - 2020/09/21
	f0,
	dtPrev,
	isFirstStep
    );

    // Limit phase change mass by availability of each specie
    dMassPC = min(mass*YPhase*YComponents, dMassPC);
 
    const scalar dMassTot = sum(dMassPC);

    // Add to cumulative phase change mass
    phaseChange.addToPhaseChangeMass(this->nParticle_*dMassTot);

    const CompositionModel<reactingCloudType>& composition =
        cloud.composition();

    forAll(dMassPC, i)
    {
        const label idc = composition.localToCarrierId(idPhase, i);

	// E. Weiand - 06/02/2020
//        const label idl = composition.globalIds(idPhase)[i];
	const label idl = i;

	// E. Weiand - 2020/03/10
	// NIST enthalpy of vaporization
	const scalar An = 50.43 * 1000.0 / 0.04607;
	const scalar alphan = -0.4475;
	const scalar betan = 0.4989;
	const scalar Tcn = 513.9;
	const scalar Tr = T/Tcn;

        const scalar dh = An*exp(-alphan*Tr)*pow(1-Tr, betan); //phaseChange.dh(idc, idl, td.pc(), Tdash);

        Sh -= dMassPC[i]*dh/dt;

    }
}



// M.P.Sitte: <calcHeatTransferAS>
template<class ParcelType>
template<class TrackCloudType>
Foam::scalar Foam::ReactingParcel<ParcelType>::calcHeatTransferAS
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt,
    const label cellI,
    const label idPhase,
    const scalar Re,
    const scalar Pr,
    const scalar Ts,
    const scalar rhos,
    const scalar nus,
    const scalar Cps,
    const scalar kappas,
    const scalar BM,
    const scalarField& dMass,
    const scalar Sh,
    scalar& dhsTrans,
    scalar& Sph
)
{
/* ---------------------------------------------------------------------------*/
// Temperature equation for a single droplet
// Abramzon & Sirginano model (using the notation from Miller et al. 1998)
//
// dTd/dt = f2*Nu/(3*PrG)*(CpG/CpL)*(1/tauD) * (TG - Td) +
//              + (LV/CpL)*(md_dot/md)
// 
//   where
// 
// f2= -md_dot/(md*BT) * (3*PrG*tauG/Nu)
//
// tauD= rhoG*d^2/(18*muG)
//
// Note:  md_dot <= 0
//        Sh     <= 0
//        Subscript G refers to "gas" and is evaluated at ref. cond. (Ts, Ys)
//
// Note 2:   WARNING,  Error in Miller et al. (1998). See Abramzon & Sirignano (1989).
//        Eqn should be,
//          dTd/dt = f2*Nu/(3*PrG)*(CpV/CpL)
//
/*----------------------------------------------------------------------------*/


    typedef typename TrackCloudType::reactingCloudType reactingCloudType;


    const scalar d = this->d();
    const scalar rho = this->rho();
    const scalar Td = this->T_;
    const scalar Tc= td.Tc();
    const scalar md = rho * pow(d,3.0)*pi/6.0;
    const scalar dmd = sum(dMass);
    const scalar md_dot= -dmd/dt;

    const CompositionModel<reactingCloudType>& composition =
          cloud.composition();

    // find label of fuel  - assumes one-component liquid
    label idl= 100;
    label idc= 100;
    forAll(dMass, i)
    {
	// E. Weiand - 06/02/2020
//          idl = composition.globalIds(idPhase)[i];
	  idl = i;

          idc = composition.localToCarrierId(idPhase, i);
    }

    // Liquid heat capacity
    // E. Weiand - 07/02/2020: faulty liquidProperties entries in default OF database.
    // Use constant property from CASE/constant/ dict for now.
    const scalar CpL= composition.liquids().properties()[idl].Cp(td.pc(), Td);
    //const scalar CpL = cloud.constProps().Cp0();

    // vapour heat capacity
    // const scalar CpV= composition.carrier().Cp(idc, td.pc(), Td); // use Td according to Abramzon & Sahzin, 
                                                                 // not right accoring to private communication with Sazhin

    // Note: composition.carrier ... does not use properties from /src/.../properties/liquidProperties/ but from /CASE/constant/...
    const scalar CpV= composition.carrier().Cp(idc, td.pc(), Ts);

    // Mean gas heat capacity @ ref. cond.
    //const scalar Cps= Cp;

    // E. Weiand - 2020/02/21: commented, see below
    // Diffusion coefficient of vapour @ ref. cond.    
//    const scalar Dab = composition.liquids().properties()[idl].D(td.pc(), Ts);

    // E. Weiand - 2020/02/21: directly read from constant properties
//    const scalar Dab = cloud.constProps().D0();
    const scalar Dfactor = cloud.constProps().D0();
    const scalar Dab = Dfactor*(8.967e-11*pow(Ts, 2.0) + 2.539e-08*Ts - 3.108e-06);	// MATLAB polynomial from Pinheiro

    // Schmidt number
    const scalar Sc = nus/Dab;

    //Pout << "Sc = " << Sc << endl;
    //Pout << "Dab = " << Dab << endl;

    // Mean Lewis number
    const scalar Le= Sc / Pr;
 
    // Sherwood number
    // Ranz Marshall: Sh= 2.0 + 0.6   * sqrt(Re) * cbrt(Sc)
    // Froessling:    Sh= 2.0 + 0.552 * sqrt(Re) * cbrt(Sc)
    // Clift:	      Sh= 1.0 + (1.0 + Re*Sc)^1/3 * f(Re) 	<---
    const scalar fRe = (Re < 1.0) ? 1.0 : pow(Re, 0.077);

    // E. Weiand - 2020/03/06: Re correction factor model as in Clift et al. 
    //const scalar Shs= 2.0 + 0.6 * sqrt(Re) * cbrt(Sc);
    const scalar Shs = 1.0 + pow(1.0 + Re*Sc, 1.0/3.0) * fRe;

    //Pout << "CpV = " << CpV << endl;
    //Pout << "Shs = " << Shs << endl;
    //Pout << "Re = " << Re << endl;
    //Pout << "Le = " << Le << endl;
    //Pout << "nus = " << nus << endl;

    // Nusselt number
    // E. Weiand - 2020/03/06: Re correction factor model as in Clift et al.
    //const scalar Nus= 2.0 + 0.6 * sqrt(Re) * cbrt(Pr);
    const scalar Nus = 1.0 + pow(1.0 + Re*Pr, 1.0/3.0) * fRe;

    //Pout << "Nus = " << Nus << endl;

    // limit BM: BM>0
    const scalar BMdash= max(BM, 1.0e-6);
    const scalar FM = pow(1.0 + BMdash,0.7) * log(1.0 + BMdash) / BMdash; 
    const scalar ShStar= 2.0 + (Shs-2.0)/FM;

    // Spalding transfer number for heat
    scalar BT;
    scalar phi= (CpV/Cps) * (Shs/Nus) / Le ;
    BT= pow((1.0+BMdash),phi) - 1.0;

    scalar NuStar= 2.0;  // Initial guess

    for(int i=1; i<=3; i++)
    {
          const scalar FT = pow(1.0 + BT, 0.7) * log(1.0 + BT)  / BT;
          NuStar= 2.0 + (Nus-2.0)/FT;
          phi= (CpV/Cps)*(ShStar/NuStar) / Le;
          BT= pow((1.0+BMdash),phi) - 1.0;
          if(BT==0.0)
          {
               BT= SMALL;
          }
    }

    //Pout << "phi = " << phi << endl;
    //Pout << "NuStar = " << NuStar << endl;
    //Pout << "BT = " << BT << endl;
    //Pout << "BM = " << BM << endl;

    // T-Equation
    //   dH/dt = d(md*h)/dt = md*dh/dt + h*dmd/dt = Q + m_dot*h_evap
    //                                                        h_evap = h_v(T_evap)
    //                        md*dh/dt = md*CpL*dTd/dt
    //   md*CpL*dTd/dt = Q + md_dot*(h_evap - h)
    //                                    h = h_d = h_l(Td)
    //
    //   Note: Sh= -dMassPC[i]*dh/dt, where dh is the latent heat, Lv= h_v - h_l


    const scalar bp = -md_dot*CpV/(md*BT*CpL);  // Abramzon Sirignano model

    /*Pout << "md_dot = " << md_dot << endl;
    Pout << "CpV = " << CpV << endl;
    Pout << "CpL = " << CpL << endl;
    Pout << "BT = " << BT << endl;
    Pout << "bp = " << bp << endl;*/

//    const scalar bp = (Nus*Cps*18.0*nus*rhos) / (3.0*Pr*CpL*rho*pow(d,2.0) + ROOTVSMALL);  // Rapid mixing model

    // Important !  Use CpV

    // E. Weiand - 19/02/2020: changed to avoid division by zero
    //const scalar ap = Tc + Sh/(CpL*md*bp);
    const scalar ap = Tc*bp + Sh/(CpL*md);

    // TODO: implement zero-evapration heat transfer (Fourier law) for shell blocked droplets
	/* Consider for modeling:
		- droplet heat-up to gas temperature?
		- micro-explosion modeling	*/

    /*Info << "Td = " << Td << endl;
    Info << "Tc = " << Tc << endl;
    Info << "ShHT = " << Sh << endl;
    Info << "md = " << md << endl;
    Info << "ap = " << ap << endl;
    Info << "dt = " << dt << endl;
    Info << "tauE(" << this->origId() << ") = " << md * CpL / (pi * d * NuStar * kappas) << endl;*/

/* -----------------------------------------------------------------------------------*/
// E. Weiand - 06/02/2020: Adapted to new OF.#
// Note: adapted to NEW integration scheme but making use of OLD variable definitions:
// bp = bp_New = bp_Old
// ap_New = bp * ap_Old
// Note 2: NO RADIATION!
/* -----------------------------------------------------------------------------------*/    
    /*IntegrationScheme<scalar>::integrationResult Tres =
        cloud.TIntegrator().integrate(Td, dt, ap*bp, bp);*/
    const scalar deltaT = cloud.TIntegrator().delta(Td, dt, ap, bp);
    scalar Tres = Td + deltaT;
   
/*    scalar f1Q = (ap - Td*bp);
    scalar dy2Q = 0.0;
    if (isFirstStepQ_ == true)
    {
	// Euler method
	dy2Q = f1Q*dt;
	isFirstStepQ_ = false;
    }
    else
    {
	// 2nd order Adams-Bashforth method
	dy2Q = 0.5*(dt + 2*dtPrevQ_)*f1Q - 0.5*dt*f0Q_;
    }
    f0Q_ = f1Q;
    dtPrevQ_ = dt;

    const scalar deltaT = dy2Q;
    scalar Tres = Td + deltaT;*/

/*    const scalar Qin = -md_dot * CpV * (Tc - Td) / BT;
    const scalar Qout = -Sh;
    const scalar dQ = Qin - Qout;
    const scalar deltaT = dQ * dt / (CpL * md);
    scalar Tres = Td + deltaT;*/

    scalar pGas = td.pc();
    //ag813 - Correction for negative pressure
    if (pGas < 0.0)
    {
       pGas = 101325.0;
    }

    scalar  Tboil = composition.liquids().properties()[idl].pvInvert(pGas);

    scalar Tnew = min(max(Tres, 240.0), 0.999*Tboil);
  
    // Heat sink/source to gas phase: dhsTrans (md_dot < 0 for evap)
    Sph = dt*(-md_dot)*CpV/BT;

    // ! Only true for analytical integration scheme !
    // E. Weiand - 07/02/2020
    // analytical integration average temperature
    const scalar expTerm = exp(min(50, -bp*dt));
    const scalar Tav = ap/bp + (Td - ap/bp) * (1-expTerm)/(bp*dt);

    // E. Weiand - 06/01/2020
//    dhsTrans += Sph*(Tres.average() - Tc);
    dhsTrans += Sph*(Tav - Tc);

    /*Pout << "Tres = " << Tres << endl;
    Pout << "Tav = " << Tav << endl;

    Pout << "Sph = " << Sph << endl;
    Pout << "dhsTrans = " << dhsTrans << endl;
    Pout << "deltaT = " << deltaT << endl;*/

      // Sph = dt*htc*As;  // previsously implemented in OF

    return Tnew;
}

template<class ParcelType>
template<class TrackCloudType>
void Foam::ReactingParcel<ParcelType>::calcPhaseChange
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt,
    const scalar Re,
    const scalar Pr,
    const scalar Ts,
    const scalar nus,
    const scalar d,
    const scalar T,
    const scalar mass,
    const label idPhase,
    const scalar YPhase,
    const scalarField& Y,
    scalarField& dMassPC,
    scalar& Sh,
    scalar& N,
    scalar& NCpW,
    scalarField& Cs
)
{
    typedef typename TrackCloudType::reactingCloudType reactingCloudType;
    const CompositionModel<reactingCloudType>& composition =
        cloud.composition();
    PhaseChangeModel<reactingCloudType>& phaseChange = cloud.phaseChange();

    if (!phaseChange.active() || (YPhase < small))
    {
        return;
    }

    scalarField X(composition.liquids().X(Y));

    scalar Tvap = phaseChange.Tvap(X);

    if (T < Tvap)
    {
        return;
    }

    const scalar TMax = phaseChange.TMax(td.pc(), X);

    const scalar Tdash = min(T, TMax);
    const scalar Tsdash = min(Ts, TMax);

    scalarField hmm(dMassPC);

    // Calculate mass transfer due to phase change
    phaseChange.calculate
    (
        dt,
        this->cell(),
        Re,
        Pr,
        d,
        nus,
        Tdash,
        Tsdash,
        td.pc(),
        td.Tc(),
        X,
        dMassPC
    );

    // Limit phase change mass by availability of each specie
    dMassPC = min(mass*YPhase*Y, dMassPC);

    const scalar dMassTot = sum(dMassPC);

    // Add to cumulative phase change mass
    phaseChange.addToPhaseChangeMass(this->nParticle_*dMassTot);

    forAll(dMassPC, i)
    {
        const label cid = composition.localToCarrierId(idPhase, i);

        const scalar dh = phaseChange.dh(cid, i, td.pc(), Tdash);
        Sh -= dMassPC[i]*dh/dt;
    }

    // Update molar emissions
    if (cloud.heatTransfer().BirdCorrection())
    {
        // Average molecular weight of carrier mix - assumes perfect gas
        const scalar Wc = td.rhoc()*RR*td.Tc()/td.pc();

        forAll(dMassPC, i)
        {
            const label cid = composition.localToCarrierId(idPhase, i);

            const scalar Cp = composition.carrier().Cp(cid, td.pc(), Tsdash);
            const scalar W = composition.carrier().Wi(cid);
            const scalar Ni = dMassPC[i]/(this->areaS(d)*dt*W);

            const scalar Dab =
                composition.liquids().properties()[i].D(td.pc(), Tsdash, Wc);

            // Molar flux of species coming from the particle (kmol/m^2/s)
            N += Ni;

            // Sum of Ni*Cpi*Wi of emission species
            NCpW += Ni*Cp*W;

            // Concentrations of emission species
            Cs[cid] += Ni*d/(2.0*Dab);
        }
    }
}


template<class ParcelType>
Foam::scalar Foam::ReactingParcel<ParcelType>::updateMassFraction
(
    const scalar mass0,
    const scalarField& dMass,
    scalarField& Y
) const
{
    scalar mass1 = mass0 - sum(dMass);

    // only update the mass fractions if the new particle mass is finite
    if (mass1 > rootVSmall)
    {
        forAll(Y, i)
        {
            Y[i] = (Y[i]*mass0 - dMass[i])/mass1;
        }
    }

    return mass1;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::ReactingParcel<ParcelType>::ReactingParcel
(
    const ReactingParcel<ParcelType>& p
)
:
    ParcelType(p),
    mass0_(p.mass0_),
    Y_(p.Y_),
    d0_(p.d0_),		// E. Weiand - 2020/03/04
    Kevap_(p.Kevap_),	// E. Weiand - 2020/03/04
    Dp_(p.Dp_),		// E. Weiand - 2020/03/04
    beta_(p.beta_),	// E. Weiand - 2020/03/04
    isDenselyPacked_(p.isDenselyPacked_),	// E. Weiand - 2020/09/29
    densMult_(p.densMult_),	// E. Weiand - 2020/09/29
    f0_(0.0),
    dtPrev_(0.0),
    dtPrevQ_(0.0),
    isFirstStep_(true),
    isFirstStepQ_(true)
{
    d0_ = this->d_;	// E. Weiand - 2020/03/04
}


template<class ParcelType>
Foam::ReactingParcel<ParcelType>::ReactingParcel
(
    const ReactingParcel<ParcelType>& p,
    const polyMesh& mesh
)
:
    ParcelType(p, mesh),
    mass0_(p.mass0_),
    Y_(p.Y_),
    d0_(p.d0_),		// E. Weiand - 2020/03/04
    Kevap_(p.Kevap_),	// E. Weiand - 2020/03/04
    Dp_(p.Dp_),		// E. Weiand - 2020/03/04
    beta_(p.beta_),	// E. Weiand - 2020/03/04
    isDenselyPacked_(p.isDenselyPacked_),	// E. Weiand - 2020/09/29
    densMult_(p.densMult_),	// E. Weiand - 2020/09/29
    f0_(0.0),
    f0Q_(0.0),
    dtPrev_(0.0),
    dtPrevQ_(0.0),
    isFirstStep_(true),
    isFirstStepQ_(true)
{
    d0_ = this->d_;	// E. Weiand - 2020/03/04
}


// * * * * * * * * * * * * *  Member Functions * * * * * * * * * * * * * * * //

template<class ParcelType>
template<class TrackCloudType>
void Foam::ReactingParcel<ParcelType>::setCellValues
(
    TrackCloudType& cloud,
    trackingData& td
)
{
    ParcelType::setCellValues(cloud, td);

    td.pc() = td.pInterp().interpolate
    (
        this->coordinates(),
        this->currentTetIndices()
    );

    if (td.pc() < cloud.constProps().pMin())
    {
        if (debug)
        {
            WarningInFunction
                << "Limiting observed pressure in cell " << this->cell()
                << " to " << cloud.constProps().pMin() <<  nl << endl;
        }

        td.pc() = cloud.constProps().pMin();
    }
}


template<class ParcelType>
template<class TrackCloudType>
void Foam::ReactingParcel<ParcelType>::cellValueSourceCorrection
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt
)
{
    scalar addedMass = 0.0;
    scalar maxMassI = 0.0;
    forAll(cloud.rhoTrans(), i)
    {
        scalar dm = cloud.rhoTrans(i)[this->cell()];
        maxMassI = max(maxMassI, mag(dm));
        addedMass += dm;
    }

    if (maxMassI < rootVSmall)
    {
        return;
    }

    const scalar massCell = this->massCell(td);

    td.rhoc() += addedMass/cloud.pMesh().cellVolumes()[this->cell()];

    const scalar massCellNew = massCell + addedMass;
    td.Uc() = (td.Uc()*massCell + cloud.UTrans()[this->cell()])/massCellNew;

    scalar CpEff = 0.0;
    forAll(cloud.rhoTrans(), i)
    {
        scalar Y = cloud.rhoTrans(i)[this->cell()]/addedMass;
        CpEff += Y*cloud.composition().carrier().Cp(i, td.pc(), td.Tc());
    }

    const scalar Cpc = td.CpInterp().psi()[this->cell()];
    td.Cpc() = (massCell*Cpc + addedMass*CpEff)/massCellNew;

    td.Tc() += cloud.hsTrans()[this->cell()]/(td.Cpc()*massCellNew);

    if (td.Tc() < cloud.constProps().TMin())
    {
        if (debug)
        {
            WarningInFunction
                << "Limiting observed temperature in cell " << this->cell()
                << " to " << cloud.constProps().TMin() <<  nl << endl;
        }

        td.Tc() = cloud.constProps().TMin();
    }
}


template<class ParcelType>
template<class TrackCloudType>
void Foam::ReactingParcel<ParcelType>::correctSurfaceValues
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar T,
    const scalarField& Cs,
    scalar& rhos,
    scalar& mus,
    scalar& Prs,
    scalar& kappas
)
{
    // No correction if total concentration of emitted species is small
    if (!cloud.heatTransfer().BirdCorrection() || (sum(Cs) < small))
    {
        return;
    }

    const SLGThermo& thermo = cloud.thermo();

    // Far field carrier  molar fractions
    scalarField Xinf(thermo.carrier().species().size());

    forAll(Xinf, i)
    {
        Xinf[i] = thermo.carrier().Y(i)[this->cell()]/thermo.carrier().Wi(i);
    }
    Xinf /= sum(Xinf);

    // Molar fraction of far field species at particle surface
    const scalar Xsff = 1.0 - min(sum(Cs)*RR*this->T_/td.pc(), 1.0);

    // Surface carrier total molar concentration
    const scalar CsTot = td.pc()/(RR*this->T_);

    // Surface carrier composition (molar fraction)
    scalarField Xs(Xinf.size());

    // Surface carrier composition (mass fraction)
    scalarField Ys(Xinf.size());

    forAll(Xs, i)
    {
        // Molar concentration of species at particle surface
        const scalar Csi = Cs[i] + Xsff*Xinf[i]*CsTot;

        Xs[i] = (2.0*Csi + Xinf[i]*CsTot)/3.0;
        Ys[i] = Xs[i]*thermo.carrier().Wi(i);
    }
    Xs /= sum(Xs);
    Ys /= sum(Ys);


    rhos = 0;
    mus = 0;
    kappas = 0;
    scalar Cps = 0;
    scalar sumYiSqrtW = 0;
    scalar sumYiCbrtW = 0;

    forAll(Ys, i)
    {
        const scalar W = thermo.carrier().Wi(i);
        const scalar sqrtW = sqrt(W);
        const scalar cbrtW = cbrt(W);

        rhos += Xs[i]*W;
        mus += Ys[i]*sqrtW*thermo.carrier().mu(i, td.pc(), T);
        kappas += Ys[i]*cbrtW*thermo.carrier().kappa(i, td.pc(), T);
        Cps += Xs[i]*thermo.carrier().Cp(i, td.pc(), T);

        sumYiSqrtW += Ys[i]*sqrtW;
        sumYiCbrtW += Ys[i]*cbrtW;
    }

    Cps = max(Cps, rootVSmall);

    rhos *= td.pc()/(RR*T);
    rhos = max(rhos, rootVSmall);

    mus /= sumYiSqrtW;
    mus = max(mus, rootVSmall);

    kappas /= sumYiCbrtW;
    kappas = max(kappas, rootVSmall);

    Prs = Cps*mus/kappas;

    //Info << "Prs = " << Prs << endl;
    //Info << "kappas = " << kappas << endl;
}


template<class ParcelType>
template<class TrackCloudType>
void Foam::ReactingParcel<ParcelType>::calc
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt
)
{
    typedef typename TrackCloudType::reactingCloudType reactingCloudType;
    const CompositionModel<reactingCloudType>& composition =
        cloud.composition();


    // Define local properties at beginning of time step
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    const scalar np0 = this->nParticle_;
    const scalar d0 = this->d_;
    const vector& U0 = this->U_;
    const scalar T0 = this->T_;
    const scalar mass0 = this->mass();


    // Calc surface values
    scalar Ts, rhos, mus, Prs, kappas;
    this->calcSurfaceValues(cloud, td, T0, Ts, rhos, mus, Prs, kappas);
    scalar Res = this->Re(rhos, U0, td.Uc(), d0, mus);


    // Sources
    // ~~~~~~~

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

    // Phase change
    // ~~~~~~~~~~~~

    // Mass transfer due to phase change
    scalarField dMassPC(Y_.size(), 0.0);

    // Molar flux of species emitted from the particle (kmol/m^2/s)
    scalar Ne = 0.0;

    // Sum Ni*Cpi*Wi of emission species
    scalar NCpW = 0.0;

    // Surface concentrations of emitted species
    scalarField Cs(composition.carrier().species().size(), 0.0);

    // Calc mass and enthalpy transfer due to phase change
    calcPhaseChange
    (
        cloud,
        td,
        dt,
        Res,
        Prs,
        Ts,
        mus/rhos,
        d0,
        T0,
        mass0,
        0,
        1.0,
        Y_,
        dMassPC,
        Sh,
        Ne,
        NCpW,
        Cs
    );


    // 2. Update the parcel properties due to change in mass
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    scalarField dMass(dMassPC);
    scalar mass1 = updateMassFraction(mass0, dMass, Y_);

    this->Cp_ = composition.Cp(0, Y_, td.pc(), T0);

    // Update particle density or diameter
    if (cloud.constProps().constantVolume())
    {
        this->rho_ = mass1/this->volume();
    }
    else
    {
        this->d_ = cbrt(mass1/this->rho_*6.0/pi);
    }

    // Remove the particle when mass falls below minimum threshold
    if (np0*mass1 < cloud.constProps().minParcelMass())
    {
        td.keepParticle = false;

        if (cloud.solution().coupled())
        {
            scalar dm = np0*mass0;

            // Absorb parcel into carrier phase
            forAll(Y_, i)
            {
                scalar dmi = dm*Y_[i];
                label gid = composition.localToCarrierId(0, i);
                scalar hs = composition.carrier().Hs(gid, td.pc(), T0);

                cloud.rhoTrans(gid)[this->cell()] += dmi;
                cloud.hsTrans()[this->cell()] += dmi*hs;
            }
            cloud.UTrans()[this->cell()] += dm*U0;

            cloud.phaseChange().addToPhaseChangeMass(np0*mass1);
        }

        return;
    }

    // Correct surface values due to emitted species
    correctSurfaceValues(cloud, td, Ts, Cs, rhos, mus, Prs, kappas);
    Res = this->Re(rhos, U0, td.Uc(), this->d(), mus);


    // 3. Compute heat- and momentum transfers
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Heat transfer
    // ~~~~~~~~~~~~~

    // Calculate new particle temperature
    this->T_ =
        this->calcHeatTransfer
        (
            cloud,
            td,
            dt,
            Res,
            Prs,
            kappas,
            NCpW,
            Sh,
            dhsTrans,
            Sph
        );

    this->Cp_ = composition.Cp(0, Y_, td.pc(), T0);


    // Motion
    // ~~~~~~

    // Erik Weiand - 26/11/2019
    // Calculate new particle Reynolds number (output var)
    const scalar Re = this->Re(td);
    this->particleRe_ = Re;

    // Erik Weiand - 26/11/2019
    // Update electric field & Lorentz force
    this->E_ = this->calcElField(cloud, td);
    const scalar Q = this->Qdens_ * pi/6.0 * pow(this->d0_, 3.0);	// charge contained in nanoparticles
    this->FLorentz_ = this->E_ * Q;

    // Calculate new particle velocity
    this->U_ =
        this->calcVelocity(cloud, td, dt, Res, mus, mass1, Su, dUTrans, Spu);


    // 4. Accumulate carrier phase source terms
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (cloud.solution().coupled())
    {
        // Transfer mass lost to carrier mass, momentum and enthalpy sources
        forAll(dMass, i)
        {
            scalar dm = np0*dMass[i];
            label gid = composition.localToCarrierId(0, i);
            scalar hs = composition.carrier().Hs(gid, td.pc(), T0);

            cloud.rhoTrans(gid)[this->cell()] += dm;
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


// * * * * * * * * * * * * * * IOStream operators  * * * * * * * * * * * * * //

#include "ReactingParcelIO.C"

// ************************************************************************* //
