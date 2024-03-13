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

#include "KinematicParcel.H"
#include "IOstreams.H"
#include "IOField.H"
#include "Cloud.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class ParcelType>
Foam::string Foam::KinematicParcel<ParcelType>::propertyList_ =
    Foam::KinematicParcel<ParcelType>::propertyList();

template<class ParcelType>
const std::size_t Foam::KinematicParcel<ParcelType>::sizeofFields_
(
    sizeof(KinematicParcel<ParcelType>)
  - offsetof(KinematicParcel<ParcelType>, active_)
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::KinematicParcel<ParcelType>::KinematicParcel
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields
)
:
    ParcelType(mesh, is, readFields),
    active_(false),
    typeId_(0),
    nParticle_(0.0),
    d_(0.0),
    dTarget_(0.0),
    dPacked_(0.0),
    U_(Zero),
    rho_(0.0),
    age_(0.0),
    tTurb_(0.0),
    UTurb_(Zero),
    Qdens_(0.0),	// E. Weiand - 27/02/2019 - electric particle charge density
    E_(Zero),		// E. Weiand - 30/10/2019 - electric field strength
    FLorentz_(Zero),	// E. Weiand - 30/10/2019 - Lorentz force
    particleRe_(0.0),	// E. Weiand - 26/11/2019 - particle Re
    localMaxCo_(GREAT)	// E. Weiand - 06/03/2020 - local max step fraction
{
    if (readFields)
    {
        if (is.format() == IOstream::ASCII)
        {
            active_ = readBool(is);
            typeId_ = readLabel(is);
            nParticle_ = readScalar(is);
            d_ = readScalar(is);
            dTarget_ = readScalar(is);
 	    dPacked_ = readScalar(is);
            is >> U_;
            rho_ = readScalar(is);
            age_ = readScalar(is);
            tTurb_ = readScalar(is);
	    Qdens_ = readScalar(is);		// E. Weiand - 27/02/2019 - electric particle charge density
            is >> UTurb_;
	    is >> E_;				// E. Weiand - 30/10/2019 - electric field strength
	    is >> FLorentz_;			// E. Weiand - 30/10/2019 - Lorentz force
	    particleRe_ = readScalar(is);	// E. Weiand - 30/10/2019 - particle Reynold number
	    localMaxCo_ = readScalar(is);	// E. Weiand - 06/03/2020 - local max step fraction
        }
        else
        {
            is.read(reinterpret_cast<char*>(&active_), sizeofFields_);
        }
    }

    // Check state of Istream
    is.check
    (
        "KinematicParcel<ParcelType>::KinematicParcel"
        "(const polyMesh&, Istream&, bool)"
    );
}


template<class ParcelType>
template<class CloudType>
void Foam::KinematicParcel<ParcelType>::readFields(CloudType& c)
{
    bool write = c.size();

    ParcelType::readFields(c);

    IOField<label> active
    (
        c.fieldIOobject("active", IOobject::MUST_READ),
        write
    );
    c.checkFieldIOobject(c, active);

    IOField<label> typeId
    (
        c.fieldIOobject("typeId", IOobject::MUST_READ),
        write
    );
    c.checkFieldIOobject(c, typeId);

    IOField<scalar> nParticle
    (
        c.fieldIOobject("nParticle", IOobject::MUST_READ),
        write
    );
    c.checkFieldIOobject(c, nParticle);

    IOField<scalar> d
    (
        c.fieldIOobject("d", IOobject::MUST_READ),
        write
    );
    c.checkFieldIOobject(c, d);

    IOField<scalar> dTarget
    (
        c.fieldIOobject("dTarget", IOobject::MUST_READ),
        write
    );
    c.checkFieldIOobject(c, dTarget);

    IOField<scalar> dPacked
    (
        c.fieldIOobject("dPacked", IOobject::MUST_READ),
        write
    );
    c.checkFieldIOobject(c, dPacked);

    IOField<vector> U
    (
        c.fieldIOobject("U", IOobject::MUST_READ),
        write
    );
    c.checkFieldIOobject(c, U);

    IOField<scalar> rho
    (
        c.fieldIOobject("rho", IOobject::MUST_READ),
        write
    );
    c.checkFieldIOobject(c, rho);

    IOField<scalar> age
    (
        c.fieldIOobject("age", IOobject::MUST_READ),
        write
    );
    c.checkFieldIOobject(c, age);

    IOField<scalar> tTurb
    (
        c.fieldIOobject("tTurb", IOobject::MUST_READ),
        write
    );
    c.checkFieldIOobject(c, tTurb);

    IOField<vector> UTurb
    (
        c.fieldIOobject("UTurb", IOobject::MUST_READ),
        write
    );
    c.checkFieldIOobject(c, UTurb);

    // E. Weiand - 27/02/2019 - electric particle charge density
    IOField<scalar> Qdens
    (
        c.fieldIOobject("Qdens", IOobject::MUST_READ),
        write
    );
    c.checkFieldIOobject(c, Qdens);

    // E. Weiand - 30/10/2019 - electric field strength
    IOField<vector> E
    (
        c.fieldIOobject("E", IOobject::MUST_READ),
        write
    );
    c.checkFieldIOobject(c, E);

    // E. Weiand - 30/10/2019 - Lorentz force
    IOField<vector> FLorentz
    (
        c.fieldIOobject("FLorentz", IOobject::MUST_READ),
        write
    );
    c.checkFieldIOobject(c, FLorentz);

    // E. Weiand - 26/11/2019 - particle Re
    IOField<scalar> particleRe
    (
        c.fieldIOobject("particleRe", IOobject::MUST_READ),
        write
    );
    c.checkFieldIOobject(c, particleRe);

    // E. Weiand - 06/03/2019 - local max step fraction
    IOField<scalar> localMaxCo
    (
        c.fieldIOobject("localMaxCo", IOobject::MUST_READ),
        write
    );
    c.checkFieldIOobject(c, localMaxCo);

    label i = 0;

    forAllIter(typename CloudType, c, iter)
    {
        KinematicParcel<ParcelType>& p = iter();

        p.active_ = active[i];
        p.typeId_ = typeId[i];
        p.nParticle_ = nParticle[i];
        p.d_ = d[i];
        p.dTarget_ = dTarget[i];
	p.dPacked_ = dPacked[i];
        p.U_ = U[i];
        p.rho_ = rho[i];
        p.age_ = age[i];
        p.tTurb_ = tTurb[i];
        p.UTurb_ = UTurb[i];
	p.Qdens_ = Qdens[i];		// E. Weiand - 27/02/2019 - electric particle charge density
	p.E_ = E[i];			// E. Weiand - 30/10/2019 - electric field strength
	p.FLorentz_ = FLorentz[i];	// E. Weiand - 30/10/2019 - Lorentz force
	p.particleRe_ = particleRe[i];	// E. Weiand - 26/11/2019 - particle Re
	p.localMaxCo_ = localMaxCo[i];	// E. Weiand - 06/03/2019 - local max step fraction

        i++;
    }
}


template<class ParcelType>
template<class CloudType>
void Foam::KinematicParcel<ParcelType>::writeFields(const CloudType& c)
{
    ParcelType::writeFields(c);

    label np = c.size();

    IOField<label> active(c.fieldIOobject("active", IOobject::NO_READ), np);
    IOField<label> typeId(c.fieldIOobject("typeId", IOobject::NO_READ), np);
    IOField<scalar> nParticle
    (
        c.fieldIOobject("nParticle", IOobject::NO_READ),
        np
    );
    IOField<scalar> d(c.fieldIOobject("d", IOobject::NO_READ), np);
    IOField<scalar> dTarget(c.fieldIOobject("dTarget", IOobject::NO_READ), np);
    IOField<scalar> dPacked(c.fieldIOobject("dPacked", IOobject::NO_READ), np);
    IOField<vector> U(c.fieldIOobject("U", IOobject::NO_READ), np);
    IOField<scalar> rho(c.fieldIOobject("rho", IOobject::NO_READ), np);
    IOField<scalar> age(c.fieldIOobject("age", IOobject::NO_READ), np);
    IOField<scalar> tTurb(c.fieldIOobject("tTurb", IOobject::NO_READ), np);
    IOField<vector> UTurb(c.fieldIOobject("UTurb", IOobject::NO_READ), np);
    IOField<scalar> Qdens(c.fieldIOobject("Qdens", IOobject::NO_READ), np);			// E. Weiand - 27/02/2019 - electric particle charge density
    IOField<vector> E(c.fieldIOobject("E", IOobject::NO_READ), np);				// E. Weiand - 18/09/2019 - electric field strength
    IOField<vector> FLorentz(c.fieldIOobject("FLorentz", IOobject::NO_READ), np);		// E. Weiand - 18/09/2019 - Lorentz force
    IOField<scalar> particleRe(c.fieldIOobject("particleRe", IOobject::NO_READ), np);		// E. Weiand - 26/11/2019 - particle Re
    IOField<scalar> localMaxCo(c.fieldIOobject("localMaxCo", IOobject::NO_READ), np);		// E. Weiand - 06/03/2019 - local max step fraction

    label i = 0;

    forAllConstIter(typename CloudType, c, iter)
    {
        const KinematicParcel<ParcelType>& p = iter();

        active[i] = p.active();
        typeId[i] = p.typeId();
        nParticle[i] = p.nParticle();
        d[i] = p.d();
        dTarget[i] = p.dTarget();
 	dPacked[i] = p.dPacked();
        U[i] = p.U();
        rho[i] = p.rho();
        age[i] = p.age();
        tTurb[i] = p.tTurb();
        UTurb[i] = p.UTurb();
	Qdens[i] = p.Qdens();		// E. Weiand - 27/02/2019 - electric particle charge density
	E[i] = p.E();			// E. Weiand - 18/09/2019 - electric field strength	
	FLorentz[i] = p.FLorentz();	// E. Weiand - 18/09/2019 - Lorentz force
	particleRe[i] = p.particleRe();	// E. Weiand - 26/11/2019 - particle Re
	localMaxCo[i] = p.localMaxCo();	// E. Weiand - 06/03/2019 - local max step fraction

        i++;
    }

    const bool write = np > 0;

    active.write(write);
    typeId.write(write);
    nParticle.write(write);
    d.write(write);
    dTarget.write(write);
    dPacked.write(write);
    U.write(write);
    rho.write(write);
    age.write(write);
    tTurb.write(write);
    UTurb.write(write);
    Qdens.write(write);		// E. Weiand - 18/09/2019 - electric particle charge
    E.write(write);		// E. Weiand - 30/10/2019 - electric field strenght
    FLorentz.write(write);	// E. Weiand - 30/10/2019 - Lorentz force
    particleRe.write(write);	// E. Weiand - 26/11/2019 - particle Re
    localMaxCo.write(write);	// E. Weiand - 06/03/2019 - local max step fraction
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class ParcelType>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const KinematicParcel<ParcelType>& p
)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const ParcelType&>(p)
            << token::SPACE << p.active()
            << token::SPACE << p.typeId()
            << token::SPACE << p.nParticle()
            << token::SPACE << p.d()
            << token::SPACE << p.dTarget()
	    << token::SPACE << p.dPacked()
            << token::SPACE << p.U()
            << token::SPACE << p.rho()
            << token::SPACE << p.age()
            << token::SPACE << p.tTurb()
            << token::SPACE << p.UTurb()
	    << token::SPACE << p.Qdens()	// E. Weiand - 27/02/2019 - electric particle charge density
	    << token::SPACE << p.E()		// E. Weiand - 30/10/2019 - electric field strenght
	    << token::SPACE << p.FLorentz()	// E. Weiand - 30/10/2019 - Lorentz force
	    << token::SPACE << p.particleRe()	// E. Weiand - 26/11/2019 - particle Re
	    << token::SPACE << p.localMaxCo();	// E. Weiand - 06/03/2019 - local max step fraction
    }
    else
    {
        os  << static_cast<const ParcelType&>(p);
        os.write
        (
            reinterpret_cast<const char*>(&p.active_),
            KinematicParcel<ParcelType>::sizeofFields_
        );
    }

    // Check state of Ostream
    os.check
    (
        "Ostream& operator<<(Ostream&, const KinematicParcel<ParcelType>&)"
    );

    return os;
}


// ************************************************************************* //
