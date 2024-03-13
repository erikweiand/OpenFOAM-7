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
#include "IOstreams.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class ParcelType>
Foam::string Foam::ReactingParcel<ParcelType>::propertyList_ =
    Foam::ReactingParcel<ParcelType>::propertyList();

template<class ParcelType>
const std::size_t Foam::ReactingParcel<ParcelType>::sizeofFields_
(
    sizeof(scalar)
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::ReactingParcel<ParcelType>::ReactingParcel
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields
)
:
    ParcelType(mesh, is, readFields),
    mass0_(0.0),
    Y_(0),
    d0_(0.0),
    Kevap_(0.0),
    Dp_(0.0),
    beta_(1.0),
    isDenselyPacked_(false),
    densMult_(1.0)
{
    if (readFields)
    {
        DynamicList<scalar> Ymix;

        if (is.format() == IOstream::ASCII)
        {
	    // E. Weiand - 2020/03/04: replaced, see below
//            is >> mass0_ >> Ymix;

	    mass0_ = readScalar(is);
	    is >> Ymix;
 	    d0_ = readScalar(is);
	    Kevap_ = readScalar(is);
	    Dp_ = readScalar(is);
	    beta_ = readScalar(is);
	    isDenselyPacked_ = readBool(is);
	    densMult_ = readScalar(is);
        }
        else
        {
            is.read(reinterpret_cast<char*>(&mass0_), sizeofFields_);
            is >> Ymix;
            is.read(reinterpret_cast<char*>(&d0_), sizeofFields_);
	    is.read(reinterpret_cast<char*>(&Kevap_), sizeofFields_);
	    is.read(reinterpret_cast<char*>(&Dp_), sizeofFields_);
	    is.read(reinterpret_cast<char*>(&beta_), sizeofFields_);
	    is.read(reinterpret_cast<char*>(&isDenselyPacked_), sizeofFields_);
	    is.read(reinterpret_cast<char*>(&densMult_), sizeofFields_);
        }

        Y_.transfer(Ymix);
    }

    // Check state of Istream
    is.check
    (
        "ReactingParcel<ParcelType>::ReactingParcel"
        "("
            "const polyMesh&, "
            "Istream&, "
            "bool"
        ")"
    );
}


template<class ParcelType>
template<class CloudType>
void Foam::ReactingParcel<ParcelType>::readFields(CloudType& c)
{
    ParcelType::readFields(c);
}


template<class ParcelType>
template<class CloudType, class CompositionType>
void Foam::ReactingParcel<ParcelType>::readFields
(
    CloudType& c,
    const CompositionType& compModel
)
{
    bool valid = c.size();

    ParcelType::readFields(c);

    IOField<scalar> mass0
    (
        c.fieldIOobject("mass0", IOobject::MUST_READ),
        valid
    );
    c.checkFieldIOobject(c, mass0);

    // E. Weiand - 2020/03/04
    IOField<scalar> d0
    (
        c.fieldIOobject("d0", IOobject::MUST_READ),
        valid
    );
    c.checkFieldIOobject(c, d0);

    // E. Weiand - 2020/03/04
    IOField<scalar> Kevap
    (
        c.fieldIOobject("Kevap", IOobject::MUST_READ),
        valid
    );
    c.checkFieldIOobject(c, Kevap);

    // E. Weiand - 2020/03/04
    IOField<scalar> Dp
    (
        c.fieldIOobject("Dp", IOobject::MUST_READ),
        valid
    );
    c.checkFieldIOobject(c, Dp);

    // E. Weiand - 2020/03/04
    IOField<scalar> beta
    (
        c.fieldIOobject("beta", IOobject::MUST_READ),
        valid
    );
    c.checkFieldIOobject(c, beta);

    // E. Weiand - 2020/09/29
    IOField<label> isDenselyPacked
    (
        c.fieldIOobject("isDenselyPacked", IOobject::MUST_READ),
        valid
    );
    c.checkFieldIOobject(c, isDenselyPacked);

    // E. Weiand - 2020/09/29
    IOField<scalar> densMult
    (
        c.fieldIOobject("densMult", IOobject::MUST_READ),
        valid
    );
    c.checkFieldIOobject(c, densMult);
    
    label i = 0;
    forAllIter(typename Cloud<ReactingParcel<ParcelType>>, c, iter)
    {
        ReactingParcel<ParcelType>& p = iter();
        p.mass0_ = mass0[i];
	p.d0_ = d0[i];			// E. Weiand - 2020/03/04
	p.Kevap_ = Kevap[i];		// E. Weiand - 2020/03/04
	p.Dp_ = Dp[i];			// E. Weiand - 2020/03/04
	p.beta_ = beta[i];		// E. Weiand - 2020/03/04
	p.isDenselyPacked_ = isDenselyPacked[i];		// E. Weiand - 2020/09/29
	p.densMult_ = densMult[i];		// E. Weiand - 2020/09/29
	
	i++;
    }

    // Get names and sizes for each Y...
    const wordList& phaseTypes = compModel.phaseTypes();
    const label nPhases = phaseTypes.size();
    wordList stateLabels(nPhases, "");
    if (compModel.nPhase() == 1)
    {
        stateLabels = compModel.stateLabels()[0];
    }


    // Set storage for each Y... for each parcel
    forAllIter(typename Cloud<ReactingParcel<ParcelType>>, c, iter)
    {
        ReactingParcel<ParcelType>& p = iter();
        p.Y_.setSize(nPhases, 0.0);
    }

    // Populate Y for each parcel
    forAll(phaseTypes, j)
    {
        IOField<scalar> Y
        (
            c.fieldIOobject
            (
                "Y" + phaseTypes[j] + stateLabels[j],
                 IOobject::MUST_READ
            ),
            valid
        );

        label i = 0;
        forAllIter(typename Cloud<ReactingParcel<ParcelType>>, c, iter)
        {
            ReactingParcel<ParcelType>& p = iter();
            p.Y_[j] = Y[i++];
        }
    }
}


template<class ParcelType>
template<class CloudType>
void Foam::ReactingParcel<ParcelType>::writeFields(const CloudType& c)
{
    ParcelType::writeFields(c);
}


template<class ParcelType>
template<class CloudType, class CompositionType>
void Foam::ReactingParcel<ParcelType>::writeFields
(
    const CloudType& c,
    const CompositionType& compModel
)
{
    ParcelType::writeFields(c);

    const label np = c.size();

    {
        IOField<scalar> mass0(c.fieldIOobject("mass0", IOobject::NO_READ), np);
        IOField<scalar> d0(c.fieldIOobject("d0", IOobject::NO_READ), np);		// E. Weiand -2020/03/04
        IOField<scalar> Kevap(c.fieldIOobject("Kevap", IOobject::NO_READ), np);		// E. Weiand -2020/03/04
        IOField<scalar> Dp(c.fieldIOobject("Dp", IOobject::NO_READ), np);		// E. Weiand -2020/03/04
        IOField<scalar> beta(c.fieldIOobject("beta", IOobject::NO_READ), np);		// E. Weiand -2020/03/04
        IOField<label> isDenselyPacked(c.fieldIOobject("isDenselyPacked", IOobject::NO_READ), np);		// E. Weiand -2020/09/29
        IOField<scalar> densMult(c.fieldIOobject("densMult", IOobject::NO_READ), np);		// E. Weiand -2020/09/29

        label i = 0;
        forAllConstIter(typename Cloud<ReactingParcel<ParcelType>>, c, iter)
        {
            const ReactingParcel<ParcelType>& p = iter();
            mass0[i] = p.mass0_;

	    d0[i] = p.d0_;		// E. Weiand - 2020/03/04	
	    Kevap[i] = p.Kevap_;	// E. Weiand - 2020/03/04
	    Dp[i] = p.Dp_;		// E. Weiand - 2020/03/04
	    beta[i] = p.beta_;		// E. Weiand - 2020/03/04
	    isDenselyPacked[i] = p.isDenselyPacked_;		// E. Weiand - 2020/09/29
	    densMult[i] = p.densMult_;		// E. Weiand - 2020/09/29

	    i++;
        }
        mass0.write(np > 0);
        d0.write(np > 0);	// E. Weiand - 2020/03/04
        Kevap.write(np > 0);	// E. Weiand - 2020/03/04
        Dp.write(np > 0);	// E. Weiand - 2020/03/04
        beta.write(np > 0);	// E. Weiand - 2020/03/04
        isDenselyPacked.write(np > 0);	// E. Weiand - 2020/09/29
        densMult.write(np > 0);	// E. Weiand - 2020/09/29

        // Write the composition fractions
        const wordList& phaseTypes = compModel.phaseTypes();
        wordList stateLabels(phaseTypes.size(), "");
        if (compModel.nPhase() == 1)
        {
            stateLabels = compModel.stateLabels()[0];
        }

        forAll(phaseTypes, j)
        {
            IOField<scalar> Y
            (
                c.fieldIOobject
                (
                    "Y" + phaseTypes[j] + stateLabels[j],
                    IOobject::NO_READ
                ),
                np
            );
            label i = 0;
            forAllConstIter
            (
                typename Cloud<ReactingParcel<ParcelType>>,
                c,
                iter
            )
            {
                const ReactingParcel<ParcelType>& p = iter();
                Y[i++] = p.Y()[j];
            }

            Y.write(np > 0);
        }
    }
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class ParcelType>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const ReactingParcel<ParcelType>& p
)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const ParcelType&>(p)
            << token::SPACE << p.mass0()
            << token::SPACE << p.Y()
	    << token::SPACE << p.d0()		// E. Weiand - 2020/03/04
	    << token::SPACE << p.Kevap()	// E. Weiand - 2020/03/04
	    << token::SPACE << p.Dp()		// E. Weiand - 2020/03/04
	    << token::SPACE << p.beta()		// E. Weiand - 2020/03/04
	    << token::SPACE << p.isDenselyPacked()	// E. Weiand - 2020/09/29
	    << token::SPACE << p.densMult();	// E. Weiand - 2020/09/29
    }
    else
    {
        os  << static_cast<const ParcelType&>(p);
        os.write
        (
            reinterpret_cast<const char*>(&p.mass0_),
            ReactingParcel<ParcelType>::sizeofFields_
        );
        os  << p.Y();
        os.write
        (
            reinterpret_cast<const char*>(&p.d0_),
            ReactingParcel<ParcelType>::sizeofFields_
        );
        os.write
        (
            reinterpret_cast<const char*>(&p.Kevap_),
            ReactingParcel<ParcelType>::sizeofFields_
        );
        os.write
        (
            reinterpret_cast<const char*>(&p.Dp_),
            ReactingParcel<ParcelType>::sizeofFields_
        );
        os.write
        (
            reinterpret_cast<const char*>(&p.beta_),
            ReactingParcel<ParcelType>::sizeofFields_
        );
        os.write
        (
            reinterpret_cast<const char*>(&p.isDenselyPacked_),
            ReactingParcel<ParcelType>::sizeofFields_
        );
        os.write
        (
            reinterpret_cast<const char*>(&p.densMult_),
            ReactingParcel<ParcelType>::sizeofFields_
        );
    }

    // Check state of Ostream
    os.check
    (
        "Ostream& operator<<(Ostream&, const ReactingParcel<ParcelType>&)"
    );

    return os;
}


// ************************************************************************* //
