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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
inline const Foam::ReactingCloud<CloudType>&
Foam::ReactingCloud<CloudType>::cloudCopy() const
{
    return cloudCopyPtr_();
}


template<class CloudType>
inline const typename CloudType::particleType::constantProperties&
Foam::ReactingCloud<CloudType>::constProps() const
{
    return constProps_;
}


template<class CloudType>
inline typename CloudType::particleType::constantProperties&
Foam::ReactingCloud<CloudType>::constProps()
{
    return constProps_;
}


template<class CloudType>
inline const Foam::CompositionModel<Foam::ReactingCloud<CloudType>>&
Foam::ReactingCloud<CloudType>::composition() const
{
    return compositionModel_;
}


template<class CloudType>
inline const Foam::PhaseChangeModel<Foam::ReactingCloud<CloudType>>&
Foam::ReactingCloud<CloudType>::phaseChange() const
{
    return phaseChangeModel_;
}


template<class CloudType>
inline Foam::PhaseChangeModel<Foam::ReactingCloud<CloudType>>&
Foam::ReactingCloud<CloudType>::phaseChange()
{
    return phaseChangeModel_();
}


template<class CloudType>
inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::ReactingCloud<CloudType>::rhoTrans(const label i)
{
    return rhoTrans_[i];
}


template<class CloudType>
inline
const Foam::PtrList<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>&
Foam::ReactingCloud<CloudType>::rhoTrans() const
{
    return rhoTrans_;
}


template<class CloudType>
inline Foam::PtrList<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>&
Foam::ReactingCloud<CloudType>::rhoTrans()
{
    return rhoTrans_;
}


// M.P. Sitte, 04/10/2017    
template<class CloudType>
inline
const Foam::autoPtr<Foam::DimensionedField<Foam::scalar, Foam::volMesh> >&
Foam::ReactingCloud<CloudType>::FsNp() const
{
    return FsNp_;
}

// M.P. Sitte, 04/10/2017  
template<class CloudType>
inline Foam::autoPtr<Foam::DimensionedField<Foam::scalar, Foam::volMesh> >&
Foam::ReactingCloud<CloudType>::FsNp()
{
    return FsNp_;
}

// M.P. Sitte, 04/10/2017    
template<class CloudType>
inline
const Foam::autoPtr<Foam::DimensionedField<Foam::scalar, Foam::volMesh> >&
Foam::ReactingCloud<CloudType>::GrFsNp() const
{
    return GrFsNp_;
}

// M.P. Sitte, 04/10/2017  
template<class CloudType>
inline Foam::autoPtr<Foam::DimensionedField<Foam::scalar, Foam::volMesh> >&
Foam::ReactingCloud<CloudType>::GrFsNp()
{
    return GrFsNp_;
}

// M.P. Sitte, 04/10/2017    
template<class CloudType>
inline
const Foam::autoPtr<Foam::DimensionedField<Foam::scalar, Foam::volMesh> >&
Foam::ReactingCloud<CloudType>::DsNp() const
{
    return DsNp_;
}

// M.P. Sitte, 04/10/2017  
template<class CloudType>
inline Foam::autoPtr<Foam::DimensionedField<Foam::scalar, Foam::volMesh> >&
Foam::ReactingCloud<CloudType>::DsNp()
{
    return DsNp_;
}



// M.P. Sitte, 28/11/2017    
template<class CloudType>
inline
const Foam::autoPtr<Foam::DimensionedField<Foam::scalar, Foam::volMesh> >&
Foam::ReactingCloud<CloudType>::rdNp() const
{
    return rdNp_;
}

// M.P. Sitte, 28/11/2017  
template<class CloudType>
inline Foam::autoPtr<Foam::DimensionedField<Foam::scalar, Foam::volMesh> >&
Foam::ReactingCloud<CloudType>::rdNp()
{
    return rdNp_;
}

// M.P. Sitte, 28/11/2017    
template<class CloudType>
inline
const Foam::autoPtr<Foam::DimensionedField<Foam::scalar, Foam::volMesh> >&
Foam::ReactingCloud<CloudType>::ndNp() const
{
    return ndNp_;
}

// M.P. Sitte, 28/11/2017  
template<class CloudType>
inline Foam::autoPtr<Foam::DimensionedField<Foam::scalar, Foam::volMesh> >&
Foam::ReactingCloud<CloudType>::ndNp()
{
    return ndNp_;
}




// M.P. Sitte, 04/10/2017    
template<class CloudType>
inline
const Foam::autoPtr<Foam::DimensionedField<Foam::scalar, Foam::volMesh> >&
Foam::ReactingCloud<CloudType>::Np() const
{
    return Np_;
}

// M.P. Sitte, 04/10/2017
template<class CloudType>
inline Foam::autoPtr<Foam::DimensionedField<Foam::scalar, Foam::volMesh> >&
Foam::ReactingCloud<CloudType>::Np()
{
    return Np_;
}


// M.P. Sitte, 04/10/2017
template<class CloudType>
inline
const Foam::autoPtr<Foam::DimensionedField<Foam::scalar, Foam::volMesh> >&
Foam::ReactingCloud<CloudType>::VolNpDt() const
{
    return VolNpDt_;
}

// M.P. Sitte, 04/10/2017
template<class CloudType>
inline Foam::autoPtr<Foam::DimensionedField<Foam::scalar, Foam::volMesh> >&
Foam::ReactingCloud<CloudType>::VolNpDt()
{
    return VolNpDt_;
}



// M.P. Sitte, 04/10/2017
template<class CloudType>
inline
const Foam::autoPtr<Foam::DimensionedField<Foam::scalar, Foam::volMesh> >&
Foam::ReactingCloud<CloudType>::NpDt() const
{
    return NpDt_;
}

// M.P. Sitte, 04/10/2017
template<class CloudType>
inline Foam::autoPtr<Foam::DimensionedField<Foam::scalar, Foam::volMesh> >&
Foam::ReactingCloud<CloudType>::NpDt()
{
    return NpDt_;
}

// M.P. Sitte, 04/10/2017
template<class CloudType>
inline
const Foam::autoPtr<Foam::DimensionedField<Foam::scalar, Foam::volMesh> >&
Foam::ReactingCloud<CloudType>::massd() const
{
    return massd_;
}

// M.P. Sitte, 04/10/2017
template<class CloudType>
inline Foam::autoPtr<Foam::DimensionedField<Foam::scalar, Foam::volMesh> >&
Foam::ReactingCloud<CloudType>::massd()
{
    return massd_;
}



template<class CloudType>
inline Foam::tmp<Foam::fvScalarMatrix> Foam::ReactingCloud<CloudType>::SYi
(
    const label i,
    volScalarField& Yi
) const
{
    if (this->solution().coupled())
    {
        if (this->solution().semiImplicit("Yi"))
        {
            tmp<volScalarField> trhoTrans
            (
                volScalarField::New
                (
                    this->name() + ":rhoTrans",
                    this->mesh(),
                    dimensionedScalar(dimMass/dimTime/dimVolume, 0)
                )
            );

            volScalarField& sourceField = trhoTrans.ref();

            sourceField.primitiveFieldRef() =
                rhoTrans_[i]/(this->db().time().deltaTValue()*this->mesh().V());

            const dimensionedScalar Yismall("Yismall", dimless, small);

            return
                fvm::Sp(neg(sourceField)*sourceField/(Yi + Yismall), Yi)
              + pos0(sourceField)*sourceField;
        }
        else
        {
            tmp<fvScalarMatrix> tfvm(new fvScalarMatrix(Yi, dimMass/dimTime));
            fvScalarMatrix& fvm = tfvm.ref();

            fvm.source() = -rhoTrans_[i]/this->db().time().deltaTValue();

            return tfvm;
        }
    }

    return tmp<fvScalarMatrix>(new fvScalarMatrix(Yi, dimMass/dimTime));
}


template<class CloudType>
inline Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::ReactingCloud<CloudType>::Srho(const label i) const
{
    tmp<volScalarField::Internal> tRhoi
    (
        volScalarField::Internal::New
        (
            this->name() + ":rhoTrans",
            this->mesh(),
            dimensionedScalar
            (
                rhoTrans_[0].dimensions()/dimTime/dimVolume,
                0
            )
        )
    );

    if (this->solution().coupled())
    {
        scalarField& rhoi = tRhoi.ref();
        rhoi = rhoTrans_[i]/(this->db().time().deltaTValue()*this->mesh().V());
    }

    return tRhoi;
}


template<class CloudType>
inline Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::ReactingCloud<CloudType>::Srho() const
{
    tmp<volScalarField::Internal> trhoTrans
    (
        volScalarField::Internal::New
        (
            this->name() + ":rhoTrans",
            this->mesh(),
            dimensionedScalar
            (
                rhoTrans_[0].dimensions()/dimTime/dimVolume,
                0
            )
        )
    );

    Info<< "  !! ReactingCloudI.H ... Srho() !!" << endl;  // M.P. Sitte, 10/10/2015

    if (this->solution().coupled())
    {
        scalarField& sourceField = trhoTrans.ref();
        forAll(rhoTrans_, i)
        {
            sourceField += rhoTrans_[i];
        }

        sourceField /= this->db().time().deltaTValue()*this->mesh().V();
    }

    return trhoTrans;
}


template<class CloudType>
inline Foam::tmp<Foam::fvScalarMatrix>
Foam::ReactingCloud<CloudType>::Srho(volScalarField& rho) const
{
    if (this->solution().coupled())
    {
        tmp<volScalarField> trhoTrans
        (
            volScalarField::New
            (
                this->name() + ":rhoTrans",
                this->mesh(),
                dimensionedScalar(dimMass/dimTime/dimVolume, 0)
            )
        );

        scalarField& sourceField = trhoTrans.ref();

        if (this->solution().semiImplicit("rho"))
        {

            forAll(rhoTrans_, i)
            {
                sourceField += rhoTrans_[i];
            }
            sourceField /= this->db().time().deltaTValue()*this->mesh().V();

            return fvm::SuSp(trhoTrans()/rho, rho);
        }
        else
        {
            tmp<fvScalarMatrix> tfvm(new fvScalarMatrix(rho, dimMass/dimTime));
            fvScalarMatrix& fvm = tfvm.ref();

            forAll(rhoTrans_, i)
            {
                sourceField += rhoTrans_[i];
            }

            fvm.source() = -trhoTrans()/this->db().time().deltaT();

            return tfvm;
        }
    }

    return tmp<fvScalarMatrix>(new fvScalarMatrix(rho, dimMass/dimTime));
}


// M.P. Sitte, 04/10/2017
template<class CloudType>
inline Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh> >
Foam::ReactingCloud<CloudType>::FmixSurf() const
{

    tmp<DimensionedField<scalar, volMesh> > tFmixS
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                this->name() + ":tFmixS",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0), 0.0  )
        )
    );

    // Average mixture fraction on the droplet surface
    scalarField& sourceField = tFmixS(); 
    sourceField += FsNp_;
    sourceField /= ( Np_ + SMALL);

    Info<< "  !! ReactingCloudI.H ... FmixSurf() !!" << endl;  // M.P. Sitte, 12/10/2015

    return tFmixS;
}



// M.P. Sitte, 04/10/2017
template<class CloudType>
inline Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh> >
Foam::ReactingCloud<CloudType>::DSurf() const
{

    tmp<DimensionedField<scalar, volMesh> > tDS
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                this->name() + ":tDS",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0), 0.0  )
        )
    );

    // Average mixture fraction on the droplet surface
    scalarField& sourceField = tDS();
    sourceField += DsNp_;
    sourceField /= ( Np_ + SMALL);

    Info<< "  !! ReactingCloudI.H ... DSurf() !!" << endl;  // M.P. Sitte, 12/10/2015

    return tDS;
}



// M.P. Sitte, 04/10/2017
template<class CloudType>
inline Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh> >
Foam::ReactingCloud<CloudType>::GrFmixSurf() const
{

    tmp<DimensionedField<scalar, volMesh> > tGrFmixS
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                this->name() + ":tGrFmixS",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0), 0.0  )
        )
    );

    // Average mixture fraction on the droplet surface
    scalarField& sourceField = tGrFmixS();
    sourceField += GrFsNp_;
    sourceField /= ( Np_ + SMALL);

    Info<< "  !! ReactingCloudI.H ... GrFmixSurf() !!" << endl;  // M.P. Sitte, 12/10/2015

    return tGrFmixS;
}


// M.P. Sitte, 28/11/2017
template<class CloudType>
inline Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh> >
Foam::ReactingCloud<CloudType>::rDrop() const
{

    tmp<DimensionedField<scalar, volMesh> > trDrop
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                this->name() + ":trDrop",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0), 0.0  )
        )
    );

    // Average mixture fraction on the droplet surface
    scalarField& sourceField = trDrop();
    sourceField += rdNp_;
    sourceField /= ( Np_ + SMALL);

    Info<< "  !! ReactingCloudI.H ... rDrop() !!" << endl;  // M.P. Sitte, 12/10/2015

    return trDrop;
}



// M.P. Sitte, 28/11/2017
template<class CloudType>
inline Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh> >
Foam::ReactingCloud<CloudType>::nDrop() const
{

    tmp<DimensionedField<scalar, volMesh> > tnDrop
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                this->name() + ":tnDrop",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0), 0.0  )
        )
    );

    // Average mixture fraction on the droplet surface
    scalarField& sourceField = tnDrop();
    sourceField += ndNp_;
    sourceField /= ( Np_ + SMALL);

    Info<< "  !! ReactingCloudI.H ... nDrop() !!" << endl;  // M.P. Sitte, 12/10/2015

    return tnDrop;
}


// M.P. Sitte, 04/10/2017
template<class CloudType>
inline Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh> >
Foam::ReactingCloud<CloudType>::VolDrop() const
{

    tmp<DimensionedField<scalar, volMesh> > tVolDrop
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                this->name() + ":tVolDrop",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimensionedScalar("zero", dimensionSet(0, 3, 0, 0, 0), 0.0  )
        )
    );

    // Average mixture fraction on the droplet surface
    scalarField& sourceField = tVolDrop();
    sourceField += VolNpDt_;
    sourceField /= ( NpDt_ + pow(SMALL,3.0));

    Info<< "  !! ReactingCloudI.H ... VolDrop() !!" << endl;  // M.P. Sitte, 12/10/2015

    return tVolDrop;
}



// M.P. Sitte, 04/10/2017
template<class CloudType>
inline Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh> >
Foam::ReactingCloud<CloudType>::massDrop() const
{

    tmp<DimensionedField<scalar, volMesh> > tmassDrop
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                this->name() + ":tmassDrop",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimensionedScalar("zero", dimensionSet(1, 0, 0, 0, 0), 0.0  )
        )
    );

    // Average mixture fraction on the droplet surface
    scalarField& sourceField = tmassDrop();
    sourceField += massd_;

    Info<< "  !! ReactingCloudI.H ... massDrop() !!" << endl;  // M.P. Sitte, 12/10/2015

    return tmassDrop;
}

// ************************************************************************* //
