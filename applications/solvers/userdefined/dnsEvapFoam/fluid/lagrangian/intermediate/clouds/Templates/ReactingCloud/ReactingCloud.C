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

#include "ReactingCloud.H"

#include "CompositionModel.H"
#include "PhaseChangeModel.H"

// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

template<class CloudType>
void Foam::ReactingCloud<CloudType>::setModels()
{
    compositionModel_.reset
    (
        CompositionModel<ReactingCloud<CloudType>>::New
        (
            this->subModelProperties(),
            *this
        ).ptr()
    );

    phaseChangeModel_.reset
    (
        PhaseChangeModel<ReactingCloud<CloudType>>::New
        (
            this->subModelProperties(),
            *this
        ).ptr()
    );
}


template<class CloudType>
void Foam::ReactingCloud<CloudType>::checkSuppliedComposition
(
    const scalarField& YSupplied,
    const scalarField& Y,
    const word& YName
)
{
    if (YSupplied.size() != Y.size())
    {
        FatalErrorInFunction
            << YName << " supplied, but size is not compatible with "
            << "parcel composition: " << nl << "    "
            << YName << "(" << YSupplied.size() << ") vs required composition "
            << YName << "(" << Y.size() << ")" << nl
            << abort(FatalError);
    }
}


template<class CloudType>
void Foam::ReactingCloud<CloudType>::cloudReset(ReactingCloud<CloudType>& c)
{
    CloudType::cloudReset(c);

    compositionModel_.reset(c.compositionModel_.ptr());
    phaseChangeModel_.reset(c.phaseChangeModel_.ptr());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ReactingCloud<CloudType>::ReactingCloud
(
    const word& cloudName,
    const volScalarField& rho,
    const volVectorField& U,
    const dimensionedVector& g,
    const SLGThermo& thermo,
    bool readFields
)
:
    CloudType(cloudName, rho, U, g, thermo, false),
    reactingCloud(),
    cloudCopyPtr_(nullptr),
    constProps_(this->particleProperties()),
    compositionModel_(nullptr),
    phaseChangeModel_(nullptr),
    rhoTrans_(thermo.carrier().species().size()),
    FsNp_                                           // M.P. Sitte, 04/10/2017
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                this->name() + ":FsNp_",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,  // READ_IF_PRESENT,
                IOobject::NO_WRITE  // AUTO_WRITE
            ),
            this->mesh(),
            dimensionedScalar("zero",  dimensionSet(0, 0, 0, 0, 0), 0.0)
        )
    ),
    GrFsNp_                                           // M.P. Sitte, 04/10/2017
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                this->name() + ":grFsNp_",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,  // READ_IF_PRESENT,
                IOobject::NO_WRITE  // AUTO_WRITE
            ),
            this->mesh(),
            dimensionedScalar("zero",  dimensionSet(0, -1, 0, 0, 0), 0.0)
        )
    ),
    DsNp_                                           // M.P. Sitte, 04/10/2017
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                this->name() + ":DsNp_",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,  // READ_IF_PRESENT,
                IOobject::NO_WRITE  // AUTO_WRITE
            ),
            this->mesh(),
            dimensionedScalar("zero",  dimensionSet(0, 2, -1, 0, 0), 0.0)
        )
    ),
    rdNp_                                           // M.P. Sitte, 28/11/2017
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                this->name() + ":rdNp_",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,  // READ_IF_PRESENT,
                IOobject::NO_WRITE  // AUTO_WRITE
            ),
            this->mesh(),
            dimensionedScalar("zero",  dimensionSet(0, 1, 0, 0, 0), 0.0)
        )
    ),
    ndNp_                                           // M.P. Sitte, 28/11/2017
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                this->name() + ":ndNp_",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,  // READ_IF_PRESENT,
                IOobject::NO_WRITE  // AUTO_WRITE
            ),
            this->mesh(),
            dimensionedScalar("zero",  dimensionSet(0, 0, 0, 0, 0), 0.0)
        )
    ),
    Np_                                             // M.P. Sitte, 04/10/2017
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                this->name() + ":Np_",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,  // READ_IF_PRESENT,
                IOobject::NO_WRITE  // AUTO_WRITE
            ),
            this->mesh(),
            dimensionedScalar("zero",  dimensionSet(0, 0, 0, 0, 0), 0.0)
        )
    ),
    VolNpDt_                                           // M.P. Sitte, 04/10/2017
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                this->name() + ":VolNpDt_",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,  // READ_IF_PRESENT,
                IOobject::NO_WRITE  // AUTO_WRITE
            ),
            this->mesh(),
            dimensionedScalar("zero",  dimensionSet(0, 3, 0, 0, 0), 0.0)
        )
    ),
    NpDt_                                             // M.P. Sitte, 04/10/2017
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                this->name() + ":NpDt_",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,  // READ_IF_PRESENT,
                IOobject::NO_WRITE  // AUTO_WRITE
            ),
            this->mesh(),
            dimensionedScalar("zero",  dimensionSet(0, 0, 0, 0, 0), 0.0)
        )
    ),
    massd_                                             // M.P. Sitte, 04/10/2017
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                this->name() + ":massd_",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,  // READ_IF_PRESENT,
                IOobject::NO_WRITE  // AUTO_WRITE
            ),
            this->mesh(),
            dimensionedScalar("zero",  dimensionSet(1, 0, 0, 0, 0), 0.0)
        )
    )
{
    if (this->solution().active())
    {
        setModels();

        if (readFields)
        {
            parcelType::readFields(*this, this->composition());
            this->deleteLostParticles();
        }
    }

    // Set storage for mass source fields and initialise to zero
    forAll(rhoTrans_, i)
    {
        const word& specieName = thermo.carrier().species()[i];
        rhoTrans_.set
        (
            i,
            new volScalarField::Internal
            (
                IOobject
                (
                    this->name() + ":rhoTrans_" + specieName,
                    this->db().time().timeName(),
                    this->db(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                this->mesh(),
                dimensionedScalar(dimMass, 0)
            )
        );
    }

    if (this->solution().resetSourcesOnStartup())
    {
        resetSourceTerms();
    }
}


template<class CloudType>
Foam::ReactingCloud<CloudType>::ReactingCloud
(
    ReactingCloud<CloudType>& c,
    const word& name
)
:
    CloudType(c, name),
    reactingCloud(),
    cloudCopyPtr_(nullptr),
    constProps_(c.constProps_),
    compositionModel_(c.compositionModel_->clone()),
    phaseChangeModel_(c.phaseChangeModel_->clone()),
    rhoTrans_(c.rhoTrans_.size()),
    FsNp_                                     // M.P. Sitte, 15/10/2015
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                name + ":FsNp_",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            c.FsNp_()
        )
    ),
    GrFsNp_                                     // M.P. Sitte, 15/10/2015
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                name + ":GrFsNp_",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            c.GrFsNp_()
        )
    ),
    DsNp_                                     // M.P. Sitte, 15/10/2015
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                name + ":DsNp_",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            c.DsNp_()
        )
    ),
    rdNp_                                     // M.P. Sitte, 28/11/2015
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                name + ":rdNp_",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            c.rdNp_()
        )
    ),
    ndNp_                                     // M.P. Sitte, 28/11/2015
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                name + ":ndNp_",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            c.ndNp_()
        )
    ),
    Np_                                     // M.P. Sitte, 16/10/2015
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                name + ":Np_",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            c.Np_()
        )
    ),
    VolNpDt_                                     // M.P. Sitte, 16/10/2015
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                name + ":VolNpDt_",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            c.VolNpDt_()
        )
    ),
    NpDt_                                     // M.P. Sitte, 16/10/2015
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                name + ":NpDt_",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            c.NpDt_()
        )
    ),
    massd_                                     // M.P. Sitte, 16/10/2015
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                name + ":massd_",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            c.massd_()
        )
    )
{
    forAll(c.rhoTrans_, i)
    {
        const word& specieName = this->thermo().carrier().species()[i];
        rhoTrans_.set
        (
            i,
            new volScalarField::Internal
            (
                IOobject
                (
                    this->name() + ":rhoTrans_" + specieName,
                    this->db().time().timeName(),
                    this->db(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                c.rhoTrans_[i]
            )
        );
    }
}


template<class CloudType>
Foam::ReactingCloud<CloudType>::ReactingCloud
(
    const fvMesh& mesh,
    const word& name,
    const ReactingCloud<CloudType>& c
)
:
    CloudType(mesh, name, c),
    reactingCloud(),
    cloudCopyPtr_(nullptr),
    constProps_(),
    compositionModel_(c.compositionModel_->clone()),
//    compositionModel_(nullptr),
    phaseChangeModel_(nullptr),
    rhoTrans_(0),
    FsNp_(0),        // M.P. Sitte, 04/10/2017
    GrFsNp_(0),      // M.P. Sitte, 04/10/2017
    DsNp_(0),        // M.P. Sitte, 04/10/2017
    rdNp_(0),        // M.P. Sitte, 28/11/2017
    ndNp_(0),        // M.P. Sitte, 28/11/2017
    Np_(0),          // M.P. Sitte, 04/10/2017
    VolNpDt_(0),     // M.P. Sitte, 04/10/2017
    NpDt_(0),        // M.P. Sitte, 04/10/2017
    massd_(0)         // M.P. Sitte
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ReactingCloud<CloudType>::~ReactingCloud()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::ReactingCloud<CloudType>::setParcelThermoProperties
(
    parcelType& parcel,
    const scalar lagrangianDt
)
{
    CloudType::setParcelThermoProperties(parcel, lagrangianDt);

    parcel.Y() = composition().YMixture0();
}


template<class CloudType>
void Foam::ReactingCloud<CloudType>::checkParcelProperties
(
    parcelType& parcel,
    const scalar lagrangianDt,
    const bool fullyDescribed
)
{
    CloudType::checkParcelProperties(parcel, lagrangianDt, fullyDescribed);

    if (fullyDescribed)
    {
        checkSuppliedComposition
        (
            parcel.Y(),
            composition().YMixture0(),
            "YMixture"
        );
    }

    // derived information - store initial mass
    parcel.mass0() = parcel.mass();
    parcel.d0() = parcel.d();		// E. Weiand - 2020/03/04
}


template<class CloudType>
void Foam::ReactingCloud<CloudType>::storeState()
{
    cloudCopyPtr_.reset
    (
        static_cast<ReactingCloud<CloudType>*>
        (
            clone(this->name() + "Copy").ptr()
        )
    );
}


template<class CloudType>
void Foam::ReactingCloud<CloudType>::restoreState()
{
    cloudReset(cloudCopyPtr_());
    cloudCopyPtr_.clear();
}


template<class CloudType>
void Foam::ReactingCloud<CloudType>::resetSourceTerms()
{
    CloudType::resetSourceTerms();
    forAll(rhoTrans_, i)
    {
        rhoTrans_[i].field() = 0.0;
    }

    // M.P. Sitte, 04/10/2017
    FsNp_().field() = 0.0;
    GrFsNp_().field() = 0.0;
    DsNp_().field() = 0.0;
    rdNp_().field() = 0.0;
    ndNp_().field() = 0.0;
    Np_().field() = 0.0;
    VolNpDt_().field() = 0.0;
    NpDt_().field() = 0.0;
    massd_().field() = 0.0;
}


template<class CloudType>
void Foam::ReactingCloud<CloudType>::relaxSources
(
    const ReactingCloud<CloudType>& cloudOldTime
)
{
    CloudType::relaxSources(cloudOldTime);

    typedef volScalarField::Internal dsfType;

    forAll(rhoTrans_, fieldi)
    {
        dsfType& rhoT = rhoTrans_[fieldi];
        const dsfType& rhoT0 = cloudOldTime.rhoTrans()[fieldi];
        this->relax(rhoT, rhoT0, "rho");
    }
}


template<class CloudType>
void Foam::ReactingCloud<CloudType>::scaleSources()
{
    CloudType::scaleSources();

    typedef volScalarField::Internal dsfType;

    forAll(rhoTrans_, fieldi)
    {
        dsfType& rhoT = rhoTrans_[fieldi];
        this->scale(rhoT, "rho");
    }
}


template<class CloudType>
void Foam::ReactingCloud<CloudType>::evolve()
{
    if (this->solution().canEvolve())
    {
        typename parcelType::trackingData td(*this);

        this->solve(*this, td);
    }
}


template<class CloudType>
void Foam::ReactingCloud<CloudType>::autoMap(const mapPolyMesh& mapper)
{
    Cloud<parcelType>::autoMap(mapper);

    this->updateMesh();
}


template<class CloudType>
void Foam::ReactingCloud<CloudType>::info()
{
    CloudType::info();

    this->phaseChange().info(Info);
}


template<class CloudType>
void Foam::ReactingCloud<CloudType>::writeFields() const
{
    if (compositionModel_.valid())
    {
        CloudType::particleType::writeFields(*this, this->composition());
    }
}


// ************************************************************************* //
