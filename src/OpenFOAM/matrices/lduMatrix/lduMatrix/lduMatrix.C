/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2024 OpenCFD Ltd.
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

#include "lduMatrix.H"
#include "IOstreams.H"
#include "Switch.H"
#include "objectRegistry.H"
#include "scalarIOField.H"
#include "Time.H"
#include "meshState.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(lduMatrix, 1);
}


const Foam::scalar Foam::lduMatrix::defaultTolerance = 1e-6;

const Foam::Enum
<
    Foam::lduMatrix::normTypes
>
Foam::lduMatrix::normTypesNames_
({
    { normTypes::NO_NORM, "none" },
    { normTypes::DEFAULT_NORM, "default" },
    { normTypes::L1_SCALED_NORM, "L1_scaled" },
});


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::lduMatrix::lduMatrix(const lduMesh& mesh)
:
    lduMesh_(mesh)
{}


Foam::lduMatrix::lduMatrix(const lduMatrix& A)
:
    lduMesh_(A.lduMesh_)
{
    if (A.lowerPtr_)
    {
        lowerPtr_ = std::make_unique<scalarField>(*(A.lowerPtr_));
    }

    if (A.diagPtr_)
    {
        diagPtr_ = std::make_unique<scalarField>(*(A.diagPtr_));
    }

    if (A.upperPtr_)
    {
        upperPtr_ = std::make_unique<scalarField>(*(A.upperPtr_));
    }
}


Foam::lduMatrix::lduMatrix(lduMatrix& A, bool reuse)
:
    lduMesh_(A.lduMesh_)
{
    if (reuse)
    {
        lowerPtr_ = std::move(A.lowerPtr_);
        diagPtr_ = std::move(A.diagPtr_);
        upperPtr_ = std::move(A.upperPtr_);
    }
    else
    {
        if (A.lowerPtr_)
        {
            lowerPtr_ = std::make_unique<scalarField>(*(A.lowerPtr_));
        }

        if (A.diagPtr_)
        {
            diagPtr_ = std::make_unique<scalarField>(*(A.diagPtr_));
        }

        if (A.upperPtr_)
        {
            upperPtr_ = std::make_unique<scalarField>(*(A.upperPtr_));
        }
    }
}


Foam::lduMatrix::lduMatrix(const lduMesh& mesh, Istream& is)
:
    lduMesh_(mesh)
{
    Switch hasLow(is);
    Switch hasDiag(is);
    Switch hasUp(is);

    if (hasLow)
    {
        lowerPtr_ = std::make_unique<scalarField>(is);
    }
    if (hasDiag)
    {
        diagPtr_ = std::make_unique<scalarField>(is);
    }
    if (hasUp)
    {
        upperPtr_ = std::make_unique<scalarField>(is);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::scalarField& Foam::lduMatrix::diag() const
{
    if (!diagPtr_)
    {
        FatalErrorInFunction
            << "diagPtr_ unallocated"
            << abort(FatalError);
    }

    return *diagPtr_;
}


Foam::scalarField& Foam::lduMatrix::diag()
{
    if (!diagPtr_)
    {
        diagPtr_ =
            std::make_unique<scalarField>(lduAddr().size(), Foam::zero{});
    }

    return *diagPtr_;
}


Foam::scalarField& Foam::lduMatrix::diag(const label size)
{
    if (!diagPtr_)
    {
        diagPtr_ = std::make_unique<scalarField>(size, Foam::zero{});
    }

    return *diagPtr_;
}


const Foam::scalarField& Foam::lduMatrix::upper() const
{
    if (upperPtr_)
    {
        return *upperPtr_;
    }
    else
    {
        if (!lowerPtr_)
        {
            FatalErrorInFunction
                << "lowerPtr_ and upperPtr_ unallocated"
                << abort(FatalError);
        }

        return *lowerPtr_;
    }
}


Foam::scalarField& Foam::lduMatrix::upper()
{
    if (!upperPtr_)
    {
        if (lowerPtr_)
        {
            upperPtr_ = std::make_unique<scalarField>(*lowerPtr_);
        }
        else
        {
            upperPtr_ =
                std::make_unique<scalarField>
                (
                    lduAddr().lowerAddr().size(),
                    Foam::zero{}
                );
        }
    }

    return *upperPtr_;
}


Foam::scalarField& Foam::lduMatrix::upper(const label nCoeffs)
{
    if (!upperPtr_)
    {
        if (lowerPtr_)
        {
            upperPtr_ = std::make_unique<scalarField>(*lowerPtr_);
        }
        else
        {
            upperPtr_ = std::make_unique<scalarField>(nCoeffs, Foam::zero{});
        }
    }

    return *upperPtr_;
}


const Foam::scalarField& Foam::lduMatrix::lower() const
{
    if (lowerPtr_)
    {
        return *lowerPtr_;
    }
    else
    {
        if (!upperPtr_)
        {
            FatalErrorInFunction
                << "lowerPtr_ and upperPtr_ unallocated"
                << abort(FatalError);
        }

        return *upperPtr_;
    }
}


Foam::scalarField& Foam::lduMatrix::lower()
{
    if (!lowerPtr_)
    {
        if (upperPtr_)
        {
            lowerPtr_ = std::make_unique<scalarField>(*upperPtr_);
        }
        else
        {
            lowerPtr_ =
                std::make_unique<scalarField>
                (
                    lduAddr().lowerAddr().size(),
                    Foam::zero{}
                );
        }
    }

    return *lowerPtr_;
}


Foam::scalarField& Foam::lduMatrix::lower(const label nCoeffs)
{
    if (!lowerPtr_)
    {
        if (upperPtr_)
        {
            lowerPtr_ = std::make_unique<scalarField>(*upperPtr_);
        }
        else
        {
            lowerPtr_ =
                std::make_unique<scalarField>(nCoeffs, Foam::zero{});
        }
    }

    return *lowerPtr_;
}


void Foam::lduMatrix::setResidualField
(
    const scalarField& residual,
    const word& fieldName,
    const bool initial
) const
{
    if (!mesh().hasDb())
    {
        return;
    }

    scalarIOField* residualPtr =
        mesh().thisDb().getObjectPtr<scalarIOField>
        (
            initial
          ? IOobject::scopedName("initialResidual", fieldName)
          : IOobject::scopedName("residual", fieldName)
        );

    if (residualPtr)
    {
        const auto* dataPtr = mesh().thisDb().findObject<meshState>("data");

        if (dataPtr)
        {
            if (initial && dataPtr->isFirstIteration())
            {
                *residualPtr = residual;
                DebugInfo
                    << "Setting residual field for first solver iteration "
                    << "for solver field: " << fieldName << endl;
            }
        }
        else
        {
            *residualPtr = residual;
            DebugInfo
                << "Setting residual field for solver field "
                << fieldName << endl;
        }
    }
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const lduMatrix& ldum)
{
    Switch hasLow = ldum.hasLower();
    Switch hasDiag = ldum.hasDiag();
    Switch hasUp = ldum.hasUpper();

    os  << hasLow << token::SPACE
        << hasDiag << token::SPACE
        << hasUp << token::SPACE;

    if (hasLow)
    {
        os  << ldum.lower();
    }

    if (hasDiag)
    {
        os  << ldum.diag();
    }

    if (hasUp)
    {
        os  << ldum.upper();
    }

    os.check(FUNCTION_NAME);

    return os;
}


Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const InfoProxy<lduMatrix>& iproxy
)
{
    const auto& ldum = *iproxy;

    Switch hasLow = ldum.hasLower();
    Switch hasDiag = ldum.hasDiag();
    Switch hasUp = ldum.hasUpper();

    os  << "Lower:" << hasLow
        << " Diag:" << hasDiag
        << " Upper:" << hasUp << endl;

    if (hasLow)
    {
        os  << "lower:" << ldum.lower().size() << endl;
    }
    if (hasDiag)
    {
        os  << "diag :" << ldum.diag().size() << endl;
    }
    if (hasUp)
    {
        os  << "upper:" << ldum.upper().size() << endl;
    }


    //if (hasLow)
    //{
    //    os  << "lower contents:" << endl;
    //    forAll(ldum.lower(), i)
    //    {
    //        os  << "i:" << i << "\t" << ldum.lower()[i] << endl;
    //    }
    //    os  << endl;
    //}
    //if (hasDiag)
    //{
    //    os  << "diag contents:" << endl;
    //    forAll(ldum.diag(), i)
    //    {
    //        os  << "i:" << i << "\t" << ldum.diag()[i] << endl;
    //    }
    //    os  << endl;
    //}
    //if (hasUp)
    //{
    //    os  << "upper contents:" << endl;
    //    forAll(ldum.upper(), i)
    //    {
    //        os  << "i:" << i << "\t" << ldum.upper()[i] << endl;
    //    }
    //    os  << endl;
    //}

    os.check(FUNCTION_NAME);

    return os;
}


// ************************************************************************* //
