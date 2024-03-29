/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::regionModels::regionModel::mapRegionPatchField
(
    const regionModel& nbrRegion,
    const label regionPatchi,
    const label nbrPatchi,
    const Field<Type>& nbrField,
    const bool flip
) const
{
    const int oldTag = UPstream::incrMsgType();

    const AMIPatchToPatchInterpolation& ami =
        interRegionAMI(nbrRegion, regionPatchi, nbrPatchi, flip);

    tmp<Field<Type>> tresult(ami.interpolateToSource(nbrField));

    UPstream::msgType(oldTag);  // Restore tag

    return tresult;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::regionModels::regionModel::mapRegionPatchField
(
    const regionModel& nbrRegion,
    const word& fieldName,
    const label regionPatchi,
    const bool flip
) const
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    const fvMesh& nbrRegionMesh = nbrRegion.regionMesh();

    if (nbrRegionMesh.foundObject<fieldType>(fieldName))
    {
        const label nbrPatchi = nbrCoupledPatchID(nbrRegion, regionPatchi);

        const int oldTag = UPstream::incrMsgType();

        const AMIPatchToPatchInterpolation& ami =
            interRegionAMI(nbrRegion, regionPatchi, nbrPatchi, flip);

        const fieldType& nbrField =
            nbrRegionMesh.lookupObject<fieldType>(fieldName);

        const Field<Type>& nbrFieldp = nbrField.boundaryField()[nbrPatchi];

        tmp<Field<Type>> tresult(ami.interpolateToSource(nbrFieldp));

        UPstream::msgType(oldTag);  // Restore tag

        return tresult;
    }
    else
    {
        const polyPatch& p = regionMesh().boundaryMesh()[regionPatchi];

        return tmp<Field<Type>>::New(p.size(), Zero);
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::regionModels::regionModel::mapRegionPatchInternalField
(
    const regionModel& nbrRegion,
    const word& fieldName,
    const label regionPatchi,
    const bool flip
) const
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    const fvMesh& nbrRegionMesh = nbrRegion.regionMesh();

    if (nbrRegionMesh.foundObject<fieldType>(fieldName))
    {
        const label nbrPatchi = nbrCoupledPatchID(nbrRegion, regionPatchi);

        const int oldTag = UPstream::incrMsgType();

        const AMIPatchToPatchInterpolation& ami =
            interRegionAMI(nbrRegion, regionPatchi, nbrPatchi, flip);

        const fieldType& nbrField =
            nbrRegionMesh.lookupObject<fieldType>(fieldName);

        const fvPatchField<Type>& nbrFieldp =
            nbrField.boundaryField()[nbrPatchi];

        tmp<Field<Type>> tresult
        (
            ami.interpolateToSource(nbrFieldp.patchInternalField())
        );

        UPstream::msgType(oldTag);  // Restore tag

        return tresult;
    }
    else
    {
        const polyPatch& p = regionMesh().boundaryMesh()[regionPatchi];

        return tmp<Field<Type>>::New(p.size(), Zero);
    }
}


template<class Type>
void Foam::regionModels::regionModel::toPrimary
(
    const label regionPatchi,
    List<Type>& regionField
) const
{
    forAll(intCoupledPatchIDs_, i)
    {
        if (intCoupledPatchIDs_[i] == regionPatchi)
        {
            const mappedPatchBase& mpb =
                refCast<const mappedPatchBase>
                (
                    regionMesh().boundaryMesh()[regionPatchi]
                );
            mpb.reverseDistribute(regionField);
            return;
        }
    }

    FatalErrorInFunction
        << "Region patch ID " << regionPatchi << " not found in region mesh"
        << abort(FatalError);
}


template<class Type>
void Foam::regionModels::regionModel::toRegion
(
    const label regionPatchi,
    List<Type>& primaryField
) const
{
    forAll(intCoupledPatchIDs_, i)
    {
        if (intCoupledPatchIDs_[i] == regionPatchi)
        {
            const mappedPatchBase& mpb =
                refCast<const mappedPatchBase>
                (
                    regionMesh().boundaryMesh()[regionPatchi]
                );
            mpb.distribute(primaryField);
            return;
        }
    }

    FatalErrorInFunction
        << "Region patch ID " << regionPatchi << " not found in region mesh"
        << abort(FatalError);
}


template<class Type, class CombineOp>
void Foam::regionModels::regionModel::toPrimary
(
    const label regionPatchi,
    List<Type>& regionField,
    const CombineOp& cop
) const
{
    forAll(intCoupledPatchIDs_, i)
    {
        if (intCoupledPatchIDs_[i] == regionPatchi)
        {
            const mappedPatchBase& mpb =
                refCast<const mappedPatchBase>
                (
                    regionMesh().boundaryMesh()[regionPatchi]
                );
            mpb.reverseDistribute(regionField, cop);
            return;
        }
    }

    FatalErrorInFunction
        << "Region patch ID " << regionPatchi << " not found in region mesh"
        << abort(FatalError);
}


template<class Type, class CombineOp>
void Foam::regionModels::regionModel::toRegion
(
    const label regionPatchi,
    List<Type>& primaryField,
    const CombineOp& cop
) const
{
    forAll(intCoupledPatchIDs_, i)
    {
        if (intCoupledPatchIDs_[i] == regionPatchi)
        {
            const mappedPatchBase& mpb =
                refCast<const mappedPatchBase>
                (
                    regionMesh().boundaryMesh()[regionPatchi]
                );
            mpb.distribute(primaryField, cop);
            return;
        }
    }

    FatalErrorInFunction
        << "Region patch ID " << regionPatchi << " not found in region mesh"
        << abort(FatalError);
}


// ************************************************************************* //
