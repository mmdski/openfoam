/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2023 OpenCFD Ltd.
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

Description
    Given a volVectorField and Fluent field identifier, write the field in
    Fluent data format


\*---------------------------------------------------------------------------*/

#include "writeFluentFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::writeFluentField
(
    const volVectorField& phi,
    const label fluentFieldIdentifier,
    Ostream& stream
)
{
    const vectorField& phiInternal = phi;

    // Writing cells
    stream
        << "(300 ("
        << fluentFieldIdentifier << " "  // Field identifier
        << "1 "                  // Zone ID: (cells=1, internal faces=2,
                                 // patch faces=patchi+10)
        << "3 "                  // Number of components (scalar=1, vector=3)
        << "0 0 "                // Unused
        << "1 " << phiInternal.size() // Start and end of list
        << ")(" << nl;

    for (const vector& val : phiInternal)
    {
        stream
            << val.x() << ' ' << val.y() << ' ' << val.z() << nl;
    }

    stream
        << "))" << nl;

    label nWrittenFaces = phiInternal.size();

    // Writing boundary faces
    forAll(phi.boundaryField(), patchi)
    {
        const vectorField& patchPhi = phi.boundaryField()[patchi];

        // Write header
        stream
            << "(300 ("
            << fluentFieldIdentifier << " "  // Field identifier
            << patchi + 10 << " "            // Zone ID: patchi+10
            << "3 "              // Number of components (scalar=1, vector=3)
            << "0 0 "            // Unused
            << nWrittenFaces + 1 << " " << nWrittenFaces + patchPhi.size()
                                 // Start and end of list
            << ")(" << nl;

        nWrittenFaces += patchPhi.size();

        for (const vector& val : patchPhi)
        {
            stream
                << val.x() << ' ' << val.y() << ' ' << val.z() << nl;
        }

        stream
            << "))" << endl;
    }
}


// ************************************************************************* //
