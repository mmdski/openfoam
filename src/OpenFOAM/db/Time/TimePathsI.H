/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenCFD Ltd.
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

inline bool Foam::TimePaths::distributed() const
{
    return distributed_;
}


inline bool Foam::TimePaths::processorCase() const
{
    return processorCase_;
}


inline const Foam::fileName& Foam::TimePaths::rootPath() const
{
    return rootPath_;
}


inline const Foam::fileName& Foam::TimePaths::globalCaseName() const
{
    return globalCaseName_;
}


inline const Foam::fileName& Foam::TimePaths::caseName() const
{
    return case_;
}


inline Foam::fileName& Foam::TimePaths::caseName()
{
    return case_;
}


inline Foam::fileName Foam::TimePaths::path() const
{
    return rootPath()/caseName();
}


inline Foam::fileName Foam::TimePaths::globalPath() const
{
    return rootPath()/globalCaseName();
}


inline Foam::fileName Foam::TimePaths::relativePath
(
    const fileName& input,
    const bool caseTag
) const
{
    return input.relative(globalPath(), caseTag);
}


inline const Foam::word& Foam::TimePaths::constant() const
{
    return constant_;
}


inline const Foam::word& Foam::TimePaths::system() const
{
    return system_;
}


inline Foam::fileName Foam::TimePaths::caseConstant() const
{
    if (processorCase_)
    {
        return ".."/constant();
    }

    return constant();
}


inline Foam::fileName Foam::TimePaths::caseSystem() const
{
    if (processorCase_)
    {
        return ".."/system();
    }

    return system();
}


inline Foam::fileName Foam::TimePaths::constantPath() const
{
    return path()/constant();
}


inline Foam::fileName Foam::TimePaths::systemPath() const
{
    return path()/system();
}


// ************************************************************************* //