/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

Application
    microstructureFoam

Authors
    Tom Flint, UoM.
    Danny Dreelan, UoM.
    Pratheek Shantraj, UoM.
    Philip Cardiff, UCD.
    Gowthaman Parivendhan, UCD.

Description
    Ray-Tracing heat source implementation with two phase incompressible VoF description of the metallic substrate and shielding gas phase.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "interfaceCompression.H"
#include "CMULES.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "noPhaseChange.H"
#include "incompressibleInterPhaseTransportModel.H"
#include "pimpleControl.H"
#include "pressureReference.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "CorrectPhi.H"
#include "fvcSmooth.H"
#include "findLocalCell.H"

#include <random>    

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    #include "createFieldRefs.H"
    #include "initCorrectPhi.H"
    #include "createUfIfPresent.H"

    if (!LTS)
    {
        #include "CourantNo.H"
        #include "setInitialDeltaT.H"
    }

    boundBox meshBox(mesh.points());
    scalar meshLx = mag( meshBox.max().component(0) - meshBox.min().component(0) );
    scalar meshLy = mag( meshBox.max().component(1) - meshBox.min().component(1) );
    scalar meshLz = mag( meshBox.max().component(2) - meshBox.min().component(2) );

    Info << "Mesh dimenions:" << endl
    << "X: " << meshBox.min().component(0) << "\tto\t" << meshBox.max().component(0) << "\t(" << meshLx << " m)" << endl
    << "Y: " << meshBox.min().component(1) << "\tto\t" << meshBox.max().component(1) << "\t(" << meshLy << " m)" << endl
    << "Z: " << meshBox.min().component(2) << "\tto\t" << meshBox.max().component(2) << "\t(" << meshLz << " m)" << endl
    << endl;

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info<< "\nStarting time loop\n" << endl;

    while (pimple.run(runTime))
    {
        #include "readControls.H"
        #include "readDyMControls.H"

        if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
        {
            #include "CourantNo.H"
            #include "alphaCourantNo.H"
            #include "setDeltaT.H"
        }

        fvModels.preUpdateMesh();

        // Store divU from the previous mesh so that it can be mapped
        // and used in correctPhi to ensure the corrected phi has the
        // same divergence
        tmp<volScalarField> divU;

        if
        (
            correctPhi
         && !isType<twoPhaseChangeModels::noPhaseChange>(phaseChange)
         && mesh.topoChanged()
        )
        {
            // Construct and register divU for correctPhi
            divU = new volScalarField
            (
                "divU0",
                fvc::div(fvc::absolute(phi, U))
            );
        }

        // Update the mesh for topology change, mesh to mesh mapping
        bool topoChanged = mesh.update();

        // Do not apply previous time-step mesh compression flux
        // if the mesh topology changed
        if (topoChanged)
        {
            talphaPhi1Corr0.clear();
        }

        runTime++;

        Info<< "Time = " << runTime.userTimeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            if (pimple.firstPimpleIter() || moveMeshOuterCorrectors)
            {
                if
                (
                    correctPhi
                 && !isType<twoPhaseChangeModels::noPhaseChange>(phaseChange)
                 && !divU.valid()
                )
                {
                    // Construct and register divU for correctPhi
                    divU = new volScalarField
                    (
                        "divU0",
                        fvc::div(fvc::absolute(phi, U))
                    );
                }

                // Move the mesh
                mesh.move();

                if (mesh.changing())
                {
                    gh = (g & mesh.C()) - ghRef;
                    ghf = (g & mesh.Cf()) - ghRef;

                    MRF.update();

                    if (correctPhi)
                    {
                        #include "correctPhi.H"
                    }

                    mixture.correct();

                    if (checkMeshCourantNo)
                    {
                        #include "meshCourantNo.H"
                    }
                }

                divU.clear();
            }

            fvModels.correct();

            surfaceScalarField rhoPhi
            (
                IOobject
                (
                    "rhoPhi",
                    runTime.timeName(),
                    mesh
                ),
                mesh,
                dimensionedScalar(dimMass/dimTime, 0)
            );

            #include "alphaControls.H"
            #include "alphaEqnSubCycle.H"

            #include "UpdateProps.H"
            #include "LaserHS.H"

            turbulence.correctPhasePhi();

            mixture.correct();

            #include "UEqn.H"
            #include "TEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }



            // List<bool> niActive (PopBal.size(),true);
            // label nNiActive(0);

            // forAll(niActive,i){
            //     // if (gMax(PopBal[i]) < SMALL)
            //     if (max(PopBal[i]).value() < SMALL)
            //     {
            //         niActive[i] = false;
            //         nNiActive++;
            //     }
            // }

            // Info << " Before nNiActive: " << nNiActive << endl;
            // FatalError << abort(FatalError);

            #   include "nucleate.H"

            // Info << "Before marking active ni" << endl;

            // Need to mark n.X fields that are active
            // Making a list, will make it on every MPI process, so we should use max rather than gMax
            List<bool> niActive (PopBal.size(),true);
            label nNiActive(0);

            forAll(niActive,i){
                // if (gMax(PopBal[i]) < SMALL)
                niActive[i] = true;
                if (max(PopBal[i]).value() < SMALL)
                {
                    niActive[i] = false;
                    nNiActive++;
                }
            }

            // Info << "After nNiActive: " << nNiActive << endl;
            // FatalError << abort(FatalError);
            // forAll()

            #include "PFEqns.H"

            Info << "Populating grain number field" << endl;

            // if (runTime.writeTime())
            // {
                forAll(PopBal,i){
                    forAll(grainNum,j){
                        // if(PopBal[i][j] > SMALL){
                        if(PopBal[i][j] > 0.5){
                            grainNum[j] = i;
                        }   
                    }
                }
            // }


            // Populating orientation fields
            if (runTime.writeTime())
            {
                forAll(qw,i){
                    label gn = grainNum[i];
                    
                    if (gn < 0)
                    {
                        qw[i] = qZero.w();
                        qv[i] = qZero.v();
                    }

                    else{
                        qw[i] = rot2[gn].w();
                        qv[i] = rot2[gn].v();
                    }
                    
                }
            }



            if (pimple.turbCorr())
            {
                turbulence.correct();
            }
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
