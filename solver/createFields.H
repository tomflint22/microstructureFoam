#include "createRDeltaT.H"

dimensionedScalar deltaN=1e-8/pow(average(mesh.V()), 1.0/3.0);

dimensionedScalar pi("pi", dimensionSet(0, 0, 0, 0, 0), constant::mathematical::pi);



Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field Temperature\n" << endl;
volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"


Info<< "Reading phaseProperties\n" << endl;
immiscibleIncompressibleTwoPhaseMixture mixture(U, phi);

autoPtr<twoPhaseChangeModel> phaseChangePtr
(
    twoPhaseChangeModel::New(mixture)
);


// Need to store rho for ddt(rho, U)
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    mixture.rho()
);
rho.oldTime();




#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"


volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + rho*gh
);

pressureReference pressureReference(p, p_rgh, pimple.dict());

if (p_rgh.needReference())
{
    p += dimensionedScalar
         (
             "p",
             p.dimensions(),
             pressureReference.refValue()
             - getRefCellValue(p, pressureReference.refCell())
         );
    p_rgh = p - rho*gh;
}

mesh.schemes().setFluxRequired(p_rgh.name());

volScalarField& alpha1(mixture.alpha1());
mesh.schemes().setFluxRequired(alpha1.name());
















IOdictionary transportProperties_metal
(
    IOobject
    (
        "physicalProperties.metal",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

IOdictionary transportProperties_gas
(
    IOobject
    (
        "physicalProperties.gas",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

IOdictionary phaseProperties
(
    IOobject
    (
        "phaseProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

dimensionedScalar cp1("cp1",dimensionSet(0, 2, -2, -1, 0),transportProperties_metal.lookup("cp"));
dimensionedScalar cp2("cp2",dimensionSet(0, 2, -2, -1, 0),transportProperties_gas.lookup("cp"));

dimensionedScalar cp1solid("cp1solid",dimensionSet(0, 2, -2, -1, 0),transportProperties_metal.lookup("cpsolid"));
dimensionedScalar cp2solid("cp2solid",dimensionSet(0, 2, -2, -1, 0),transportProperties_gas.lookup("cpsolid"));

dimensionedScalar kappa1("kappa1",dimensionSet(1, 1, -3, -1, 0),transportProperties_metal.lookup("kappa"));
dimensionedScalar kappa2("kappa2",dimensionSet(1, 1, -3, -1, 0),transportProperties_gas.lookup("kappa"));

dimensionedScalar kappa1solid("kappa1solid",dimensionSet(1, 1, -3, -1, 0),transportProperties_metal.lookup("kappasolid"));
dimensionedScalar kappa2solid("kappa2solid",dimensionSet(1, 1, -3, -1, 0),transportProperties_gas.lookup("kappasolid"));

dimensionedScalar Tsolidus1("Tsolidus1",dimensionSet(0, 0, 0, 1, 0),transportProperties_metal.lookup("Tsolidus"));
dimensionedScalar Tsolidus2("Tsolidus2",dimensionSet(0, 0, 0, 1, 0),transportProperties_gas.lookup("Tsolidus"));

dimensionedScalar Tliquidus1("Tliquidus1",dimensionSet(0, 0, 0, 1, 0),transportProperties_metal.lookup("Tliquidus"));
dimensionedScalar Tliquidus2("Tliquidus2",dimensionSet(0, 0, 0, 1, 0),transportProperties_gas.lookup("Tliquidus"));

dimensionedScalar LatentHeat1("LatentHeat1",dimensionSet(0, 2, -2, 0, 0),transportProperties_metal.lookup("LatentHeat"));
dimensionedScalar LatentHeat2("LatentHeat2",dimensionSet(0, 2, -2, 0, 0),transportProperties_gas.lookup("LatentHeat"));

dimensionedScalar beta1("beta1",dimensionSet(0, 0, 0, -1, 0),transportProperties_metal.lookup("beta"));
dimensionedScalar beta2("beta2",dimensionSet(0, 0, 0, -1, 0),transportProperties_gas.lookup("beta"));

dimensionedScalar Marangoni_Constant("Marangoni_Constant",dimensionSet(1, 0, -2, -1, 0),phaseProperties.lookup("dsigmadT"));
dimensionedScalar p0("p0",dimensionSet(1,-1,-2,0,0,0,0),phaseProperties.lookup("p0"));//atmospheric pressure
dimensionedScalar Tvap("Tvap",dimensionSet(0, 0, 0, 1, 0),phaseProperties.lookup("Tvap"));
dimensionedScalar Mm("Mm",dimensionSet(1, 0, 0, 0, -1),phaseProperties.lookup("Mm"));
dimensionedScalar LatentHeatVap("LatentHeatVap",dimensionSet(0, 2, -2, 0, 0),phaseProperties.lookup("LatentHeatVap"));
dimensionedScalar R ("R", dimensionSet(1,2,-2,-1,-1,0,0),8.314);


volScalarField cp
(
    IOobject
    (
        "cp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("cp",dimensionSet(0, 2, -2, -1, 0),0.0)
);
cp=alpha1*cp1 + (1.0-alpha1)*cp2;


volScalarField kappa
(
    IOobject
    (
        "kappa",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("kappa",dimensionSet(1, 1, -3, -1, 0),0.0)
);

kappa=alpha1*kappa1 + (1.0-alpha1)*kappa2;


volScalarField TSolidus
(
    IOobject
    (
        "TSolidus",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("TSolidus",dimensionSet(0, 0, 0, 1, 0),0.0)
);
TSolidus=alpha1*Tsolidus1 + (1.0-alpha1)*Tsolidus2;

volScalarField TLiquidus
(
    IOobject
    (
        "TLiquidus",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("TLiquidus",dimensionSet(0, 0, 0, 1, 0),1.0)
);
TLiquidus=alpha1*Tliquidus1 + (1.0-alpha1)*Tliquidus2;

volScalarField Tu
(
    IOobject
    (
        "Tu",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    T - TLiquidus
);


volScalarField DC
(
    IOobject
    (
        "DC",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("DC",dimensionSet(1,-3,-1,0,0),1.0e14)
);



volVectorField n_filtered
(
    IOobject
    (
        "n_filtered",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("n_filtered",dimensionSet(0, 0, 0, 0, 0),vector::zero)
);

volScalarField epsilon1mask
(
    IOobject
    (
        "epsilon1mask",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("epsilon1mask",dimensionSet(0,0,0,0,0),0.0)
);


volVectorField gradT
(
    IOobject
    (
        "gradT",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("gradT",dimensionSet(0, -1, 0, 1, 0),vector::zero)
);

volScalarField LatentHeat
(
    IOobject
    (
        "LatentHeat",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("LatentHeat",dimensionSet(0, 2, -2, 0, 0),0.0)
);
LatentHeat=alpha1*LatentHeat1 + (1.0-alpha1)*LatentHeat2;

volScalarField yDim
(
    IOobject
    (
        "yDim",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("yDim", dimensionSet(0, 1, 0, 0, 0), 1.0)
);

const faceList & ff = mesh.faces();
const pointField & pp = mesh.points();

forAll( mesh.C(), celli)
{
    // vector XYZ = mesh.C()[celli];
    // xcoord[celli]=XYZ.x();
    // zcoord[celli]=XYZ.z();

    const cell & cc = mesh.cells()[celli];
    labelList pLabels(cc.labels(ff));
    pointField pLocal(pLabels.size(), vector::zero);

    forAll (pLabels, pointi)
    {
        pLocal[pointi] = pp[pLabels[pointi]];
    }

    // xDim[celli] = Foam::max(pLocal & vector(1,0,0)) - Foam::min(pLocal & vector(1,0,0));
    yDim[celli] = Foam::max(pLocal & vector(0,1,0)) - Foam::min(pLocal & vector(0,1,0));
    // zDim[celli] = Foam::max(pLocal & vector(0,0,1)) - Foam::min(pLocal & vector(0,0,1));
}

// xDim.correctBoundaryConditions();
yDim.correctBoundaryConditions();
// zDim.correctBoundaryConditions();



volScalarField TRHS
(
    IOobject
    (
        "TRHS",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("TRHS",dimensionSet(1, -1, -3, 0, 0),scalar(0.0))
);
TRHS.write();

volScalarField Tcorr
(
    IOobject
    (
        "Tcorr",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("Tcorr",dimensionSet(0, 0, 0, 1, 0),0.0)
);

volScalarField beta
(
    IOobject
    (
        "beta",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("beta",dimensionSet(0, 0, 0, -1, 0),0.0)
);
beta=alpha1*beta1 + (1.0-alpha1)*beta2;

volScalarField rhok
(
    IOobject
    (
        "rhok",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("rhok",dimensionSet(0,0,0,0,0),0.0)
);

volScalarField Num_divU
(
    IOobject
    (
        "Num_divU",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("Num_divU",dimensionSet(0,0,-1,0,0),0.0)
);


volVectorField Marangoni
(
    IOobject
    (
        "Marangoni",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("Marangoni",dimensionSet(1, -2, -2, 0, 0),vector::zero)
);

volScalarField pVap
(
    IOobject
    (
        "pVap",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("0", dimensionSet(1,-1,-2,0,0,0,0), 0.0)
); // Recoil Pressure

volScalarField Qv
(
    IOobject
    (
        "Qv",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("0", dimensionSet(1,0,-3,0,0,0,0), 0.0)
);// Evaporation Cooling


volScalarField ddte1
(
    IOobject
    (
        "ddte1",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("ddte1",dimensionSet(0,0,-1,0,0),0.0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField damper
(
    IOobject
    (
        "damper",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("damper",dimless,1.0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField thermalDamper
(
    IOobject
    (
        "thermalDamper",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("thermalDamper",dimless,1.0),
    zeroGradientFvPatchScalarField::typeName
);



volScalarField Ray_Number
(
    IOobject
    (
        "Ray_Number",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("Ray_Number",dimensionSet(0, 0, 0, -0, 0),-1.0)
);

volScalarField Deposition
(
    IOobject
    (
        "Deposition",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("Deposition",dimensionSet(1, -1, -3, -0, 0),-1.0)
);

volScalarField Ray_Q
(
    IOobject
    (
        "Ray_Q",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("Ray_Q",dimensionSet(1, 0, -3, 0, 0),scalar(0.0))
);

volScalarField Laser_boundary
(
    IOobject
    (
        "Laser_boundary",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Laser_boundary=fvc::average(Laser_boundary);



volScalarField alpha_filtered
(
    IOobject
    (
        "alpha_filtered",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("alpha_filtered",dimensionSet(0, 0, 0, -0, 0),0.0)
);


volScalarField ErrorTrack
(
    IOobject
    (
        "ErrorTrack",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("ErrorTrack",dimensionSet(0, 0, 0, -0, 0),0.0)
);

volScalarField ErrorTrack2
(
    IOobject
    (
        "ErrorTrack2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("ErrorTrack2",dimensionSet(0, 0, 0, -0, 0),0.0)
);


    volScalarField GBs
    (
        IOobject
        (
            "GBs",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("GBs",dimensionSet(0, 0, 0, 0, 0),0.0)
    );
        volScalarField Sum_Phases
    (
        IOobject
        (
            "Sum_Phases",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("Sum_Phases",dimensionSet(0, 0, 0, 0, 0),0.0)
    );




        volScalarField Lg
    (
        IOobject
        (
            "Lg",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("Lg",dimensionSet(-1, 1, 1, 0, 0),0.0)
    );

        volScalarField Lp
    (
        IOobject
        (
            "Lp",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("Lp",dimensionSet(-1, 1, 1, 0, 0),0.0)
    );



    IOdictionary PhaseFieldProperties
    (
        IOobject
        (
            "PhaseFieldProperties",    // dictionary name
            runTime.constant(),     // dict is found in "constant"
            mesh,                   // registry for the dict
            IOobject::MUST_READ,    // must exist, otherwise failure
            IOobject::NO_WRITE      // dict is only read by the solver
        )
    );

#   include "setNucleationSites.H"
// This will set ni_fieldsNeeded to the number of nucleation sites assigned.

// label N_Ori(readLabel(PhaseFieldProperties.lookup("N_Phases")));
    // Need to get number of tesselated fields.
label N_Seeds(readLabel(PhaseFieldProperties.lookup("N_Seeds")));
label N_Ori = ni_fieldsNeeded + N_Seeds;

Info<<"Number of Phase Field Order Parameters: "<< N_Ori << nl
<< "\tSeeds: " << N_Seeds << nl
<< "\tNuc Sites: " << ni_fieldsNeeded << nl
<<endl;

Info<<"Number of Phase Field Order Parameters: "<<N_Ori<<endl;

// abort();






PtrList<volScalarField> PopBal(N_Ori);
PtrList<volScalarField> Sum_parts(N_Ori);
List<bool> niActive (PopBal.size(),true);

    for(label count = 0; count < N_Ori; count++)
{

    word name_ni ("n." + name(count));

    PopBal.set
    (
        count,
        new volScalarField 
    (
            IOobject
            (
                name_ni,
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,//IOobject::MUST_READ,
                // IOobject::AUTO_WRITE
                IOobject::NO_WRITE
            ),
            mesh,
	    dimensionedScalar("ni",dimensionSet(0, 0, 0, 0, 0),0.0),
	    zeroGradientFvPatchScalarField::typeName
        )
    );
}
//Phase Field Stuff for Non-Conserved Parameters
volScalarField xi
(
    IOobject
    (
        "xi",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("xi",dimensionSet(0, 0, 0, 0, 0),0.0),
    zeroGradientFvPatchScalarField::typeName
);
// xi.write();

// // Melt any cells above liquidus
// forAll( mesh.C(), celli)
// {
//     if(T[celli]>TLiquidus[celli])
//     {
//         for(label count = 0; count < N_Ori; count++)
//         {
//             PopBal[count][celli] = 0.0;
//         }
//     }
//     xi[celli] = 0.0; // TOM: Should this be included
// }

volScalarField grainNum
(
    IOobject
    (
        "grainNum",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("grainNum",dimensionSet(0, 0, 0, 0, 0),-1.0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField maxNiVal
(
    IOobject
    (
        "maxNiVal",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("maxNiVal",dimensionSet(0, 0, 0, 0, 0),0.0),
    zeroGradientFvPatchScalarField::typeName
);

scalar grainNumThreshold( readScalar(PhaseFieldProperties.lookup("grainNumThreshold")));

Switch write_ni_all( PhaseFieldProperties.lookupOrDefault<Switch>("write_ni_all",true) );
Switch write_ni_active( PhaseFieldProperties.lookupOrDefault<Switch>("write_ni_active",true) );


// Set grain num field to indicate which field occupies the cell
grainNum = -1; // Assume that no grains occupy any cell
maxNiVal = 0;

forAll(grainNum,i){
    forAll(PopBal,j){
        if(
            (PopBal[j][i] > grainNumThreshold)
        &&  (PopBal[j][i] > maxNiVal[i])
        )
        {
            maxNiVal[i] = PopBal[j][i];
            grainNum[i] = j;
        }   
    }
}


// Get max field number that currently exists
label maxNiIndex(-1);

forAll(PopBal,i){
    volScalarField& ni = PopBal[i];
    scalar maxVal = gMax(ni);
    Info << "PopBal[" <<i << "] max = " << maxVal << endl;
    if (maxVal > SMALL){ // field is non-zero
        // ni_active[i] = true; // turn on writing of field
        if (i > maxNiIndex)
        {
            maxNiIndex = i;
        }
    }
}

// the maxNiIndex that was found was specific to the processor
// Find the global maxNiIndex by maxOp
reduce(maxNiIndex,maxOp<label>());


volScalarField phi_PF
(
    IOobject
    (
        "phi_PF",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("phi_PF",dimensionSet(0, 0, 0, 0, 0),0.0)
);

volScalarField epsilon1
(
    IOobject
    (
        "epsilon1",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("epsilon1",dimensionSet(0,0,0,0,0),0.0)
);

epsilon1=1.0-xi;

Random randObj2(clock::getTime());


double phi1[N_Ori]={0.0};

double phi2[N_Ori]={0.0};

double phi3[N_Ori]={0.0};

for(label count_j = 0; count_j < N_Ori; count_j++){

phi1[count_j]=randObj2.scalar01()*2.0*pi.value();
// phi2[count_j]=randObj2.scalar01()*2.0*pi.value();
// phi3[count_j]=randObj2.scalar01()*2.0*pi.value();

}

// phi1[0]=0.0;
// phi2[0]=0.0;
// phi3[0]=0.0;

// phi1[1]=0.25*pi.value();//0.0;
// phi2[1]=0.0;//0.25*pi.value();
// phi3[1]=0.0;

// phi1[2]=0.125*pi.value();//0.0;
// phi2[2]=0.0;
// phi3[2]=0.0;



// phi1[N_Ori-1]=0.0;
// phi2[N_Ori-1]=0.0;
// phi3[N_Ori-1]=0.0;



PtrList<volVectorField> gradi(N_Ori);

for(label count = 0; count < N_Ori; count++)
{

    word name_gradi ("gradi." + name(count));

    gradi.set
    (
        count,
        new volVectorField
        (
            IOobject
            (
                name_gradi,
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,//IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedVector("gradi",dimensionSet(0, 0, 0, 0, 0),vector::zero), // NEED TO SET DIMENSIONS TO CHEM POT AFTER I DEBUG
            // zeroGradientFvPatchVectorField::typeName
             fixedValueFvPatchVectorField::typeName
        )
    );
}


tensor rot[N_Ori]={tensor::zero};

for (label i = 0; i < N_Ori; i++)
{


 tensor temprot(
    Foam::cos(phi2[i])*Foam::cos(phi3[i]), 
    (Foam::sin(phi1[i])*Foam::sin(phi2[i])*Foam::cos(phi3[i]))-(Foam::cos(phi1[i])*Foam::sin(phi3[i])),//-Foam::cos(phi3)*Foam::sin(phi2), 
    (Foam::cos(phi1[i])*Foam::sin(phi2[i])*Foam::cos(phi3[i]))+(Foam::sin(phi1[i])*Foam::sin(phi3[i])),//Foam::sin(phi2)*Foam::sin(phi3),
   
    Foam::cos(phi2[i])*Foam::sin(phi3[i]), 
    (Foam::sin(phi1[i])*Foam::sin(phi2[i])*Foam::sin(phi3[i]))+(Foam::cos(phi1[i])*Foam::cos(phi3[i])), 
    (Foam::cos(phi1[i])*Foam::sin(phi2[i])*Foam::sin(phi3[i]))-(Foam::sin(phi1[i])*Foam::cos(phi3[i])),
    
    -Foam::sin(phi2[i]), 
    Foam::sin(phi1[i])*Foam::cos(phi2[i]), 
    Foam::cos(phi1[i])*Foam::cos(phi2[i])
);//https://en.wikipedia.org/wiki/Rotation_matrix - 
//"represents an extrinsic rotation whose (improper) Euler angles are α, β, γ, about axes x, y, z" - phi1 ==alpha, phi2==beta,phi3==gamma

rot[i]=temprot;



}



// Need to check if all procs have the same seed for the random gen
scalar randGenSeed = clock::getTime();
reduce(randGenSeed,minOp<scalar>());
Random randObjq(randGenSeed);

quaternion rot2[N_Ori]={quaternion::zero};

for (label i = 0; i < N_Ori; i++)
{

    scalar X0 = randObjq.scalar01();
    scalar X1 = randObjq.scalar01();
    scalar X2 = randObjq.scalar01();

    scalar theta1 = 2*M_PI*X1;
    scalar theta2 = 2*M_PI*X2;
    scalar s1 = Foam::sin(theta1);
    scalar c1 = Foam::cos(theta1);
    scalar s2 = Foam::sin(theta2);
    scalar c2 = Foam::cos(theta2);
    scalar r1 = Foam::sqrt(1-X0);
    scalar r2 = Foam::sqrt(X0);

    
    // q.normalize();

    // if(mesh.nGeometricD() == 2){
    // quaternion q(s1*r1,vector(c1*r1,s2*r2,c2*r2));
    // rot2[i]=q;
    // }

    // else{
    quaternion q(s1*r1,vector(c1*r1,s2*r2,c2*r2));
    rot2[i]=q;
    // }

    // Info<<q<<endl;

}

quaternion qZero = {quaternion::zero};

Info << "Creating orientation fields\n\tScalar component qw" << endl;
volScalarField qw 
(
    IOobject
    (
        "qw",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    qZero.w()
);

// Info << "\tVector component qv" << endl;
volVectorField qv
(
    IOobject
    (
        "qv",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    qZero.v()
);

// Info<<rot2[0]<<endl;
// Info<<rot2[1]<<endl;
// Info<<rot2[2]<<endl;

// exit(0);


PtrList<volScalarField> phi_Mob(N_Ori);


    for(label count = 0; count < N_Ori; count++)
{

    word name_ni_mob ("n." + name(count)+ "mob");

    phi_Mob.set
    (
        count,
        new volScalarField 
    (
            IOobject
            (
                name_ni_mob,
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,//IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh,
	    dimensionedScalar("phi_Mob",dimensionSet(0, 0, 0, 0, 0),0.0),
	    zeroGradientFvPatchScalarField::typeName
        )
    );
}

volScalarField numerator
(
    IOobject
    (
        "numerator",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("numerator",dimensionSet(0,0,0,0,0),0.0)
);


volScalarField denominator
(
    IOobject
    (
        "denominator",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("denominator",dimensionSet(0,0,0,0,0),0.0)
);

volScalarField Aniso_Mob
(
    IOobject
    (
        "Aniso_Mob",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("Aniso_Mob",dimensionSet(0,0,0,0,0),0.0)
);

dimensionedScalar lengthcorr ("lengthcorr", dimensionSet(0,1,0,0,0,0,0),1.0);
dimensionedScalar deltaN2=lengthcorr*(1e-8/pow(average(mesh.V()), 1.0/3.0));



dimensionedScalar Aniso_strength("Aniso_strength", dimensionSet(0, 0, 0, 0, 0), readScalar(PhaseFieldProperties.lookup("Aniso_strength")));





#include "createAlphaFluxes.H"

// Construct incompressible turbulence model
incompressibleInterPhaseTransportModel turbulence
(
    U,
    phi,
    alphaPhi1,
    mixture
);


#include "createMRF.H"
#include "createFvModels.H"
#include "createFvConstraints.H"


#include "UpdateProps.H"
