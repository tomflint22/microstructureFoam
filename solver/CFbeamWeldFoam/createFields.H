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

#   include "setNucleationSites.H"

// label N_Ori(nucCount);
#include "createRDeltaT.H"

dimensionedScalar deltaN=1e-8/pow(average(mesh.V()), 1.0/3.0);

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

#include "createPhi.H"

Info<< "Reading transportProperties\n" << endl;
immiscibleIncompressibleTwoPhaseMixture mixture(U, phi);

volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());

const dimensionedScalar& rho1 = mixture.rho1();
const dimensionedScalar& rho2 = mixture.rho2();

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

dictionary phase1 = transportProperties.subDict("phase1");
dictionary phase2 = transportProperties.subDict("phase2");

dimensionedScalar cp1("cp1",dimensionSet(0, 2, -2, -1, 0),phase1.lookup("cp"));
dimensionedScalar cp2("cp2",dimensionSet(0, 2, -2, -1, 0),phase2.lookup("cp"));

dimensionedScalar cp1solid("cp1solid",dimensionSet(0, 2, -2, -1, 0),phase1.lookup("cpsolid"));
dimensionedScalar cp2solid("cp2solid",dimensionSet(0, 2, -2, -1, 0),phase2.lookup("cpsolid"));

dimensionedScalar kappa1("kappa1",dimensionSet(1, 1, -3, -1, 0),phase1.lookup("kappa"));
dimensionedScalar kappa2("kappa2",dimensionSet(1, 1, -3, -1, 0),phase2.lookup("kappa"));

dimensionedScalar kappa1solid("kappa1solid",dimensionSet(1, 1, -3, -1, 0),phase1.lookup("kappasolid"));
dimensionedScalar kappa2solid("kappa2solid",dimensionSet(1, 1, -3, -1, 0),phase2.lookup("kappasolid"));

dimensionedScalar Tsolidus1("Tsolidus1",dimensionSet(0, 0, 0, 1, 0),phase1.lookup("Tsolidus"));
dimensionedScalar Tsolidus2("Tsolidus2",dimensionSet(0, 0, 0, 1, 0),phase2.lookup("Tsolidus"));

dimensionedScalar Tliquidus1("Tliquidus1",dimensionSet(0, 0, 0, 1, 0),phase1.lookup("Tliquidus"));
dimensionedScalar Tliquidus2("Tliquidus2",dimensionSet(0, 0, 0, 1, 0),phase2.lookup("Tliquidus"));

dimensionedScalar LatentHeat1("LatentHeat1",dimensionSet(0, 2, -2, 0, 0),phase1.lookup("LatentHeat"));
dimensionedScalar LatentHeat2("LatentHeat2",dimensionSet(0, 2, -2, 0, 0),phase2.lookup("LatentHeat"));

dimensionedScalar beta1("beta1",dimensionSet(0, 0, 0, -1, 0),phase1.lookup("beta"));
dimensionedScalar beta2("beta2",dimensionSet(0, 0, 0, -1, 0),phase2.lookup("beta"));

dimensionedScalar Marangoni_Constant("Marangoni_Constant",dimensionSet(1, 0, -2, -1, 0),transportProperties.lookup("dsigmadT"));
dimensionedScalar p0("p0",dimensionSet(1,-1,-2,0,0,0,0),transportProperties.lookup("p0"));//atmospheric pressure
dimensionedScalar Tvap("Tvap",dimensionSet(0, 0, 0, 1, 0),transportProperties.lookup("Tvap"));
dimensionedScalar Mm("Mm",dimensionSet(1, 0, 0, 0, -1),transportProperties.lookup("Mm"));
dimensionedScalar LatentHeatVap("LatentHeatVap",dimensionSet(0, 2, -2, 0, 0),transportProperties.lookup("LatentHeatVap"));
dimensionedScalar R ("R", dimensionSet(1,2,-2,-1,-1,0,0),8.314);

dimensionedScalar TSmooth("TSmooth",dimensionSet(0, 0, 0, 1, 0),50.0);


Info<< "Reading field Temperature\n" << endl;
volScalarField Temperature
(
    IOobject
    (
        "Temperature",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


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
TSolidus=alpha1*Tsolidus1 + alpha2*Tsolidus2;

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
TLiquidus=alpha1*Tliquidus1 + alpha2*Tliquidus2;

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
    TLiquidus - Temperature
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




volVectorField nneps1
(
    IOobject
    (
        "nneps1",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
        mesh,
         dimensionedVector("nneps1",dimensionSet(0, 0, 0, 0, 0),vector::zero)
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
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("LatentHeat",dimensionSet(0, 2, -2, 0, 0),0.0)
);
LatentHeat=alpha1*LatentHeat1 + alpha2*LatentHeat2;

volScalarField xcoord
(
    IOobject
    (
        "xcoord",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("xcoord", dimensionSet(0, 1, 0, 0, 0), 0.0)
);

volScalarField zcoord
(
    IOobject
    (
        "zcoord",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
        mesh,
         dimensionedScalar("zcoord", dimensionSet(0, 1, 0, 0, 0), 0.0)
);

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

volScalarField xDim
(
    IOobject
    (
        "xDim",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("xDim", dimensionSet(0, 1, 0, 0, 0), 1.0)
);


volScalarField zDim
(
    IOobject
    (
        "zDim",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zDim", dimensionSet(0, 1, 0, 0, 0), 1.0)
);

const faceList & ff = mesh.faces();
const pointField & pp = mesh.points();

forAll( mesh.C(), celli)
{
    vector XYZ = mesh.C()[celli];
    xcoord[celli]=XYZ.x();
    zcoord[celli]=XYZ.z();

    const cell & cc = mesh.cells()[celli];
    labelList pLabels(cc.labels(ff));
    pointField pLocal(pLabels.size(), vector::zero);

    forAll (pLabels, pointi)
    {
        pLocal[pointi] = pp[pLabels[pointi]];
    }

    xDim[celli] = Foam::max(pLocal & vector(1,0,0)) - Foam::min(pLocal & vector(1,0,0));
    yDim[celli] = Foam::max(pLocal & vector(0,1,0)) - Foam::min(pLocal & vector(0,1,0));
    zDim[celli] = Foam::max(pLocal & vector(0,0,1)) - Foam::min(pLocal & vector(0,0,1));
}

xDim.correctBoundaryConditions();
yDim.correctBoundaryConditions();
zDim.correctBoundaryConditions();

volScalarField cellarea
(
    IOobject
    (
        "cellarea",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
        mesh,
         dimensionedScalar("cellarea", dimensionSet(0, 2, 0, 0, 0), 0.0)
);

forAll( mesh.C(), celli)
{
    cellarea[celli]=(xDim[celli]*zDim[celli]);
}

cellarea.correctBoundaryConditions();

volScalarField sourceTerm
(
    IOobject
    (
        "sourceTerm",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("sourceTerm",dimensionSet(1, -1, -3, 0, 0),scalar(0.0))
);
sourceTerm.write();

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
beta=alpha1*beta1 + alpha2*beta2;

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

// Need to store rho for ddt(rho, U)
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
    ),
    alpha1*rho1 + alpha2*rho2
);
rho.oldTime();

// Mass flux
surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rho)*phi
);

volScalarField ViscousDissipation
(
    IOobject
    (
        "ViscousDissipation",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
        mesh,
        dimensionedScalar("ViscousDissipation",dimensionSet(1, -1, -3, 0, 0),scalar(0.0))
);

volScalarField BeamProfile
(
    IOobject
    (
        "BeamProfile",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
        mesh,
        dimensionedScalar("BeamProfile",dimensionSet(1, -1, -3, 0, 0),scalar(0.0))
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

// Construct incompressible turbulence model
autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, mixture)
);

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

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
    p_rgh = p - rho*gh;
}

mesh.setFluxRequired(p_rgh.name());
mesh.setFluxRequired(alpha1.name());

#include "createMRF.H"
#include "createFvOptions.H"

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

//----------- SET UP PHASE FIELD PARAMETER FIELDS -------------//
// Need to first check how many grains exist, if tesselation is used
// to set initial grain structure.

label N_Ori(ni_fieldsNeeded); // No. fields needed for nucleation sites

// Read number of seeds set for tesselation
label N_Seeds(readLabel(PhaseFieldProperties.lookup("N_Seeds")));


// Check if last seeded field exists
word lastSeededFieldName("n." + name(N_Seeds-1)); // -1 as seeds start at 0
volScalarField tempNi
(
    IOobject
    (
        lastSeededFieldName,
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    -1.0 // dummy value, if field isn't read
);

if (gMax(tempNi) > 0) // Field exists, and has a grain in it
{
    // Assume the lower fields are needed also.
    N_Ori += N_Seeds;
}


// Create / read grainNum field. Reading if present, so that simulation
// can be restarted.
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
    -1.0 // No grains is -1
);

// Getting the maximum value of this field (if it was read) will give
// the number of grains present
label numGrainsInitial(gMax(grainNum)+1); // +1 as grains start at n.0

if (numGrainsInitial > 0) // If grains are present
{
    // Ignore seeded tesselation as this should be accounted for in grainNum
    // Num fields needed = nSites + number of grains already present
    // ni_fieldsNeeded are the nSites set. No sites are set in already solid domain.
    N_Ori = ni_fieldsNeeded + numGrainsInitial; 

}

Info<<"Number of Phase Field Order Parameters: "<<N_Ori<<endl;


reduce(N_Ori,maxOp<label>()); // Make sure all procs have the up to date N_Ori
PtrList<volScalarField> PopBal(N_Ori);
PtrList<volScalarField> Sum_parts(N_Ori);

Info << "After PtrList" << endl;

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
                IOobject::NO_WRITE
            ),
            mesh,
	    dimensionedScalar("ni",dimensionSet(0, 0, 0, 0, 0),0.0),
	    zeroGradientFvPatchScalarField::typeName
        )
    );
    // Info << "Field " <<name_ni << "is set." << endl;
}


// Get max field number that currently exists.
label maxNiIndex(-1);

forAll(PopBal,i){
    scalar maxVal = gMax(PopBal[i]);
    Info << "PopBal[" <<i << "] max = " << maxVal << endl;
    if (maxVal > SMALL){ // field is non-zero
        // ni_active[i] = true; // turn on writing of field
        if (i > maxNiIndex)
        {
            maxNiIndex = i;
        }
    }
}


Info << "\nmaxNiIndex: " << maxNiIndex << nl << endl; 

//Phase Field Stuff for Non-Conserved Parameters
forAll( mesh.C(), celli)
{
    if(Temperature[celli]>TLiquidus[celli]){

        for(label count = 0; count < N_Ori; count++)
        {
            PopBal[count][celli] = 0.0;
        }
    }
}


Info << "Phase field parameter fields are set up." << endl;

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
cp=alpha1*cp1 + alpha2*cp2;

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

kappa=alpha1*kappa1 + alpha2*kappa2;

#include "UpdateProps.H"
