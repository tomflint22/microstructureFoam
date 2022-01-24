Info<< "Reading Nucleation Parameters from nucleationDict\n" << endl;
IOdictionary nucleationDict
(
    IOobject
    (
        "nucleationDict",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

// Bulk Nucleation Parameters
dimensionedScalar Tu_mean
(
    nucleationDict.lookup("Tu_mean")
);

dimensionedScalar Tu_stdev
(
    nucleationDict.lookup("Tu_stdev")
);

dimensionedScalar n_max // grains per m3
(
    nucleationDict.lookup("n_max")
);
// Check if mesh is 2D, if so scale volumetric n_max
// 0.7236012546 = sqrt(pi/6)
if(mesh.nGeometricD() == 2){
    Info << "Mesh is 2D. Converting 3D n_max = " << n_max.value() << " sites/m3 to 2D n_max = ";
    n_max.value() = Foam::pow( n_max.value() / 0.7236012546 , (2.0/3.0) ) ;
    Info << n_max.value() << " sites/m2" << endl;
}

Info << "Setting nucleation sites with activation undercoolings based on Gaussian distribution" << endl;

// BULK NUCLEATION
// n_max is the maximum density of nucleation sites possible.
scalar meshVolume = gSum(mesh.V());

// Total possible number of nuc sites is meshVolume * n_max
scalar N_max ( /*(int )*/ n_max.value() * meshVolume);

label totalCellsGlobal(mesh.nCells());
reduce(totalCellsGlobal,sumOp<label>());

// Proportion of sites that can nucleate
scalar pSite = scalar(N_max) / totalCellsGlobal;


Info << "Total mesh volume: " << meshVolume 
<< "\tn_max: " << n_max << " sites/m2" 
<< "\tN_max: " << N_max << " sites (" << pSite * 100.0 << "%% of cells)" 
<< endl;

//- 2 options for setting nucleation sites:
//  1:  Set a single cell as a site, with a nucleation undercooling that would be associated with
//      a particular diameter. Assume diameters are based on a Gaussian distribution.
//
//  2:  Set a group of cells as the site, which are nucleated once they become undercooled.
//      This is probably the most accurate, and growth should be staggered such that sites will
//      only grow once the undercooling matches the undercooling needed for thermodynamically 
//      favourable growth for that size circle/sphere.

// Case 1: set point sites with a Gaussian distribution of nucleation undercoolings


// int randSeedOffset = 0;

// Uniform distribution - for selection of nucleation sites
std::random_device rd;  // Will be used to obtain a seed for the random number engine
std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
std::uniform_real_distribution<> rand01(0.0, 1.0);

// Gaussian distribution - for assignment of nucTu
std::default_random_engine generator;
std::normal_distribution<double> Tu_normal(Tu_mean.value(),Tu_stdev.value());

Info << "Random distributions set up. Probability that a cell will be nucleated: " << pSite 
<< " (" << N_max << " out of " << mesh.nCells() << ")" << endl;

/*
// // Create lists for nucleation sites
label listSize(N_max);

labelList nucCells (listSize,-1);
// List<scalar> nucRadius (listSize);
List<bool> nucActive (listSize, false);
List<scalar> nucTu (listSize,scalar(-100));
*/

// //- Using dynamic lists instead
// DynamicList<label> nucCells (0);
// DynamicList<bool> nucActive (0);
// DynamicList<scalar> nucTu(0);

//- Use volScalarFields for all each nucleation variable. 
// It was be a very sparse field, but it's the easiest way to split it up over MPI

Info<< "Creating nucleation undercooling field (nucTu)" << endl;
volScalarField nucTu
(
    IOobject
    (
        "nucTu",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE // Could be manually written instead
    ),
    mesh,
    dimensionedScalar("nucTu",dimensionSet(0, 0, 0, 1, 0),-100.0) // -100 is a dummy value
);


Info << "Marking nucleation sites (nucSite):\t"
    << "-1 = Not a site\t"
    << "0 = Site that is NOT active\t"
    << "1 = Sit that IS active."
    << endl;

volScalarField nucSite
(
    IOobject
    (
        "nucSite",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE // Could be manually written instead
    ),
    mesh,
    dimensionedScalar("nucSite",dimensionSet(0, 0, 0, 1, 0),-1)
);


//- Create counters for number of fields needed
Info << "Creating counters for number of fields needed" << endl;
label ni_fieldsNeeded(0);
label numNucSitesSet(0);
// These will be different counters for each processor.
// Need to sum them from all processors on the master processor,
// Then broadcast it back out, so that each processor is using the
// same number of phase field parameter fields.

// Read threshold value for xi, above which no nucleation sites should be set.
scalar xiThreshold(readScalar(nucleationDict.lookup("xiThreshold")));

// Info << "pSite: " << pSite << "\tN_max: " << N_max << "nCell: " << mesh.nCells() 
// << "\tN_max/nCells" << N_max/mesh.nCells() << endl;
// abort();

forAll(nucTu, i){
    const scalar rand0to1 = rand01(gen);

    if (
            (rand0to1 < pSite)
        &&  (xi[i] < xiThreshold) // Not already solid.
        )
    {
        scalar GaussZ = Tu_normal(generator);
        scalar GaussTu = (GaussZ * Tu_stdev.value()) + Tu_mean.value();

        if (GaussTu < 0)
        {
            GaussTu = 0;
        }

        // nucCells[nucCount] = label(i); // Add cell to list of nucleation sites
        // nucActive[nucCount] = true; // Hasn't been nucleated yet, so it is active
        // nucTu[nucCount] =GaussTu;
        
        //- Dynamic lists
        // nucCells.append( label(i) ); // Add cell to list of nucleation sites
        // nucActive.append( true ); // Hasn't been nucleated yet, so it is active
        // nucTu.append( GaussTu );

        //- Sparsely populated volScalarFields
        // nucCells[i]( label(i) ); // Add cell to list of nucleation sites
        nucSite[i] = 1.0; // 1 = active nuc site
        nucTu[i] = GaussTu;

        ni_fieldsNeeded++;
        numNucSitesSet++;
        // // nucCells[nucCount] = i;
        // nucRadius[nucCount] = nucCount * 1e-6; // Need to change this
        // nucActive[nucCount] = true;
        // nucCount++;
    }
}

// Need to synchronise ni_fieldsNeeded
int nProcs = Pstream::nProcs();
int myProc = Pstream::myProcNo();

// Info << "Current processor number: " << myProc << " out of " << nProcs << endl;

Info << "\nCarrying out reduce operation to sum up all number of ni_needed across processors" << endl;

int ni_fieldsNeeded_myProc = ni_fieldsNeeded;
reduce(ni_fieldsNeeded, sumOp<int>());


for(int i = 0; i < nProcs; i++){
    if (i == myProc)
    {
        printf("Proc %d has %d sites and thinks the total is %d\n", myProc, ni_fieldsNeeded_myProc, ni_fieldsNeeded);
    }
}


Switch checkThatNucSitesExist
(
    nucleationDict.lookup("checkThatNucSitesExist")
);

if (checkThatNucSitesExist)
{
    if (gMax(nucTu) == 0)
    {
        FatalError << "No nucleation sites are set, check nucleation parameters and xiThreshold used in nucleationDict." << nl
        << "Domain could be already solid, if this is intentional, set checkThatNucSitesExist to false." << abort(FatalError);
    }
}

Info << "Number of potential nucleation sites: " << numNucSitesSet << endl;

// Printing nucSite
// Info << "\n--NucSites--  size: " << numNucSitesSet << "\tnCells: " << mesh.nCells() << endl;
// int counter = 0;
// forAll(nucSite, i){
//     if (nucSite[i] == 1)
//     {
//         Info << "Site " << counter << "\tCell:" << i  << "\tActive: " << nucSite[i]  << "\tnucTu: " << nucTu[i] << endl;
//         counter++;
//     }
// }
// Info << "\n\n";

nucTu.write();
// abort();
