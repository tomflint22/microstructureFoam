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

Info << "Setting nucleation sites with temperatures based on Gaussian distribution" << endl;

// BULK NUCLEATION
// n_max is the maximum density of nucleation sites possible.
scalar meshVolume = gSum(mesh.V());

// Total possible number of nuc sites is meshVolume * n_max
scalar N_max ( /*(int )*/ n_max.value() * meshVolume);

// Proportion of sites that can nucleate
scalar pSite = scalar(N_max) / scalar(mesh.nCells());


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


int randSeedOffset = 0;


Info<< "Creating nucleation undercooling field (TuNuc)" << endl;
volScalarField TuNuc
(
    IOobject
    (
        "TuNuc",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE // Could be manually written instead
    ),
    mesh,
    dimensionedScalar("TuNuc",dimensionSet(0, 0, 0, 1, 0),-100.0) // -100 is a dummy value
);

// Uniform distribution - for selection of nucleation sites
std::random_device rd;  // Will be used to obtain a seed for the random number engine
std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
std::uniform_real_distribution<> rand01(0.0, 1.0);

// Gaussian distribution - for assignment of TuNuc
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

//- Using dynamic lists instead
DynamicList<label> nucCells (0);
DynamicList<bool> nucActive (0);
DynamicList<scalar> nucTu(0);




int nucCount = 0;
int numNucSitesSet(0);

// Info << "pSite: " << pSite << "\tN_max: " << N_max << "nCell: " << mesh.nCells() 
// << "\tN_max/nCells" << N_max/mesh.nCells() << endl;
// abort();

forAll(TuNuc, i){
    const scalar rand0to1 = rand01(gen);

    if (
            (rand0to1 < pSite)
        &&  (xi[i] < SMALL) // Not already solid.
        )
    {
        scalar GaussZ = Tu_normal(generator);
        scalar GaussTu = (GaussZ * Tu_stdev.value()) + Tu_mean.value();

        if (GaussTu < 0)
        {
            GaussTu = 0;
        }

        TuNuc[i] = GaussTu;
        // nucCells[nucCount] = label(i); // Add cell to list of nucleation sites
        // nucActive[nucCount] = true; // Hasn't been nucleated yet, so it is active
        // nucTu[nucCount] =GaussTu;
        

        nucCells.append( label(i) ); // Add cell to list of nucleation sites
        nucActive.append( true ); // Hasn't been nucleated yet, so it is active
        nucTu.append( GaussTu );

        nucCount++;
        numNucSitesSet++;
        // // nucCells[nucCount] = i;
        // nucRadius[nucCount] = nucCount * 1e-6; // Need to change this
        // nucActive[nucCount] = true;
        // nucCount++;
    }
}

if (numNucSitesSet == 0)
{
    FatalError << "No nucleation sites are set, check nucleation parameters!" << abort(FatalError);
}

Info << "Number of potential nucleation sites: " << numNucSitesSet << endl;

// Printing nucSite
Info << "\n--NucSites--  size: " << nucCells.size() << "\tnCells: " << mesh.nCells() << endl;
forAll(nucCells, i){
    Info << "Site " << i << "\tCell:" << nucCells[i] <<"\tactive: " << nucActive[i] << "\tTu: " << TuNuc[nucCells[i]] << endl;
}
Info << "\n\n";
// abort();
