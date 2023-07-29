/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Version:  10
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/
FoamFile
{
    format      ascii;
    class       dictionary;
    location    "constant";
    object      physicalProperties.water;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

viscosityModel  constant;

nu              2.5e-06;

rho             8000;

    cp  600;
    cpsolid 600.0;
    kappa  60.0;
	kappasolid  26.0; 
	Tsolidus 1658;
	Tliquidus 1723;
    LatentHeat 2.5e4;//2.5e5;
    beta    2.0e-5;


// ************************************************************************* //
