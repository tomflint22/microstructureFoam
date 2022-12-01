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

rho             7966;

    cp  520;
    cpsolid 520.0;
    kappa  10000.0;
	kappasolid  10000.0; 
	Tsolidus 1613;
	Tliquidus 1623.15;
    LatentHeat 1.0;//2.5e5;
    beta    2.0e-5;


// ************************************************************************* //
