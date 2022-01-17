Info << "Populating grain number field" << endl;

volScalarField grainNum
(
    IOobject
    (
        "grainNum",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    -1.0 // No grains is -1
);

forAll(PopBal,i){
    forAll(grainNum,j){
        if(PopBal[i][j] > 0.5){
            grainNum[j] = i;
        }   
    }
}