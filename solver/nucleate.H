// reduce(blank0, sumOp<int>()); // Just to sync processes
// Info << "\n\nNucleating fields based on temperature. Current number of grains: " << maxNiIndex << endl;
// Info << "\tN_Seeds: " << N_Seeds << "\tN_Ori: " << N_Ori << "\tPopBal.size = " << PopBal.size() << endl;

Tu = TLiquidus - T;
Info << "Tu max: " << gMax(Tu) << "\tmin: " << gMin(Tu) << "\tavg: " << gAverage(Tu) << endl;

Info << "TLiquidus max: " << gMax(TLiquidus) << "\tmin: " << gMin(TLiquidus) << "\tavg: " << gAverage(TLiquidus) << endl;


for(int i = 0; i < nProcs; i++){
    if (i == myProc)
    {
        printf("Proc %d ready for nucleation. N_Ori: %d\n", myProc,N_Ori);
    }
}

// Need to activate or deactive nucleation sites based on their xi value
forAll(nucSite,i){
    // Deactivate any sites that are now within solid
    if (
            ( mag(nucSite[i]-1.0) < SMALL) // i.e. == 1
        // &&  ( xi[i] > SMALL) // has some solid in it
        // &&  ( xi[i] > xiThreshold) // has some solid in it
        &&  ( grainNum[i] > -0.5) // has some solid in it
        )
    {
        nucSite[i] = 0.0;
    }

    // Activate any sites that are now within liquid (melted)
    if (
            (mag(nucSite[i]) < SMALL) // Currently inactive site (i.e. == 0)
        // &&  (xi[i] < SMALL) // xi indicates liquid
        // &&  (xi[i] < xiThreshold) // xi indicates liquid
        // &&  ( grainNum[i] < -0.5) // has no solid in it
        && (T[i] > TLiquidus[i])
        )
    {
        nucSite[i] = 1.0;
    }
}

// Use a dynamic list on each processor to record cells that should be nucleated.
DynamicList<label> cellsToBeNuced (0);
// label numSitesToBeNuced(0);

forAll(nucSite,i){
    if (nucSite[i] > 0){

        // Check if undercooling has been reached.
        if (Tu[i] >= nucTu[i])
        {
            // Cell is to be nucleated

            cellsToBeNuced.append(i);
            xi[i] = 1.0; // Set global field
            nucSite[i] = 0; // 0 = deactived nucleation site
        }
    }
}


// Now we have a list of cells that should be nucleated on each processor.
// Need to cycle through processors and nucleate all cells on the processor,
// updating maxNiIndex each time. When moving to next processor, pass on maxNiIndex

for(int i = 0; i < nProcs; i++){
    if (i == myProc)
    {
        printf("Proc %d wants to nucleate %d cells\n", myProc, cellsToBeNuced.size());
        forAll(cellsToBeNuced, j){
            label nucCellI = cellsToBeNuced[j];
            maxNiIndex++; // Increment, so that it nucleates it in the next available field
            PopBal[maxNiIndex][nucCellI] = 1.0;
            printf("Cell %d was nucleated my proc %d on field n.%d\n", nucCellI, myProc, maxNiIndex);
        }
    }
    // Need to sync this new maxNiIndex with the other processes
    reduce(maxNiIndex, maxOp<label>());   
}

Info << "end of nucleation " << endl;