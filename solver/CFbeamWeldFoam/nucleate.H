// reduce(blank0, sumOp<int>()); // Just to sync processes
Info << "\n\nNucleating fields based on temperature. Current number of grains: " << maxNiIndex << endl;
Info << "\tN_Seeds: " << N_Seeds << "\tN_Ori: " << N_Ori << "\tPopBal.size = " << PopBal.size() << endl;

Tu = TLiquidus - Temperature;
Info << "Tu max: " << gMax(Tu) << "\tmin: " << gMin(Tu) << "\tavg: " << gAverage(Tu) << endl;


for(int i = 0; i < nProcs; i++){
    if (i == myProc)
    {
        printf("Proc %d ready for nucleation. N_Ori: %d\n", myProc,N_Ori);
    }
}

// Use a dynamic list on each processor to record cells that should be nucleated.
DynamicList<label> cellsToBeNuced (0);

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