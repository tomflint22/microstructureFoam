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

forAll(nucSiteField,i){
    if (nucSiteField[i] > -100){ // dummy value = -100

        if (
                (Tu[i] >= nucTu[i]) // Nucleation undercooling has been reached.
            &&  (xi[i] < xiThreshold) // xi is below the threshold, to stop nucleation inside existing solid.
            )
        {
            // Cell is to be nucleated, or re nucleated.
            xi[i] = 1.0; // Set global field

            // Check if site hasn't been nucleated yet
            // It will need to be assigned a new field
            if (nucSiteField[i] == -1)
            {
                cellsToBeNuced.append(i);    
            }

            else if (nucSiteField[i] > -1) // It was already nucleated, and assigned a field
            {
                // nucSiteField should be equal to the number of the field e.g. n.17, that was set during the first nucleation event
                label fieldNumber(nucSiteField[i]);

                // Set that field to 1.0
                PopBal[fieldNumber][i] = 1.0;
            }
            
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
            nucSiteField[nucCellI] = maxNiIndex; // Set nucSiteField number. Used for re-nucleation.
            printf("Cell %d was nucleated my proc %d on field n.%d\n", nucCellI, myProc, maxNiIndex);
        }
    }
    // Need to sync this new maxNiIndex with the other processes
    reduce(maxNiIndex, maxOp<label>());   
}

