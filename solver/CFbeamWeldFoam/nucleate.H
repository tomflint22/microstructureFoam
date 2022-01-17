Info << "\n\nNucleating fields based on temperature. Current number of grains: " << maxNiIndex << endl;
Info << "\tN_Seeds: " << N_Seeds << "\tN_Ori: " << N_Ori << "\tPopBal.size = " << PopBal.size() << endl;

Tu = TLiquidus - Temperature;
Info << "Tu max: " << gMax(Tu) << "\tmin: " << gMin(Tu) << "\tavg: " << gAverage(Tu) << endl;

forAll(nucCells,i){
    if (nucActive[i] == true){
        label cenCellI = nucCells[i];

        // Check if undercooling has been reached.
        if (Tu[cenCellI] >= TuNuc[cenCellI])
        {
            // Cell is to be nucleated
            maxNiIndex++; // Increment the number of active PF fields
            PopBal[maxNiIndex][cenCellI] = 1.0; // Set individual field
            xi[cenCellI] = 1.0; // Set global field
            nucActive[i] = false; // Set site to false, to prevent multiple nucleation of the same site
            
            ni_active[maxNiIndex] = true; // Set field to active
        }
        // vector nucCentre = mesh.C()[cenCellI];
        // // scalar R = nucRadius[i];

        // Info << "Nucleating cell: " << highestIndexActive +1 << endl;
        // // Cycle through all cells and and check if the distance is less than R
        // forAll(mesh.C(),j){
        //     scalar d = Foam::mag(nucCentre - mesh.C()[j]);
        //     if (d <= R){
        //         // Cell is to be nucleated
        //         xi[j] = 1.0;
        //         PopBal[highestIndexActive+1][j] = 1.0;
        //     }
        // }    
    }
}