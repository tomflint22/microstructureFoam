Info << "Populating grain number field" << endl;

forAll(PopBal,i){
    forAll(grainNum,j){
        if(PopBal[i][j] > 0.5){
            grainNum[j] = i;
        }   
    }
}