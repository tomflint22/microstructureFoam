Info << "Populating grain number field" << endl;

// forAll(PopBal,i){
//     forAll(grainNum,j){
//         if(PopBal[i][j] > 0.5){
//             grainNum[j] = i;
//         }   
//     }
// }

// For every cell
forAll(grainNum, i){
    // Find field number with max value
    label fieldNum(-1);
    scalar maxVal = (xiThreshold); // Set to threshold, ignore if field is less than this

    forAll(PopBal,j){
        if (PopBal[j][i] > maxVal)
        {
            maxVal = PopBal[j][i];
            fieldNum = j;
        }
    }

    grainNum[i] = fieldNum;
}