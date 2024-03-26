// Read fields to be written
IOdictionary writeDict
(
    IOobject
    (
        "writeDict",    // dictionary name
        runTime.constant(),     // dict is found in "constant"
        mesh,                   // registry for the dict
        IOobject::MUST_READ,    // must exist, otherwise failure
        IOobject::NO_WRITE      // dict is only read by the solver
    )
);

// Printing of all fields are off by default, printing all available fields to the user.
wordList volScalarFieldNames = T.db().objectRegistry::sortedNames("volScalarField");
wordList volVectorFieldNames = T.db().objectRegistry::sortedNames("volVectorField");
Info << "Available volScalarFields:" << nl << volScalarFieldNames << endl;
Info << "Available volVectorFields:" << nl << volVectorFieldNames << endl;

DynamicList<word> volScalarFieldsToWrite(writeDict.lookupOrDefault<List<word>>("volScalarFieldsToWrite",List<word>(0)));
DynamicList<word> volVectorFieldsToWrite(writeDict.lookupOrDefault<List<word>>("volVectorFieldsToWrite",List<word>(0)));

// Include simple switch to add default fields
if(writeDict.lookupOrDefault<Switch>("writeDefaultFields",true)){
	Info <<"Forcing writing of default fields: T, grainNum, xi, qw, qv, alpha.metal, U" << endl;
	Info << "This can be switched off by adding writeDefaultFields false; to constant/writeDict" << endl;
    // Need to cycle through existing list and check for existance so that no duplicates
    
    // volScalarFields
    bool T_found(false);
    bool grainNum_found(false);
    bool xi_found(false);
    bool qw_found(false);
    bool alpha_metal_found(false);
    // Feel free to add more, if desired

    forAll(volScalarFieldsToWrite,i){
    	if(volScalarFieldsToWrite[i]=="T") 			{T_found = true;}
    	if(volScalarFieldsToWrite[i]=="grainNum")	{grainNum_found = true;}
    	if(volScalarFieldsToWrite[i]=="xi")			{xi_found = true;}
    	if(volScalarFieldsToWrite[i]=="qw")			{qw_found = true;}
    	if(volScalarFieldsToWrite[i]=="alpha.metal")			{alpha_metal_found = true;}
    }

  	// volVectorFields
	bool U_found(false);
	bool qv_found(false);
	// Feel free to add more, if desired

    forAll(volVectorFieldsToWrite,i){
    	if(volVectorFieldsToWrite[i]=="U") 			{U_found = true;}
		if(volVectorFieldsToWrite[i]=="qv")			{qv_found = true;}
    }

	if(T_found == false){volScalarFieldsToWrite.append("T");}
	if(grainNum_found == false){volScalarFieldsToWrite.append("grainNum");}
	if(xi_found == false){volScalarFieldsToWrite.append("xi");}
	if(qw_found == false){volScalarFieldsToWrite.append("qw");}
	if(alpha_metal_found == false){volScalarFieldsToWrite.append("alpha.metal");}
	
	if(U_found == false){volVectorFieldsToWrite.append("U");}
	if(qv_found == false){volVectorFieldsToWrite.append("qv");}
}

if(volScalarFieldsToWrite.size() == 0){
    Info << "WARNING!! No volScalarFields are to be written."<<nl
    << " Please add desired fields to write from the lists above to volScalarFieldsToWrite list in constant/writeDict" << endl;
}
else{
    Info << "volScalarFields specified to write:" << nl
    << volScalarFieldsToWrite << nl << endl;
}

if(volVectorFieldsToWrite.size() == 0){
    Info << "WARNING!! No volVectorFields are to be written."<<nl
    << " Please add desired fields to write from the lists above to volVectorFieldsToWrite list in constant/writeDict" << endl;
}
else{
    Info << "volVectorFields specified to write:" << nl
    << volVectorFieldsToWrite << nl << endl;
}


Switch write_ni_all( writeDict.lookupOrDefault<Switch>("write_ni_all",true) );
Switch write_ni_active( writeDict.lookupOrDefault<Switch>("write_ni_active",true) );


Info << "write_ni_active: " << write_ni_active << "\twrite_ni_all: " << write_ni_all << endl;