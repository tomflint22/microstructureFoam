{


dimensionedScalar Kg("Kg", dimensionSet(0, 2, -1, 0, 0), /*1.0e-6*/readScalar(PhaseFieldProperties.lookup("Kg")));
dimensionedScalar gamma_pf("gamma_pf", dimensionSet(0, 0, 0, 0, 0), /*1.5*/readScalar(PhaseFieldProperties.lookup("gamma_pf")));
dimensionedScalar m_g("m_g", dimensionSet(0, 0, -1, 0, 0), /*4500*/readScalar(PhaseFieldProperties.lookup("m_g")));
dimensionedScalar Lg("Lg", dimensionSet(0, 0, 0, 0, 0), /*1.2*/readScalar(PhaseFieldProperties.lookup("Lg")));

dimensionedScalar m_p("m_p", dimensionSet(0, 0, -1, 0, 0), /*4500.0*/readScalar(PhaseFieldProperties.lookup("m_p")));
dimensionedScalar Lp("Lp", dimensionSet(0, 0, 0, 0, 0), /*1.2*/readScalar(PhaseFieldProperties.lookup("Lp")));
dimensionedScalar Kp("Kp", dimensionSet(0, 2, -1, 0, 0), /*1.0e-6*/readScalar(PhaseFieldProperties.lookup("Kp")));







Sum_Phases*=0.0;

for(label count = 0; count < N_Ori; count++)
{





volScalarField& ni = PopBal[count];
volScalarField test12 = (ni*ni*ni) - ni + 2.0*ni*(1.0-xi)*(1.0-xi);

    Sum_Phases+=ni;

for(label count2 = 0; count2 < N_Ori; count2++)
{
if(count != count2){
volScalarField& nj = PopBal[count2];
test12+=2.0*gamma_pf*ni*nj*nj;
}


}

Sum_parts.set(count,new volScalarField(test12));//
        // volScalarField& Su=Sum_parts[count];



}








GBs*=0.0;

for(label count = 0; count < N_Ori; count++)
{

    volScalarField& ni = PopBal[count];
    volScalarField& Su=Sum_parts[count];

GBs+=ni*ni;
    
solve
                (
                    fvm::ddt(ni) 
                    ==
                    -Lg*m_g*Su
                     +Lg*fvm::laplacian(Kg,ni)
                );
}


phi_PF=0.5*(1.0-(Foam::tanh(1e3*((Temperature/TLiquidus)-1.0))));


solve(
fvm::ddt(xi)
==
-Lp*m_p*((2.0*xi*(1.0-phi_PF))-(2.0*phi_PF*(1.0-xi))) + Lp*m_g*(2.0*(1.0-xi)*GBs)
+Lp*fvm::laplacian(Kp,xi)
);









}
