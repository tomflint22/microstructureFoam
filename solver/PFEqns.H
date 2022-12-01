{



dimensionedScalar gamma_pf("gamma_pf", dimensionSet(0, 0, 0, 0, 0), readScalar(PhaseFieldProperties.lookup("gamma_pf")));

dimensionedScalar sigma_P("sigma_P", dimensionSet(1, 0, -2, 0, 0),readScalar(PhaseFieldProperties.lookup("sigma_P")));
dimensionedScalar sigma_G("sigma_G", dimensionSet(1, 0, -2, 0, 0),readScalar(PhaseFieldProperties.lookup("sigma_G")));
dimensionedScalar lp("lp", dimensionSet(0, 1, 0, 0, 0),readScalar(PhaseFieldProperties.lookup("lp")));
dimensionedScalar lg("lg", dimensionSet(0, 1, 0, 0, 0),readScalar(PhaseFieldProperties.lookup("lg")));

dimensionedScalar deltaf_g("deltaf_g", dimensionSet(0, 0, 0, 0, 0),readScalar(PhaseFieldProperties.lookup("deltaf_g")));
dimensionedScalar deltaf_p("deltaf_p", dimensionSet(0, 0, 0, 0, 0),readScalar(PhaseFieldProperties.lookup("deltaf_p")));

dimensionedScalar D0("D0", dimensionSet(-1, 1, 1, 0, 0),readScalar(PhaseFieldProperties.lookup("D0")));
dimensionedScalar Qg("Qg", dimensionSet(1, 2, -2, 0, -1),readScalar(PhaseFieldProperties.lookup("Qg")));


dimensionedScalar m_p = (3.0/4.0)*((sigma_P)/(deltaf_p*lp));
// dimensionedScalar m_p("m_p", dimensionSet(0, 0, -1, 0, 0), (3.0/4.0)*((sigma_P.value())/(deltaf_p.value()*lp.value())));

dimensionedScalar m_g= (3.0/4.0)*((sigma_G)/(deltaf_g*lg));
dimensionedScalar Kg = sigma_G*lg;
dimensionedScalar Kp = (3.0/4.0)*sigma_P*lp;


// dimensionedScalar Lg("Lg", dimensionSet(-1, 1, 1, 0, 0), readScalar(PhaseFieldProperties.lookup("Lg")));

dimensionedScalar gasconstant("gasconstant",dimensionSet(1, 2, -2, -1, -1),scalar(8.314));













//Anisotropy_Stuff


for(label count = 0; count < N_Ori; count++)
{
    // if (niActive[count])
    // {
volScalarField& ni = PopBal[count];

// Info<<rot[count]<<endl;

gradi[count] = (rot2[count].R() & (((fvc::grad(ni))/(mag(fvc::grad(ni))+deltaN))));
gradi[count].correctBoundaryConditions();
    // }
}



for(label count = 0; count < N_Ori; count++)//find the number of non zero phases at each point in the domain
{

// if (niActive[count])
//     {
volScalarField n_x=(gradi[count].component(0));
volScalarField n_y=(gradi[count].component(1));
volScalarField n_z=(gradi[count].component(2));


phi_Mob[count]=
                            1.0+(Aniso_strength.value()*(


                            (4.0*(Foam::pow(n_x,4.0)+
                            Foam::pow(n_y,4.0)+
                            Foam::pow(n_z,4.0)))-3.0
                            
                            ));
phi_Mob[count].correctBoundaryConditions();
    // }
}






numerator *=0.0;
denominator *=0.0;
Aniso_Mob*=0.0;
for(label i = 0; i < N_Ori; i++)//find the number of non zero phases at each point in the domain
{
    // if (niActive[i])
    // {

volScalarField niav = fvc::average(PopBal[i]);
volScalarField xiav = fvc::average(xi);

numerator+=phi_Mob[i]*niav*niav*(1.0-xi)*(1.0-xi);
denominator+=niav*niav*(1.0-xi)*(1.0-xi);

    // }
}

Aniso_Mob=numerator/(denominator+deltaN2);

// Aniso_Mob=max(min(Aniso_Mob,2.0),0.1);












//Lg = D0*Foam::exp((-Qg)/(gasconstant*T));
Lg = (GBs+Aniso_Mob)*(D0*Foam::exp((-Qg)/(gasconstant*T)))*1.0;

// dimensionedScalar Lp("Lp", dimensionSet(-1, 1, 1, 0, 0), readScalar(PhaseFieldProperties.lookup("Lp")));

// dimensionedScalar Lp = 1.0*D0*Foam::exp((-Qg)/(gasconstant*Tliquidus1));

Lp = Aniso_Mob*(1.0*D0*Foam::exp((-Qg)/(gasconstant*Tliquidus1)));


// Info<<m_p.value()<<endl;
// Info<<m_g.value()<<endl;
// Info<<Kg.value()<<endl;
// Info<<Kp.value()<<endl;
// Info<<Lp.value()<<endl;
// Info<<Lg.value()<<endl;
// Info<<endl;
// Info<<gamma_pf.value()<<endl;
// Info<<sigma_P.value()<<endl;
// Info<<sigma_G.value()<<endl;
// Info<<lp.value()<<endl;
// Info<<lg.value()<<endl;
// Info<<deltaf_g.value()<<endl;
// Info<<deltaf_p.value()<<endl;


Sum_Phases*=0.0;
GBs*=0.0;

for(label count = 0; count < N_Ori; count++)
{
    if (niActive[count])
    {
        volScalarField& ni = PopBal[count];

        GBs+=ni*ni;

        volScalarField test12 = (ni*ni*ni) - ni + 2.0*ni*(1.0-xi)*(1.0-xi);

        Sum_Phases+=ni;
    // volScalarField test12addition = 0.0*ni;//NEED TO SET THIS TO ZERO BETTER
        for(label count2 = 0; count2 < N_Ori; count2++)
        {
            if (niActive[count2])
            {
                if(count != count2){
                    volScalarField& nj = PopBal[count2];
                    // test12addition += nj*nj;
                    test12+=2.0*gamma_pf*ni*nj*nj;
                }
            }
        }
        // test12+=2.0*gamma_pf*ni*test12addition;

        Sum_parts.set(count,new volScalarField(test12));//
        // volScalarField& Su=Sum_parts[count];
    }
}



for(label count = 0; count < N_Ori; count++)
{

    if (niActive[count])
    {
        volScalarField& ni = PopBal[count];
        volScalarField& Su=Sum_parts[count];
        
        solve
        (
            fvm::ddt(ni) 
            ==
            -Lg*m_g*Su
            +Lg*fvm::laplacian(Kg,ni)
        );
    }
}


// <<<<<<< Updated upstream
// phi_PF=0.5*(1.0-(Foam::tanh(1e5*((T/TLiquidus)-1.0))));
// // =======
phi_PF=0.5*(1.0-(Foam::tanh(thetaInterface*((T/TLiquidus)-1.0))));
// >>>>>>> Stashed changes


solve(
    fvm::ddt(xi)
    ==
    -Lp*m_p*((2.0*xi*(1.0-phi_PF))-(2.0*phi_PF*(1.0-xi))) + Lp*m_g*(2.0*(1.0-xi)*GBs)
    +Lp*fvm::laplacian(Kp,xi)
);









}
