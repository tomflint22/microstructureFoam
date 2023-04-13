
cp=(epsilon1*(alpha1*cp1 + (1.0-alpha1)*cp2)) + ((1.0-epsilon1)*(alpha1*cp1solid + (1.0-alpha1)*cp2solid));
kappa=(epsilon1*(alpha1*kappa1 + (1.0-alpha1)*kappa2)) + ((1.0-epsilon1)*(alpha1*kappa1solid + (1.0-alpha1)*kappa2solid));
TSolidus=alpha1*Tsolidus1 + (1.0-alpha1)*Tsolidus2;
TLiquidus=alpha1*Tliquidus1 + (1.0-alpha1)*Tliquidus2;
// LatentHeat=alpha1*LatentHeat1 + (1.0-alpha1)*LatentHeat2;
beta=alpha1*beta1 + (1.0-alpha1)*beta2;

rhok=(1.0-max(epsilon1*(beta)*(T-TSolidus),0.0));
rhok.correctBoundaryConditions();



alpha_filtered=alpha1;

forAll( mesh.C(), celli)
{
    if(alpha1[celli]>0.99)
    {

        alpha_filtered[celli]=1.0;

    }
    else
    {

        if(alpha1[celli]<0.01)
        {

            alpha_filtered[celli]=0.0;

        }
        else
        {
            alpha_filtered[celli]=alpha1[celli];
        }

    }
}
alpha_filtered.correctBoundaryConditions();

LatentHeat=alpha_filtered*LatentHeat1 + (1.0-alpha_filtered)*LatentHeat2;

const volVectorField gradepsilon1(fvc::grad(alpha_filtered));

n_filtered = fvc::grad(alpha_filtered)/((mag(gradepsilon1) + deltaN));

const volScalarField e1temp(fvc::average(epsilon1));



forAll( mesh.C(), celli)
{
    if(e1temp[celli]<=0.95)
    {
        epsilon1mask[celli]=0.0;
    }
    else
    {
        epsilon1mask[celli]=epsilon1[celli];
    }

}

epsilon1mask.correctBoundaryConditions();








//Anisotropy_Stuff


for(label count = 0; count < N_Ori; count++)
{
volScalarField& ni = PopBal[count];

// Info<<rot[count]<<endl;

gradi[count] = (rot2[count].R() & (((fvc::grad(ni))/(mag(fvc::grad(ni))+deltaN))));
gradi[count].correctBoundaryConditions();
}



for(label count = 0; count < N_Ori; count++)//find the number of non zero phases at each point in the domain
{

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

}






numerator *=0.0;
denominator *=0.0;
Aniso_Mob*=0.0;
for(label i = 0; i < N_Ori; i++)//find the number of non zero phases at each point in the domain
{

volScalarField niav = fvc::average(PopBal[i]);
volScalarField xiav = fvc::average(xi);

numerator+=phi_Mob[i]*niav*niav*(1.0-xiav)*(1.0-xiav);
denominator+=niav*niav*(1.0-xiav)*(1.0-xiav);
    // for(label j = 0; j < N_Ori; j++)//find the number of non zero phases at each point in the domain
    // {
    //    numerator+=phi_Mob[i]


    // }

}

Aniso_Mob=numerator/(denominator+deltaN2);
















gh=epsilon1mask*(g & mesh.C());
ghf= fvc::interpolate(epsilon1mask)*(g & mesh.Cf());
