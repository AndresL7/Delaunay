#include "GalCos_variables.h"

void periodic_boundary_corrections(int IDhalo)
{
  
  float Xpos,Ypos,Zpos,x,y,z,Plenght,dist;
  int xcounter,ycounter,zcounter,i,ipart,*indexpot=NULL;
  float *potencial=NULL;
  int NPARTICLES,HaloIDCenter,k,dk,j,jpart;
  
  /*********************************************************/
  /* CHOOSING A CLOSE CANDIDATE TO THE CENTER OF THA HALO */
  
  NPARTICLES = 20;
  
  potencial = (float *) malloc((size_t) NPARTICLES*sizeof(float));
  if(potencial == NULL)
    {
      printf("Problem of memory allocation\n");
      exit(0);
    }
  
  indexpot = (int *) malloc((size_t) NPARTICLES*sizeof(int));
  if(indexpot == NULL)
    {
      printf("Problem of memory allocation\n");
      exit(0);
    }
  
  for(i=0; i<NPARTICLES; i++)
    {
      potencial[i] = 0;
      indexpot[i] = 0;
    }
  
  dk = (int)  (Halos[IDhalo].Nmembers/20);
  k=0;
  for(i=0; i<NPARTICLES; i++)
    {
      
      ipart = Halos[IDhalo].Halo_particles[k];
      Particle[ipart].EP = 0;
      potencial[i] = 0;
      
      for(j=(Halos[IDhalo].Nmembers-1); j>=0; j--)
	{
	  
	  jpart = Halos[IDhalo].Halo_particles[j];
	  
	  if(ipart != jpart)
	    {
	      dist = distance(Particle[ipart].pos[0],Particle[ipart].pos[1],Particle[ipart].pos[2],
			      Particle[jpart].pos[0],Particle[jpart].pos[1],Particle[jpart].pos[2]);
	      
	      Particle[ipart].EP = Particle[ipart].EP - G_INTERNAL_UNITS*(Particle[jpart].mass/GRAV_SOFT)*grav_soft_spline(dist,GRAV_SOFT);
	    }
        }
      
      potencial[i] = Particle[ipart].EP;
      indexpot[i] = ipart;
      Particle[ipart].EP = 0;
      k=k+dk;
      
    }
  
  gsl_fisort(NPARTICLES,potencial,indexpot);
  
  HaloIDCenter = indexpot[0];
  
  free(potencial);
  free(indexpot);
  
  Xpos = Particle[HaloIDCenter].pos[0];
  Ypos = Particle[HaloIDCenter].pos[1];
  Zpos = Particle[HaloIDCenter].pos[2];
  
  Plenght = pow(0.5*BoxSize,2);
  
  for(i=0; i<Halos[IDhalo].Nmembers; i++)
    {
      ipart = Halos[IDhalo].Halo_particles[i];
      
      x = Particle[ipart].pos[0];
      if(pow(Xpos-x,2) > Plenght)
	{
	  
	  if(Xpos > x)
	    Particle[ipart].pos[0] = Particle[ipart].pos[0] + BoxSize;
	  else
	    Particle[ipart].pos[0] = Particle[ipart].pos[0] - BoxSize;
	  
	}
      
      y = Particle[ipart].pos[1];
      if(pow(Ypos-y,2) > Plenght)
	{
	  
	  if(Ypos > y)
	    Particle[ipart].pos[1] = Particle[ipart].pos[1] + BoxSize;
	  else
	    Particle[ipart].pos[1] = Particle[ipart].pos[1] - BoxSize;
	  
	}
      
      z = Particle[ipart].pos[2];
      if(pow(Zpos-z,2) > Plenght)
	{
	  
	  if(Zpos > z)
	    Particle[ipart].pos[2] = Particle[ipart].pos[2] + BoxSize;
	  else
	    Particle[ipart].pos[2] = Particle[ipart].pos[2] - BoxSize;
	  
	}
      
    }
  
}

void Halo_center_mass(struct halo *Haloinp)
{
    
  float VXcm,VYcm,VZcm,M,*radius;
  float Xcenter,Ycenter,Zcenter,VXcenter,VYcenter,VZcenter;
  int k,ipart,*counter_radius;
  
  int FLAG;
  float RhoCrit,TotalMass,density;
  
  radius = (float *) malloc((size_t) Haloinp->Nmembers*sizeof(float));
  counter_radius = (int *) malloc((size_t) Haloinp->Nmembers*sizeof(int));

  M=0;
  for(k=0; k<Haloinp->Nmembers; k++)
    {
      ipart = Haloinp->Halo_particles[k];
      radius[k] = distance(Particle[ipart].pos[0],Particle[ipart].pos[1],Particle[ipart].pos[2],0,0,0);
      counter_radius[k] = ipart;
      M += Particle[ipart].mass;
    }
  
  gsl_fisort(Haloinp->Nmembers,radius,counter_radius);
  
  Haloinp->Radius = radius[Haloinp->Nmembers-1];
  
  /// Para calcular la masa y radio virial
  
  RhoCrit = Virial_criterion();
  
  FLAG = 0;
  TotalMass = 20*Particle[10].mass;
  for(k=20; k<Haloinp->Nmembers; k++)
    {
      
      TotalMass += Particle[10].mass;
      density = (3.0*TotalMass)/(4.0*M_PI*radius[k]*radius[k]*radius[k]);
      
      if(TotalMass <= M/2.0)
	Haloinp->Hmr = radius[k];
      
      if( density <= RhoCrit )
	{
	  Haloinp->Rvir=radius[k];
	  Haloinp->Mvir=TotalMass;
	  FLAG=1;
	  break;
	}
      
    }
  
  
  if(FLAG == 0)
    {
      Haloinp->Mvir = M;
      Haloinp->Rvir = radius[Haloinp->Nmembers-1];
    }
  
  Haloinp->mass = M;
  
  free(radius);
  free(counter_radius);
  
}


//This computes the virial temperature as in galics I, eq. 3.3
/*
  float virial_temperature(float Vc)
  {
  
  float temp,Vcaux;
  
  Vcaux=(Vc*VELOCITY_INTERNAL_UNITS)/100000.0; // the equation to compute the velocity uses the velocity in km/s
  temp=35.9*(Vcaux*Vcaux);
  return temp;
  
  }
*/

/*
  float Total_energy(struct halo *Haloinp)
  {
  
  int i,iparticle;
  float TotalEPotential;
  float TotalEk,TotalEnergy;
  
  /// COMPUTING GRAVITATIONAL POTENTIAL AT PARTICLE I POTENTIAL
  
  TotalEPotential=0.0;
  TotalEk=0.0;
  for(i=0; i<Haloinp->Nmembers; i++)
  {
  iparticle=Haloinp->Halo_particles[i];
  TotalEPotential=TotalEPotential + Particle[iparticle].EP;
  //TotalEk=TotalEk + Particle[iparticle].EK;
  }
  
  TotalEPotential=TotalEPotential/2.0;
  
  Haloinp->TotalEk=TotalEk;
  Haloinp->TotalEP=TotalEPotential;
  TotalEnergy=TotalEk + TotalEPotential;
  
  return TotalEnergy;
  
  }
*/

/*
  float Total_angular_momentum(struct halo *Haloinp)
  {
  
  int ipart,i;
  float Lxtot,Lytot,Lztot,Lx,Ly,Lz,Total_angular;
  
  Lxtot=0.0;
  Lytot=0.0;
  Lztot=0.0;
  
  for(i=0; i<Haloinp->Nmembers; i++)
  {
  
  ipart=Haloinp->Halo_particles[i];
  
  Lx = Particle[ipart].pos[1]*Particle[ipart].vel[2] - Particle[ipart].pos[2]*Particle[ipart].vel[1];
  Ly = Particle[ipart].pos[2]*Particle[ipart].vel[0] - Particle[ipart].pos[0]*Particle[ipart].vel[2];
  Lz = Particle[ipart].pos[0]*Particle[ipart].vel[1] - Particle[ipart].pos[1]*Particle[ipart].vel[0];
      
  Lxtot = Lxtot + Particle[ipart].mass*Lx;
  Lytot = Lytot + Particle[ipart].mass*Ly;
  Lztot = Lztot + Particle[ipart].mass*Lz;
  
  }
  
  Total_angular=sqrt(Lxtot*Lxtot + Lytot*Lytot + Lztot*Lztot);
  
  return Total_angular;

  }
  
  
  float Lambda_parameter(float L, float ETOT, float M)
  {
  
  float LAMDA;
  LAMDA=( L*sqrt(fabs(ETOT)) ) / (G_INTERNAL_UNITS*pow(M,2.5));
  return LAMDA;
  
  }
*/

void GalCos_Halo_prop(struct halo *Haloinp)
{
  
  Halo_center_mass(Haloinp); 
  //Haloinp->Vvir=sqrt(G_INTERNAL_UNITS*Haloinp->Mvir/Haloinp->Rvir);
  //Haloinp->Tvir=virial_temperature(Haloinp->Vvir);
  
  //Haloinp->TotalEnergy=Total_energy(Haloinp);
  //Haloinp->TotalAngularMom=Total_angular_momentum(Haloinp);
  //Haloinp->Lambda_spin=Lambda_parameter(Haloinp->TotalAngularMom,Haloinp->TotalEnergy,Haloinp->mass);
  
}
