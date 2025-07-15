#include "GalCos_variables.h"

int Radial_profiles(void)
{
  
  int k,h,j,*N_h,ipart,NITER,counter_halos;
  float r,dr,*V_s,*V_p,dist,*Ntotal_h,Mparticle,RMAX,*radius,*volumen,counter_volume;
  FILE *pf;
  
  RMAX=50;
  dr=0.5;
  
  NITER = (int) (RMAX/dr);
  
  printf("%d bins\n",NITER);
  
  N_h=(int *) malloc((size_t) NITER*sizeof(int));
  
  V_s=(float *) malloc((size_t) NITER*sizeof(float));
  V_p=(float *) malloc((size_t) NITER*sizeof(float));
  Ntotal_h=(float *) malloc((size_t) NITER*sizeof(float));
  radius=(float *) malloc((size_t) NITER*sizeof(float));
  volumen=(float *) malloc((size_t) NITER*sizeof(float));
  
  Mparticle=Particle[10].mass;
  
  pf=fopen("Voldist.dat","w");
  
  k=0;
  
  for(r=0; r<RMAX; r=r+dr)
    {
      
      radius[k] = r;
      V_s[k] = 0;
      Ntotal_h[k] = 0;
      
      for(h=0; h<NCLUSTERS; h++)
	{
	  
	  if((Halos[h].Nmembers > 75) && (Halos[h].Nmembers < 250))
	    {
	      
	      V_p[k]=0.0;
	      N_h[k]=0;
	      
	      /* CONTRIBUTION FROM DOMAIN ONLY */
	      
	      for(j=0; j<Halos[h].NDomain_particles; j++)
		{
		  
		  ipart = Halos[h].Domain_particles[j];
		  dist = distance(Particle[ipart].pos[0],Particle[ipart].pos[1],Particle[ipart].pos[2],
				  0,0,0);
		  
		  dist=dist/Halos[h].Rvir;
		  //dist=dist/Halos[h].Radius;
		  
		  if((dist >= r) && (dist < (r+dr)))
		    {
		      N_h[k] = N_h[k] + 1;
		      V_p[k] = V_p[k] + Particle[ipart].Volume;
		    }
		  
		}
	      
	      /* CONTRIBUTION FROM HALO ITSELF ONLY */
	      
	      counter_halos=0;
	      counter_volume=0;
	      for(j=0; j<Halos[h].Nmembers; j++)
		{
		  
		  ipart = Halos[h].Halo_particles[j];
		dist = distance(Particle[ipart].pos[0],Particle[ipart].pos[1],Particle[ipart].pos[2],
				0,0,0);
		
		dist=dist/Halos[h].Rvir;
		
		if((dist >= r) && (dist < (r+dr)))
		  {
		    N_h[k] = N_h[k] + 1;
		    V_p[k] = V_p[k] + Particle[ipart].Volume;
		    
		    counter_halos=counter_halos+1;
		    counter_volume=counter_volume+Particle[ipart].Volume;
		  }
		
		}
	      
	      fprintf(pf,"%g %g %d\n",radius[k],V_p[k],N_h[k]);
	      
	      V_s[k] = V_s[k] + V_p[k];
	      Ntotal_h[k] = Ntotal_h[k] + N_h[k]*Mparticle;
	      volumen[k]=(4.0/3.0)*M_PI*(pow((r+dr),3) - pow(r,3));
	    }
	  
	}
      
      k++; 
    }
  
  fclose(pf);
  
  for(k=0; k<NITER; k++)
    printf("%f %f %f %f\n",radius[k],V_s[k],volumen[k],Ntotal_h[k]);
  
  return 0;
  
}
