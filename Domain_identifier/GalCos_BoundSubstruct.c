#include<malloc.h>
#include "GalCos_variables.h"

//#include "GalCos_treebuild_aux.c"

int GalCos_BoundSubstruct(int *inputs,int *sizearray,int *icenter)
{
  
  int i,ipart,*indexek=NULL,FLAG,NPARTICLES,ToErase;
  float vx,vy,vz,*TotalEnergy=NULL;
  float VXcm,VYcm,VZcm,Xcm,Ycm,Zcm;
  
  TREENODEPTR root=NULL;
    
  NPARTICLES=*sizearray;
  
  if(NPARTICLES >= (MINIMUM_NSUBSTRUCT-1))
    {
      
      /* 
	 Building BH tree: Geometric decomposition.
	 Computing properties of nodes: Mass, CM, moments.
      */
      
      root=Build_tree(inputs,NPARTICLES);
      
      /////////////////////////////////////////////////////////////////////
      ///              COMPUTING GRAVITATIONAL POTENTIAL
      /////////////////////////////////////////////////////////////////////
      
      for(i=0; i<NPARTICLES; i++)
	{
	  ipart=inputs[i];
	  Particle[ipart].EP=0;
	  GalCos_tree_force(ipart,root);
	}
      
      GalCos_Free_BHtree(&root);	  
      
      /////////////////////////////////////////////////////////////////////
            
      Xcm = Particle[*icenter].pos[0];
      Ycm = Particle[*icenter].pos[1];
      Zcm = Particle[*icenter].pos[2];
      
      VXcm = Particle[*icenter].vel[0];
      VYcm = Particle[*icenter].vel[1];
      VZcm = Particle[*icenter].vel[2];
      
      ///////////////////**********************
      
      for(i=0; i<NPARTICLES; i++)
	{
	  ipart=inputs[i];
	  
	  Particle[ipart].pos[0] = Particle[ipart].pos[0] - Xcm;
	  Particle[ipart].pos[1] = Particle[ipart].pos[1] - Ycm;
	  Particle[ipart].pos[2] = Particle[ipart].pos[2] - Zcm;
	  
	  Particle[ipart].vel[0] = Particle[ipart].vel[0] - VXcm;
	  Particle[ipart].vel[1] = Particle[ipart].vel[1] - VYcm;
	  Particle[ipart].vel[2] = Particle[ipart].vel[2] - VZcm;
	  
	  vx = Particle[ipart].vel[0]*Particle[ipart].vel[0];
	  vy = Particle[ipart].vel[1]*Particle[ipart].vel[1];
	  vz = Particle[ipart].vel[2]*Particle[ipart].vel[2];
	  
	  Particle[ipart].EK = Particle[ipart].mass*(vx + vy + vz)*0.5;
	  	  
	}
      
      
      TotalEnergy=(float *) malloc((size_t) NPARTICLES*sizeof(float));
      if(TotalEnergy == NULL){
	printf("memory\n");
	exit(0);
      }
      
      indexek=(int *) malloc((size_t) NPARTICLES*sizeof(int));
      if(indexek == NULL){
	printf("memory\n");
	exit(0);
      }
      
      for(i=0; i<NPARTICLES; i++)
	{
	  ipart=inputs[i];
	  TotalEnergy[i] = Particle[ipart].EP + Particle[ipart].EK;
	  indexek[i]=ipart;
	}
      
            
      gsl_fisort(NPARTICLES,TotalEnergy,indexek);
      
      for(i=0; i<NPARTICLES; i++)
	{
	  inputs[i]=indexek[i];
	}
      
      
      FLAG=0;
      for(i=0; i<NPARTICLES; i++)
	{
	  if(TotalEnergy[i] >= 0.0)
	    {
	      FLAG=FLAG+1;
	    }
	}
      
      
      free(TotalEnergy);
      ToErase=(int) (FLAG/4);

      if(FLAG != 0)
	{
	  
	  free(indexek);
	  	  
	  if( FLAG > 20)
	    {
	      (*sizearray)=(*sizearray)-ToErase;
	    }
	  else
	    {
	      (*sizearray)=(*sizearray)-1;
	    }
	  
	  inputs=(int *) realloc(inputs,(size_t) (*sizearray)*sizeof(int));
	  GalCos_BoundSubstruct(inputs,sizearray,icenter);
	}
      else
	{
	  free(indexek);
	  
	  return 0;
	}
      
    }
  
  return 0;
}
