#include<stdlib.h>
#include "GalCos_variables.h"

/*
  If counter = 0 there is no overlap between 
  If counter = 1, then there is overlap and the node will be opened.
*/

/* Looking for intersection between the search sphere and the node*/

int GalCos_SPHintersect_node(float Lsearch,float *pos,float *posnode,float sizenode)
{
  
  int counter=0;
  float dist,Shalf,a;
  
  dist = distance(pos[0],pos[1],pos[2],posnode[0],posnode[1],posnode[2]);
  Shalf=sizenode/2.0;
  a = sqrt(2.0*Shalf*Shalf);
  
  if( dist <= (Lsearch + 2*a))
    {
      counter=1;
    }
  
  return counter;
}

int GalCos_node_inside_node(float Lsearch,float *pos,float *posnode,float sizenode)
{
  
  int counter=0;
  float dist,Shalf,a;
  
  dist = distance(pos[0],pos[1],pos[2],posnode[0],posnode[1],posnode[2]);
  Shalf = sizenode/2.0;
  a = sqrt(2.0*Shalf*Shalf);
  
  if(Lsearch >= (dist+a))
    {
      counter=1;
    }
  
  
  if(sizenode >= Lsearch)
    {
      
      if( ((pos[0]+Lsearch) <= (posnode[0]+Shalf)) && ((pos[0]-Lsearch) >= (posnode[0]-Shalf)) )
	{
	  if( ((pos[1]+Lsearch) <= (posnode[1]+Shalf)) && ((pos[1]-Lsearch) >= (posnode[1]-Shalf)) )
	    {
	      if( ((pos[2]+Lsearch) <= (posnode[2]+Shalf)) && ((pos[2]-Lsearch) >= (posnode[2]-Shalf)) )
		{
		  counter=1;
		}
	    }
	}
    }
  
  return counter;
  
}

void NGB_driver(int ipart,TREENODEPTR *p,float Lsearch)
{
  
  int flag_sup=0,i,jpart;
  float dist,pos[3];
  
  if((*p) != NULL)
    {
      
      if((*p)->tag != EMPTY_FLAG) 
	{
	  
	  pos[0] = Particle[ipart].pos[0];
	  pos[1] = Particle[ipart].pos[1];
	  pos[2] = Particle[ipart].pos[2];
	  
	  flag_sup = flag_sup + GalCos_SPHintersect_node(Lsearch,pos,(*p)->pos,(*p)->size);
	  flag_sup = flag_sup + GalCos_node_inside_node(Lsearch,pos,(*p)->pos,(*p)->size);
	  
	  if(flag_sup > 0)
	    {
	      
	      if((*p)->tag == 1)
		{
		  
		  for(i=0; i<(*p)->Nsons; i++)
		    {
		      NGB_driver(ipart,&((*p)->sons[i]),Lsearch);
		    }
		}
	      else if(((*p)->tag == 0) && ((*p)->IDparticle != Particle[ipart].id))
		{
		  
		  jpart=(*p)->IDparticle;
		  dist = distance(Particle[ipart].pos[0],Particle[ipart].pos[1],Particle[ipart].pos[2],
				  Particle[jpart].pos[0],Particle[jpart].pos[1],Particle[jpart].pos[2]);
		  
		  if(dist <= Lsearch)
		    {
		      encola(&NGB_cand,jpart);
		      encola_float(&NGB_DISTcand,dist);
		    }
		}
	    }
	  	  
	}
      
    }
  
}

void call_func(int ipart,TREENODEPTR *p,float Lsearch)
{

  NGB_driver(ipart,&(*p),Lsearch);
    
  if(NGB_cand.tail < NGB_MAX)
    {
      
      Lsearch = Lsearch + Lsearch/2.0;
      
      reinitialize_cola(&NGB_cand);
      reinitialize_cola_float(&NGB_DISTcand);
                  
      call_func(ipart,&(*p),Lsearch);
            
    }
    
}

void GalCos_SPH_NGB(int ipart,int *SPH_NGB,float *SPH_NGB_DISTANCES,int NPARTICLES,TREENODEPTR *p,int *inputs)
{
  
  int i;
  float N_dens,Lsearch;
      
  N_dens = 1.0*NPARTICLES/pow((*p)->size,3);
  Lsearch = 0.01*pow(1.0*NGB_MAX/N_dens,1.0/3.0);
    
  load_cola(&NGB_cand);
  load_cola_float(&NGB_DISTcand);
  
  call_func(ipart,&(*p),Lsearch);
    
  gsl_fisort(NGB_cand.tail,NGB_DISTcand.inputs,NGB_cand.inputs);
  
  for(i=0; i<NGB_MAX; i++)
    {
      SPH_NGB[i] = NGB_cand.inputs[i];
      SPH_NGB_DISTANCES[i] = NGB_DISTcand.inputs[i];
    }
  
  reinitialize_cola(&NGB_cand);
  reinitialize_cola_float(&NGB_DISTcand);
  
}

