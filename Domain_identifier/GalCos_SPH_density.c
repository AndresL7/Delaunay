#include<stdlib.h>
#include "GalCos_variables.h"

/*
  If counter = 0 there is no overlap between 
  If counter = 1, then there is overlap and the node will be opened.
*/

/*
  int GalCos_node_inside_node(float Lsearch,float *pos,float *posnode,float sizenode)
  {
  
  int counter=0;
  float Xleft_node,Yleft_node,Zleft_node,Xleft,Yleft,Zleft;
  float Xright_node,Yright_node,Zright_node,Xright,Yright,Zright;
  
  Xleft_node = posnode[0] - sizenode/2.0;
  Yleft_node = posnode[1] - sizenode/2.0;
  Zleft_node = posnode[2] - sizenode/2.0;
  
  Xright_node = posnode[0] + sizenode/2.0;
  Yright_node = posnode[1] + sizenode/2.0;
  Zright_node = posnode[2] + sizenode/2.0;
  
  Xleft = pos[0] - Lsearch/2.0;
  Yleft = pos[1] - Lsearch/2.0;
  Zleft = pos[2] - Lsearch/2.0;
  
  Xright = pos[0] + Lsearch/2.0;
  Yright = pos[1] + Lsearch/2.0;
  Zright = pos[2] + Lsearch/2.0;
  
  if( (Xleft_node <= Xleft) && (Xright_node >= Xright) )
  {
  
  if( (Yleft_node <= Yleft) && (Yright_node >= Yright) )
  {
	  
  if( (Zleft_node <= Zleft) && (Zright_node >= Zright) )
  {
  counter=1;
  }
  }
  }
  
  return counter;
  
  }
*/

/*
  Verifies if there is some superposition between the segments Xa and
  Xb Xa is defined by the coordinates Xamin and Xamax, as well as for
  the segment Xb.
  
  Routine returns -1 if there is no superposition at all and 1 if
  there is sone superposition.
  
*/

 /*
   int intersect_segments(float Xamin, float Xamax, float Xbmin, float Xbmax)
   {
   
   int flag=-1;
   
   if((Xamin <= Xbmin) && (Xbmin <= Xamax))
   flag=1;
   
   if((Xamin <= Xbmax) && (Xbmax <= Xamax))
   flag=1;
   
   if((Xbmin <= Xamin) && (Xamin <= Xbmax))
   flag=1;
   
   if((Xbmin <= Xamax) && (Xamax <= Xbmax))
   flag=1;
   
   return flag;
   
   }

   
   
   int GalCos_SPHintersect_node(float Lsearch,float *pos,float *posnode,float sizenode)
   {
   
   int counter=0;
   float Xleft_node,Yleft_node,Zleft_node,Xleft,Yleft,Zleft;
   float Xright_node,Yright_node,Zright_node,Xright,Yright,Zright;
   FILE *pf;
   
   Xleft_node = posnode[0] - sizenode/2.0;
   Yleft_node = posnode[1] - sizenode/2.0;
   Zleft_node = posnode[2] - sizenode/2.0;
   
   Xright_node = posnode[0] + sizenode/2.0;
   Yright_node = posnode[1] + sizenode/2.0;
   Zright_node = posnode[2] + sizenode/2.0;
   
   Xleft = pos[0] - Lsearch/2.0;
   Yleft = pos[1] - Lsearch/2.0;
   Zleft = pos[2] - Lsearch/2.0;
  
   Xright = pos[0] + Lsearch/2.0;
   Yright = pos[1] + Lsearch/2.0;
   Zright = pos[2] + Lsearch/2.0;
   
   counter = counter + intersect_segments(Xleft_node,Xright_node,Xleft,Xright);
   counter = counter + intersect_segments(Yleft_node,Yright_node,Yleft,Yright);
   counter = counter + intersect_segments(Zleft_node,Zright_node,Zleft,Zright);
   
   if(counter == 3)
   {
   counter=1;
   }
   else
   counter=0;
  
   return counter;
   
   }
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
		  
		  //dist = distance(Particle[ipart].pos[0],Particle[ipart].pos[1],Particle[ipart].pos[2],
		  //(*p)->pos[0],(*p)->pos[1],(*p)->pos[2]);

		  jpart=(*p)->IDparticle;
		  dist = distance(Particle[ipart].pos[0],Particle[ipart].pos[1],Particle[ipart].pos[2],
				  Particle[jpart].pos[0],Particle[jpart].pos[1],Particle[jpart].pos[2]);
		  
		  if(dist <= Lsearch)
		    {
		      encola(&NGB_cand,jpart);
		      encola_float(&NGB_DISTcand,dist);
		      //printf("Adding neighbour #%d -- %d dist=%f!! %d %f\n",NGB_cand.tail,(*p)->IDparticle,dist,NGB_cand.tail,Lsearch);
		    }
		}
	    }
	  
	  /*
	    if((*p)->tag == 1) 
	    {
	    //printf("Opening node\n");
	    
	    for(i=0; i<(*p)->Nsons; i++)
	    {
	    NGB_driver(ipart,&((*p)->sons[i]),Lsearch);
	    }
	    }
	  */
	  
	}
      
    }
  
}

void call_func(int ipart,TREENODEPTR *p,float Lsearch, int NGB_MAX)
{

  NGB_driver(ipart,&(*p),Lsearch);
  
  //printf("Antes de repetir again %d %f\n",NGB_cand.tail,Lsearch);
  //getchar();

  if(NGB_cand.tail < NGB_MAX)
    {
      //printf("Antes de repetir again %d %f\n",NGB_cand.tail,Lsearch);
      Lsearch = Lsearch + Lsearch/2.0;
      
      reinitialize_cola(&NGB_cand);
      reinitialize_cola_float(&NGB_DISTcand);
      
      
      //NGB_driver(ipart,p,Lsearch);
      call_func(ipart,&(*p),Lsearch,NGB_MAX);
      //printf("Despues de repetir again %d %f\n",NGB_cand.tail,Lsearch);
      //getchar();
      
    }
  //printf("Despues de repetir again %d %f\n",NGB_cand.tail,Lsearch);

  
}

void GalCos_SPH_NGB(int ipart,int *SPH_NGB,float *SPH_NGB_DISTANCES,int NPARTICLES,TREENODEPTR *p,int *inputs,int NGB_MAX)
{
  
  int i;
  float N_dens,Lsearch;
      
  N_dens = 1.0*NPARTICLES/pow((*p)->size,3);
  Lsearch = 0.01*pow(1.0*NGB_MAX/N_dens,1.0/3.0);
    
  load_cola(&NGB_cand);
  load_cola_float(&NGB_DISTcand);
  
  call_func(ipart,&(*p),Lsearch,NGB_MAX);
    
  gsl_fisort(NGB_cand.tail,NGB_DISTcand.inputs,NGB_cand.inputs);
  
  for(i=0; i<NGB_MAX; i++)
    {
      SPH_NGB[i] = NGB_cand.inputs[i];
      SPH_NGB_DISTANCES[i] = NGB_DISTcand.inputs[i];
    }
  
  reinitialize_cola(&NGB_cand);
  reinitialize_cola_float(&NGB_DISTcand);
  
}

