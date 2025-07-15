#include "GalCos_variables.h"

float BoxHalf;

float ngb_periodic(float x)
{
  while(x > BoxHalf)
    x-=BoxSize;
  while(x < -BoxHalf)
    x+=BoxSize;
  return x;
}


float periodic_distance(float *posa,float *posb)
{
  
  float dist,Npos[3],Plenght;
  
  dist = distance(posa[0],posa[1],posa[2],posb[0],posb[1],posb[2]);
  
  if(dist <= 0.5*BoxSize)
    return dist;
  
  // Including corrections for periodic boundary conditions /
  
  Plenght = pow(0.5*BoxSize,2);
  
  Npos[0] = posb[0];
  Npos[1] = posb[1];
  Npos[2] = posb[2];
  
  if(pow(posb[0]-posa[0],2) > Plenght)
    {
      if(posa[0] > posb[0])
	Npos[0] = posb[0] + BoxSize;
      else
	Npos[0] = posb[0] - BoxSize;
    }
  
  if(pow(posb[1]-posa[1],2) > Plenght)
    {
      if(posa[1] > posb[1])
	Npos[1] = posb[1] + BoxSize;
      else
	Npos[1] = posb[1] - BoxSize;
    }
  
  if(pow(posb[2]-posa[2],2) > Plenght)
    {
      if(posa[2] > posb[2])
	Npos[2] = posb[2] + BoxSize;
      else
	Npos[2] = posb[2] - BoxSize;
    }
  
  dist = distance(posa[0],posa[1],posa[2],Npos[0],Npos[1],Npos[2]);
  
  return dist;
}


int GalCos_domain(int ipart)
{
  
  int i,counter;
  float dist,mindist,pos[3],ppos[3],x,y,z,Plenght;
  float Xpos,Ypos,Zpos,dx,dy,dz;
  
  counter = EMPTY_FLAG;
  mindist = 1000000*BoxSize;
  
  pos[0] = Particle[ipart].pos[0];
  pos[1] = Particle[ipart].pos[1];
  pos[2] = Particle[ipart].pos[2];
  
  BoxHalf = 0.5*BoxSize;
  
  for(i=0; i<NCLUSTERS; i++)
    {
      
      
      //Xpos = Halos[i].pos[0];
      //Ypos = Halos[i].pos[1];
      //Zpos = Halos[i].pos[2];
      
      /* For periodic boundary corrections */
      /*
	dx = pos[0] - Xpos;
	dy = pos[1] - Ypos;
	dz = pos[2] - Zpos;
	
	dx = ngb_periodic(dx);
	dy = ngb_periodic(dy);
	dz = ngb_periodic(dz);
	
	dist = sqrt(dx*dx + dy*dy + dz*dz);
      */
      dist = distance(Particle[ipart].pos[0],Particle[ipart].pos[1],Particle[ipart].pos[2],
		      Halos[i].pos[0],Halos[i].pos[1],Halos[i].pos[2]);
      
	//dist = periodic_distance(Halos[i].pos,Particle[ipart].pos);
	
      dist = dist/Halos[i].Rvir;
      
      if(dist < mindist)
      {
	  mindist = dist;
	  counter = i;
      }
      
    }
  
  Particle[ipart].Cluster_ID = counter;
  
  Halos[counter].NDomain_particles++;
  Halos[counter].Domain_particles = realloc(Halos[counter].Domain_particles,(size_t) Halos[counter].NDomain_particles*sizeof(int));
  Halos[counter].Domain_particles[Halos[counter].NDomain_particles-1] = Particle[ipart].Oid;
  
  return 0;
}
