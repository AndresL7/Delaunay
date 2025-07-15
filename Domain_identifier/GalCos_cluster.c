#include "GalCos_variables.h"

void GalCos_BuildCluster(struct COLA cola,int Ncluster,int HaloIDCenter)
{
  
  int k;
  
  Halos[Ncluster].ID_CenterHalo=HaloIDCenter;
  Halos[Ncluster].IDcluster=Ncluster;
  Halos[Ncluster].Nmembers=cola.tail;
  Halos[Ncluster].NumberOfSubhalos=0;
  
  Halos[Ncluster].Halo_particles= (int *) malloc((size_t) Halos[Ncluster].Nmembers*sizeof(int));
  
  Halos[Ncluster].Domain_particles=NULL;
  Halos[Ncluster].NDomain_particles=0;

  for(k=0; k<Halos[Ncluster].Nmembers; k++){
    Halos[Ncluster].Halo_particles[k]=cola.inputs[k];       
    Particle[Halos[Ncluster].Halo_particles[k]].Cluster_ID=Halos[Ncluster].IDcluster; 
  }
  
}
