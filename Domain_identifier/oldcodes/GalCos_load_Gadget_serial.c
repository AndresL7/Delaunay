#include "GalCos_variables.h"

#define MAXPATHLEN 256

struct gadget_head
{
  int      npart[6];
  double   mass[6];
  double   time;
  double   redshift;
  int      flag_sfr;
  int      flag_feedback;
  int      npartTotal[6];
  int      flag_cooling;
  int      num_files;
  double   BoxSize;
  double   Omega0;
  double   OmegaLambda;
  double   HubbleParam; 
  char     fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8];  /* fills to 
							      256 Bytes */
};

struct gadget_head header1;

int GalCos_preload_Gadget(char *infile)
{
  
  char error[2*MAXPATHLEN+14];
  int dummi, i;
  
  
  if((fp_inp=fopen(infile,"r"))==NULL)    
    {
      sprintf(error,"read_gadget cannot open %s",infile);
      exit(FAILURE);
    }
  
  fread(&dummi,sizeof(dummi),1,fp_inp);
  fread(&header1,sizeof(header1),1,fp_inp);
  fread(&dummi,sizeof(dummi),1,fp_inp);
    
  OMEGA_MATTER = header1.Omega0;
  OMEGALAMBDA = header1.OmegaLambda;
  HUBBLEPARAM = header1.HubbleParam;
  COSMIC_TIME = header1.time;
  REDSHIFT = header1.redshift;

  Npart_Total=0;
  for(i=0; i<6; i++)
    {
      if(task == 0) printf(" * Header nall[%d] is: %d \n",i,header1.npartTotal[i]); fflush(stdout);
      Npart_Total = Npart_Total + header1.npartTotal[i];
    }
  
  if(task == 0) printf("\n"); 
  
  for(i=0; i<6; i++)
    { 
      if((header1.npart[i] != 0) && (header1.mass[i] != 0))
	if(task == 0) printf(" * The mass of each particle is %d es %g\n",i,header1.mass[i]); fflush(stdout);
      
      if((header1.npart[i] != 0) && (header1.mass[i] == 0))
	if(task == 0) printf(" * There are individual mases for this particle set %d\n",i); fflush(stdout);
    }     
  
  if(task == 0) 
    {
      printf("\n"); fflush(stdout);
      printf(" * Frame's Time... %g\n",header1.time); fflush(stdout);
      printf(" * Redshift... %g\n",header1.redshift); fflush(stdout);
      printf(" * Flagsfr... %d\n",header1.flag_sfr); fflush(stdout);
      printf(" * Flagfed... %d\n",header1.flag_feedback); fflush(stdout);
      printf(" * Flagcool... %d\n",header1.flag_cooling); fflush(stdout);
      printf(" * numfiles... %d\n",header1.num_files); fflush(stdout);
      printf(" * Boxsize... %g\n",header1.BoxSize); fflush(stdout);
      printf(" * Omega0... %g\n",header1.Omega0); fflush(stdout);
      printf(" * OmageLa... %g\n",header1.OmegaLambda); fflush(stdout);
      printf(" * Hubbleparam... %g\n",header1.HubbleParam); fflush(stdout);
    }

  BoxSize = header1.BoxSize;
  
  return 0;
  
}


int GalCos_load_Gadget(int *unclust, int *sacum, char *infile)
{
  
  int dummi,i,NPH,j,counter,auxint,sorted_groups;
  struct part dp;
  char buff[100];
  FILE *aux_pf;
  
  NPH = Npart_Total;
  
  if(Particle != NULL) free(Particle);
  
  if(NPH != 0)
    {
      
      Particle = (struct part *) malloc((size_t) Npart_Total*sizeof(struct part));
      
      if(Particle == NULL){
	printf("No memory available to load particles \n");
	exit(0);
      }
      
    }
  
  counter=0;
  for(i=0; i<6; i++)
    {
      if((header1.npart[i] != 0) && (header1.mass[i] != 0))
	counter++;
    }
  
  if(counter != 1)
    {
      printf("ERROR Multiple mass distribution or istake reading file\n"); 
      exit(0);
    }
  
  dp.mass = 0.0;
  
  for(i=0; i<6; i++)
    { 
      if((header1.npart[i] != 0) && (header1.mass[i] != 0.0))
	dp.mass=header1.mass[i];
    }  
  
  
  //////////////////// positions
  
  fread(&dummi,sizeof(dummi),1,fp_inp);
  
  sprintf(buff,"%s%s",infile,".grp");
  aux_pf = fopen(buff,"r");
  fscanf(aux_pf,"%d",&auxint);
  
  (*sacum) = 0;
  for(i=0; i<NPH; i++) 
    {
      
      fread(&dp.pos[0],sizeof(float),3,fp_inp);
      fscanf(aux_pf,"%d",&sorted_groups);
      
      if(sorted_groups > (*sacum))
	(*sacum) = sorted_groups;
      if(sorted_groups == 0)
	(*unclust) = (*unclust) + 1;
      
      if(sorted_groups > 0)
	Particle[i].Cluster_ID = sorted_groups-1;
      else
	Particle[i].Cluster_ID = EMPTY_FLAG;
      
      for(j=0; j<3; j++) 
	Particle[i].pos[j] = COSMIC_TIME*dp.pos[j];
      
      Particle[i].mass = dp.mass;
      Particle[i].id = i;
      
    }
  
  fclose(aux_pf);
  fread(&dummi,sizeof(dummi),1,fp_inp);
    
  //////////////////// velocities
  
  fread(&dummi,sizeof(dummi),1,fp_inp);
  
  for(i=0; i<NPH; i++) 
  {
    
    fread(&dp.vel[0],sizeof(float),3,fp_inp);
    
    for(j=0; j<3; j++) 
      Particle[i].vel[j] = sqrt(COSMIC_TIME)*dp.vel[j]; 
    
  }
  
  fread(&dummi,sizeof(dummi),1,fp_inp);
  
  fclose(fp_inp);
    
  return 0;  
  
}
