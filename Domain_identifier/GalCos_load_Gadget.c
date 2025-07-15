#include "GalCos_variables.h"
#include "hdf5.h"


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
  char     fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8]; 
};

struct gadget_head header1;

/* Preload Gadget old
int GalCos_preload_Gadget(char *infile)
{
  
  char error[2*MAXPATHLEN+14];
  int dummi, i, counter;
  FILE *fp_inp=NULL;

  hid_t file_id, dataset_id, dataspace_id;
  
  if((fp_inp=fopen(infile,"r"))==NULL)    
    {
      printf("read_gadget cannot open %s",infile);
      exit(FAILURE);
    }
  
  fread(&dummi,sizeof(dummi),1,fp_inp);
  fread(&header1,sizeof(header1),1,fp_inp);
  fread(&dummi,sizeof(dummi),1,fp_inp);
  
  fclose(fp_inp);
  
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
  

  counter=0;
  for(i=0; i<6; i++)
    {
      if((header1.npart[i] != 0) && (header1.mass[i] != 0))
	counter++;
    }
  
  if(counter != 1)
    {
      printf("ERROR Multiple mass distribution or mistake reading file\n"); 
      exit(0);
    }
  
  PARTMASS = 0.0;
  
  for(i=0; i<6; i++)
    { 
      if((header1.npart[i] != 0) && (header1.mass[i] != 0.0))
	PARTMASS = header1.mass[i];
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
      printf(" * Particle mass... %g\n",PARTMASS); fflush(stdout);
    }

  BoxSize = header1.BoxSize;
  
  return 0;
  
}
*/

/////////////Preload Gadget new hdf5 ////////////////////////////////////////////////////////////////////////////
int GalCos_preload_Gadget(char *infile)
{
  hid_t file_id, header_group, attr_id;
  herr_t status;
  hsize_t dims[2];

  printf("%s",infile);
  if((file_id = H5Fopen(infile, H5F_ACC_RDONLY, H5P_DEFAULT)) < 0) {
    printf("Error opening file %s\n", infile);
    exit(FAILURE);
  }

  if((header_group = H5Gopen(file_id, "/Header", H5P_DEFAULT)) < 0) {
    printf("Error opening Header group in %s\n", infile);
    exit(FAILURE);
  }

  int npart[6];
  if((attr_id = H5Aopen(header_group, "NumPart_Total", H5P_DEFAULT)) < 0) {
    printf("Error opening NumPart_Total attribute in %s\n", infile);
    exit(FAILURE);
  }
  status = H5Aread(attr_id, H5T_STD_U32LE, &npart[0]); // number of particles in each category
  H5Aclose(attr_id);
  Npart_Total = npart[1]; //number of particles in the second category (DM particles)  
  if(task == 0) printf("Total number of particles: %d\n", Npart_Total); fflush(stdout);

  if((attr_id = H5Aopen(header_group, "MassTable", H5P_DEFAULT)) < 0) {
    printf("Error opening MassTable attribute in %s\n", infile);
    exit(FAILURE);
  }
  double mass_table[6];
  status = H5Aread(attr_id, H5T_IEEE_F64LE, &mass_table[0]); // mass of particles in each category
  H5Aclose(attr_id);

  PARTMASS = mass_table[1]; // mass of DM particles

  if((attr_id = H5Aopen(header_group, "Time", H5P_DEFAULT)) < 0) {
    printf("Error opening Time attribute in %s\n", infile);
    exit(FAILURE);
  }
  status = H5Aread(attr_id, H5T_IEEE_F64LE, &COSMIC_TIME); // cosmic time
  H5Aclose(attr_id);

  if((attr_id = H5Aopen(header_group, "Redshift", H5P_DEFAULT)) < 0) {
    printf("Error opening Redshift attribute in %s\n", infile);
    exit(FAILURE);
  }
  status = H5Aread(attr_id, H5T_IEEE_F64LE, &REDSHIFT); // redshift
  H5Aclose(attr_id);

  if((attr_id = H5Aopen(header_group, "Flag_Sfr", H5P_DEFAULT)) < 0) {
    printf("Error opening FlagSfr attribute in %s\n", infile);
    exit(FAILURE);
  }
  status = H5Aread(attr_id, H5T_STD_I32LE, &FlagSfr); // star formation flag
  H5Aclose(attr_id);

  if((attr_id = H5Aopen(header_group, "Flag_Feedback", H5P_DEFAULT)) < 0) {
    printf("Error opening FlagFeedback attribute in %s\n", infile);
    exit(FAILURE);
  }
  status = H5Aread(attr_id, H5T_STD_I32LE, &FlagFeedback); // feedback flag
  H5Aclose(attr_id);

  if((attr_id = H5Aopen(header_group, "Flag_Cooling", H5P_DEFAULT)) < 0) {
    printf("Error opening FlagCooling attribute in %s\n", infile);
    exit(FAILURE);
  }
  status = H5Aread(attr_id, H5T_STD_I32LE, &FlagCooling); // cooling flag
  H5Aclose(attr_id);

  if((attr_id = H5Aopen(header_group, "NumFilesPerSnapshot", H5P_DEFAULT)) < 0) {
    printf("Error opening NumFiles attribute in %s\n", infile);
    exit(FAILURE);
  }
  status = H5Aread(attr_id, H5T_STD_I32LE, &NumFiles); // number of files
  H5Aclose(attr_id);

  if((attr_id = H5Aopen(header_group, "BoxSize", H5P_DEFAULT)) < 0) {
    printf("Error opening BoxSize attribute in %s\n", infile);
    exit(FAILURE);
  }
  status = H5Aread(attr_id, H5T_IEEE_F64LE, &BoxSize); // box size
  H5Aclose(attr_id);

  if((attr_id = H5Aopen(header_group, "Omega0", H5P_DEFAULT)) < 0) {
    printf("Error opening Omega0 attribute in %s\n", infile);
    exit(FAILURE);
  }
  status = H5Aread(attr_id, H5T_IEEE_F64LE, &OMEGA_MATTER); // matter density parameter
  H5Aclose(attr_id);

  if((attr_id = H5Aopen(header_group, "OmegaLambda", H5P_DEFAULT)) < 0) {
    printf("Error opening OmegaLambda attribute in %s\n", infile);
    exit(FAILURE);
  }
  status = H5Aread(attr_id, H5T_IEEE_F64LE, &OMEGALAMBDA); // dark energy density parameter
  H5Aclose(attr_id);

  if((attr_id = H5Aopen(header_group, "HubbleParam", H5P_DEFAULT)) < 0) {
    printf("Error opening HubbleParam attribute in %s\n", infile);
    exit(FAILURE);
  }
  status = H5Aread(attr_id, H5T_IEEE_F64LE, &HUBBLEPARAM); // Hubble parameter
  H5Aclose(attr_id);

  if(task == 0) {
    printf("* Frame's Time... %g\n", COSMIC_TIME); fflush(stdout);
    printf("* Redshift... %g\n", REDSHIFT); fflush(stdout);
    printf("* FlagSfr... %d\n", FlagSfr); fflush(stdout);
    printf("* FlagFeedback... %d\n", FlagFeedback); fflush(stdout);
    printf("* FlagCooling... %d\n", FlagCooling); fflush(stdout);
    printf("* NumFiles... %d\n", NumFiles); fflush(stdout);
    printf("* BoxSize... %g\n", BoxSize); fflush(stdout);
    printf("* Omega0... %g\n", OMEGA_MATTER); fflush(stdout);
    printf("* OmegaLambda... %g\n", OMEGALAMBDA); fflush(stdout);
    printf("* HubbleParam... %g\n", HUBBLEPARAM); fflush(stdout);
    printf("* Particle mass... %g\n", PARTMASS); fflush(stdout);
    printf("====================================================\n");
  }
  H5Gclose(header_group);
  H5Fclose(file_id);

  return 0;
}

// Preload Catalog Gadget
int GalCos_preload_Catalog(char *cat_file)
{
  hid_t file_id, header_group, attr_id;
  herr_t status;
  hid_t group_id;

  char buff[100];
  int i, j, k, sorted_groups, auxint, counter, maxsend, maxrecv;
  struct part dp;
  FILE *aux_pf = NULL;
  FILE *fp_inp = NULL;

  // Open the catalog file
  if((file_id = H5Fopen(cat_file, H5F_ACC_RDONLY, H5P_DEFAULT)) < 0) {
    printf("Error opening catalog file %s\n", cat_file);
    exit(FAILURE);
  }

  // Open the group containing the halo data
  if((group_id = H5Gopen(file_id, "/Group", H5P_DEFAULT)) < 0) {
    printf("Error opening Group in %s\n", cat_file);
    exit(FAILURE);
  }

  //Header
  if((header_group = H5Gopen(file_id, "/Header", H5P_DEFAULT)) < 0) {
    printf("Error opening Header group in %s\n", cat_file);
    exit(FAILURE);
  }

  if((attr_id = H5Aopen(header_group, "Ngroups_Total", H5P_DEFAULT)) < 0) {
    printf("Error opening NgroupsTotal attribute in %s\n", cat_file);
    exit(FAILURE);
  }
  status = H5Aread(attr_id, H5T_NATIVE_INT, &NCLUSTERS); // total number of clusters
  H5Aclose(attr_id);

  if((attr_id = H5Aopen(header_group, "Nids_ThisFile", H5P_DEFAULT)) < 0) {
    printf("Error opening NidsThisFile attribute in %s\n", cat_file);
    exit(FAILURE);
  }
  status = H5Aread(attr_id, H5T_NATIVE_INT, &Npart_clustered); // number of clustered particles
  H5Aclose(attr_id);

  Npart_unclustered = Npart_Total - Npart_clustered; // number of unclustered particles

  if(task == 0) {
    // printf("Total number of clusters: %d\n", NCLUSTERS); fflush(stdout);
    // printf("Number of clustered particles: %d\n", Npart_clustered); fflush(stdout);
    printf("Number of unclustered particles: %d\n", Npart_unclustered); fflush(stdout);
    printf("====================================================\n");
  }

  H5Gclose(group_id);
  H5Fclose(file_id);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Load Gadget New
int GalCos_load_Gadget(int *unclust, int *sacum, char *infile)
{
  // New loading mechanism for Gadget files
  hid_t file_id, group_id, dataset_id, dataspace_id;
  herr_t status;
  hsize_t dims[2];
  int i, j, k, sorted_groups, auxint, counter, maxsend, maxrecv;
  struct part dp;
  char buff[100];
  long int bloksize;
  FILE *aux_pf = NULL;
  FILE *fp_inp = NULL;

  // Open the input HDF5 file
  // if((file_id = H5Fopen(infile, H5F_ACC_RDONLY, H5P_DEFAULT)) < 0) {
  //   printf("read_gadget cannot open %s\n", infile);
  //   exit(FAILURE);
  // }

  // Unclustered particles and clusters
  (*unclust) = Npart_unclustered; // set unclustered particles
  (*sacum) = NCLUSTERS; // set the number of clusters

  if (task == 0) {
    domain_info[task].istart = 0;
    domain_info[task].iend = Npart_Total;

    domain_info[task].Nparts_per_node = Npart_Total - (*unclust);
    printf("There are %d clustered particles\n", domain_info[task].Nparts_per_node); fflush(stdout);
    printf("There are %d FOF groups\n", (*sacum)); fflush(stdout);
  }
  

  // creo que ya no necesito esto
  if (Particle != NULL) 
    free(Particle);
  
  if (domain_info[task].Nparts_per_node > 0) {
    Particle = (struct part *) malloc((size_t) domain_info[task].Nparts_per_node * sizeof(struct part));
    if (Particle == NULL) {
      printf("No memory available to load particles\n");
      exit(0);
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
  
  return 0;  
  
}

////////////////////////////////// Load Gadget Old
/*
int GalCos_load_Gadget(int *unclust, int *sacum, char *infile)
{
  // loading mechanism for Gadget files
  int dummi,i,j,counter,auxint,sorted_groups,maxsend,maxrecv,k;
  struct part dp;
  char buff[100];
  long int bloksize;
  FILE *aux_pf=NULL;
  FILE *fp_inp=NULL;
  
  // Open the input file
  if((fp_inp=fopen(infile,"r"))==NULL)    
    {
      printf("read_gadget cannot open %s",infile);
      exit(FAILURE);
    }

  // loading in task 0 
  if(task == 0)
    {
      
      sprintf(buff,"%s%s",infile,".grp"); // auxiliary file with group information ESTO YA NO EXISTE EN EL NUEVO METODO
      aux_pf = fopen(buff,"r"); // open the auxiliary file
      fscanf(aux_pf,"%d",&auxint); // read the first integer (total number of particles)
      
      if(auxint != Npart_Total) // check if the number of particles matches
	{
	  printf("Error reading files\n");
	  MPI_Finalize();
	  exit(0);
	}
      
      (*unclust) = 0; // initialize unclust counter
      (*sacum) = 0; // initialize group counter
      
      //
      for(k=0; k<Npart_Total; k++) 
	{
	  fscanf(aux_pf,"%d",&sorted_groups);
	  
	  if(sorted_groups > (*sacum)) // if the current group is larger than the previous maximum
	    (*sacum) = sorted_groups; // update the maximum group number

	  // if the particles do not belong to any group, increment unclust counter (sorted_groups == 0 means no group)
    //esto ya no toca hacerlo porque ya se la cantidad de particulas que no pertenecen a ningun grupo
    if(sorted_groups == 0)
	    (*unclust) = (*unclust) + 1;
	}
      
      fclose(aux_pf);
      
      domain_info[task].istart = 0;
      domain_info[task].iend = Npart_Total;
      
      domain_info[task].Nparts_per_node = Npart_Total - (*unclust);
      printf("There are %d clustered particles\n",domain_info[task].Nparts_per_node); fflush(stdout);
      printf("There are %d FOF groups\n",(*sacum)); fflush(stdout);
      
      if(Particle != NULL) 
	    free(Particle);
      
      if(domain_info[task].Nparts_per_node > 0)
	{
	  
	  Particle = (struct part *) malloc((size_t) domain_info[task].Nparts_per_node*sizeof(struct part));
	  if(Particle == NULL){
	    printf("No memory available to load particles \n");
	    exit(0);
	  }
	  
	}
      
  //NO NECESITO NADA DE ESTO, ESTO SE USA PARA CALCULAR LA


  //////////////////// Positions
      
      sprintf(buff,"%s%s",infile,".grp");
      aux_pf = fopen(buff,"r");
      fscanf(aux_pf,"%d",&auxint);
      
      bloksize = sizeof(int)*3 + sizeof(struct gadget_head);
      fseek(fp_inp,bloksize,SEEK_SET);
      
      i=0;
      for(k=0; k<Npart_Total; k++) 
	{
	  
	  fread(&dp.pos[0],sizeof(float),3,fp_inp);
	  fscanf(aux_pf,"%d",&sorted_groups);
	  
	  if(sorted_groups > 0)
	    {
	      Particle[i].Cluster_ID = sorted_groups-1;
	      
	      for(j=0; j<3; j++)
		Particle[i].pos[j] = dp.pos[j];
	      //Particle[i].pos[j] = COSMIC_TIME*dp.pos[j];
	      
	      Particle[i].mass = PARTMASS;
	      Particle[i].id = i;
	      Particle[i].Oid = k;
	      
	      i++;
	    }
	  
	}
      
      printf("i=%d %d\n",i,domain_info[task].Nparts_per_node); fflush(stdout);
      
      fclose(aux_pf);
      fread(&dummi,sizeof(dummi),1,fp_inp);
      
      //////////////////// velocities
      
      fread(&dummi,sizeof(dummi),1,fp_inp);
      
      i=0;
      for(k=0; k<Npart_Total; k++) 
	{
	  fread(&dp.vel[0],sizeof(float),3,fp_inp);
	  
	  if((Particle[i].Oid == k) && (Particle[i].Cluster_ID >= 0))
	    {
	      for(j=0; j<3; j++) 
		Particle[i].vel[j] = dp.vel[j]; 
	      //Particle[i].vel[j] = sqrt(COSMIC_TIME)*dp.vel[j]; 
	      
	      i++;
	    }
	  
	}
      
      printf("i=%d %d\n",i,domain_info[task].Nparts_per_node); fflush(stdout);
      fread(&dummi,sizeof(dummi),1,fp_inp);
      
    }
    
  
  if(task != 0)
    {
	
      if(Particle != NULL) 
	free(Particle);
      
      if(domain_info[task].Nparts_per_node != 0)
	{
	  
	  Particle = (struct part *) malloc((size_t) domain_info[task].Nparts_per_node*sizeof(struct part));
	  if(Particle == NULL){
	    printf("No memory available to load particles \n");
	    exit(0);
	  }
	  
	}
      
      //////////////////// Positions
      
      sprintf(buff,"%s%s",infile,".grp");
      aux_pf = fopen(buff,"r");
      fscanf(aux_pf,"%d",&auxint);
      
      // Positioning the pointers on the rigth position 
      
      for(k=0; k<domain_info[task].istart; k++)
	fscanf(aux_pf,"%d",&auxint);
      
      bloksize = sizeof(int)*3 + sizeof(struct gadget_head);
      fseek(fp_inp,bloksize+sizeof(float)*3*domain_info[task].istart,SEEK_SET);
      
      i=0;
      for(k=domain_info[task].istart; k<domain_info[task].iend; k++) 
      {
	  
	  fread(&dp.pos[0],sizeof(float),3,fp_inp);
	  fscanf(aux_pf,"%d",&sorted_groups);
	  
	  if(sorted_groups > 0)
	    Particle[i].Cluster_ID = sorted_groups-1;
	  else
	    Particle[i].Cluster_ID = EMPTY_FLAG;
	  
	  for(j=0; j<3; j++)
	    Particle[i].pos[j] = dp.pos[j];
	  //Particle[i].pos[j] = COSMIC_TIME*dp.pos[j];
	  
	  Particle[i].mass = PARTMASS;
	  Particle[i].id = i;
	  Particle[i].Oid = k;
	  
	  i++;
      }
      
      fclose(aux_pf);
      fread(&dummi,sizeof(dummi),1,fp_inp);
      
      //////////////////// velocities
      
      bloksize = sizeof(int)*5 + sizeof(struct gadget_head) + sizeof(float)*3*Npart_Total;
      fseek(fp_inp,bloksize + sizeof(float)*3*domain_info[task].istart,SEEK_SET);
      
      i=0;
      for(k=domain_info[task].istart; k<domain_info[task].iend; k++) 
	{
	  fread(&dp.vel[0],sizeof(float),3,fp_inp);
	  
	  for(j=0; j<3; j++)
	    Particle[i].vel[j] = dp.vel[j];
	  //Particle[i].vel[j] = sqrt(COSMIC_TIME)*dp.vel[j];
	  
	  i++;
	}
      
      fread(&dummi,sizeof(dummi),1,fp_inp);
      
      (*unclust) = 0;
      (*sacum) = 0;
      
    }
  
  fclose(fp_inp);
  MPI_Barrier(MPI_COMM_WORLD);
  
  return 0;  
  
}


int reload_task0(char *infile)
{
  
  int dummi,i,j,counter,auxint,sorted_groups,maxsend,maxrecv,k;
  struct part dp;
  char buff[100];
  FILE *aux_pf=NULL;
  FILE *fp_inp=NULL;
  
  long int bloksize;
  
  if((fp_inp=fopen(infile,"r"))==NULL)    
    {
      printf("read_gadget cannot open %s",infile);
      exit(FAILURE);
    }
  

  if(domain_info[task].Nparts_per_node != 0)
    {
      
      Particle = (struct part *) malloc((size_t) domain_info[task].Nparts_per_node*sizeof(struct part));
      if(Particle == NULL){
	printf("No memory available to load particles \n");
	exit(0);
      }
      
    }
  
  //////////////////// Positions
  
  
  sprintf(buff,"%s%s",infile,".grp");
  aux_pf = fopen(buff,"r");
  fscanf(aux_pf,"%d",&auxint);
  
  // Positioning the pointers on the rigth position 
  
  //for(k=0; k<istart; k++) 
  //fscanf(aux_pf,"%d",&auxint);
  
  bloksize = sizeof(int)*3 + sizeof(struct gadget_head);
  fseek(fp_inp,bloksize+sizeof(float)*3*domain_info[task].istart,SEEK_SET);
  
  i=0;
  for(k=domain_info[task].istart; k<domain_info[task].iend; k++) 
    {
      fread(&dp.pos[0],sizeof(float),3,fp_inp);
      fscanf(aux_pf,"%d",&sorted_groups);
      
      if(sorted_groups > 0)
	Particle[i].Cluster_ID = sorted_groups-1;
      else
	Particle[i].Cluster_ID = EMPTY_FLAG;
      
      for(j=0; j<3; j++) 
	Particle[i].pos[j] = dp.pos[j];
      //Particle[i].pos[j] = COSMIC_TIME*dp.pos[j];
      
      Particle[i].mass = PARTMASS;
      Particle[i].id = i;
      Particle[i].Oid = k;
	  
      i++;
    }
  
  fclose(aux_pf);
  fread(&dummi,sizeof(dummi),1,fp_inp);
  
  //////////////////// velocities
  
  bloksize = sizeof(int)*5 + sizeof(struct gadget_head) + sizeof(float)*3*Npart_Total;
  fseek(fp_inp,bloksize+sizeof(float)*3*domain_info[task].istart,SEEK_SET);
  fread(&dummi,sizeof(dummi),1,fp_inp);
  
  i=0;
  for(k=domain_info[task].istart; k<domain_info[task].iend; k++) 
    {
      
      fread(&dp.vel[0],sizeof(float),3,fp_inp);
      
      for(j=0; j<3; j++)
	Particle[i].vel[j] = dp.vel[j];
      //Particle[i].vel[j] = sqrt(COSMIC_TIME)*dp.vel[j]; 
      
      i++;
    }
  
  fread(&dummi,sizeof(dummi),1,fp_inp);
  
  fclose(fp_inp);
  return 0;  
  
}
*/
