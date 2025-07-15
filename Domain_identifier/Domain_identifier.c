/* 
   This program takes a snapshot of a simulations and the list of FOF
   halos. With that, the code identifies the domain PARTICLES of every
   FOF halo in that snapshot.
*/

#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<time.h>
#include<mpi.h>

#include "GalCos_variables.h"
#include "GalCos_variables.c"

#define MAX_FILENAME 100

int task,Number_Processes;

// structure to store the domain information of every processor
// Nparts_per_node: number of particles in the node
// istart: first particle in the node
// iend: last particle in the node
struct domain
{
    int Nparts_per_node, istart, iend;
} *domain_info;

void usage_main(void)
{
    
    if(task == 0)
    {
	printf("Domain_identifier = Identifies the domain of FOF halos.\n"); fflush(stdout);
	printf("Juan Carlos Munoz C.\n"); fflush(stdout);
	printf("Usage:./Domain_identifier <Gadget_snapshoth_file> <Gadget_catalog_file> <param_file>\n"); fflush(stdout);
    }
    
    MPI_Finalize();
    exit(0);
}

#include "GalCos_routines.c"
#include "GalCos_load_Gadget.c"

#include "GalCos_units.c"
#include "GalCos_cola.c"
#include "GalCos_clusterBound.c"

#include "GalCos_Halo_prop.c"
#include "GalCos_domain.c"

#include "GalCos_PartitionParticles.c"

// structure to store the information of every halo
struct data
{
    float mass;
    float pos[3];
    float Rvir;
    int Nmembers;
    int NDomain_particles;
    int IDcluster;
};

// Function to write a rescue file with the information of the halos
// and the particles of the halos in a binary file for backup purposes
// storing the halo mass, position, radius, number of members,
// number of domain particles, and the ID of the cluster
int write_rescue(char *infile)
{
    int i;
    FILE *pf;
    char buf[80];
    struct data aux;
    
    sprintf(buf,"%s%s",infile,".rescue");
    pf=fopen(buf,"wb");
    
    fwrite(&NCLUSTERS,sizeof(int),1,pf);
    
    for(i=0; i<NCLUSTERS; i++)
    {
	aux.mass = Halos[i].mass;
	aux.pos[0] = Halos[i].pos[0];
	aux.pos[1] = Halos[i].pos[1];
	aux.pos[2] = Halos[i].pos[2];
	aux.Nmembers = Halos[i].Nmembers; // number of particles in the halo	
	aux.NDomain_particles = Halos[i].NDomain_particles; // number of particles in the domain of halo
	aux.IDcluster = Halos[i].IDcluster; 
	aux.Rvir = Halos[i].Rvir;
	
	fwrite(&aux,sizeof(struct data),1,pf);
	fwrite(&(Halos[i].Halo_particles[0]),sizeof(int),Halos[i].Nmembers,pf);
	fwrite(&(Halos[i].Domain_particles[0]),sizeof(int),Halos[i].NDomain_particles,pf);
    }
    
    fclose(pf);
    
    return 0;
}


int main(int argc, char *argv[])
{
    
    int i,Tini,counter,iadvance,NGROUPS=0,UNCLUSTERED,Number_Processes;
    int info[3],ihalo,Oistart,Oiend,ONparts_per_node,ipart,Nparts_in_node;
    char *param_file=NULL;
    char buf[80],bufhalo[100];
    FILE *aux_pf=NULL;
    FILE *pf=NULL;
    
    int COLLECT_TAG=1,Nparticles_in_Node,k,j,*Info_domains=NULL,sendbuff;
    int NDomain_parts_per_proc,Parts_in_Domain_per_node,l,m,Old_NDomain_particles;
    MPI_Status status;
  
    struct halo *auxHalos=NULL;
    
	// Initializing MPI
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&task);
    MPI_Comm_size(MPI_COMM_WORLD,&Number_Processes);

    if(task == 0)
	printf(" * Running on %d processors\n",Number_Processes); fflush(stdout);
    
    if(argc != 4)
	usage_main();
    
    data_prefix = argv[1];       // Ej: "data_"
    catalog_prefix = argv[2];    // Ej: "catalog_"

	// Reading input parameters
    char infile[MAX_FILENAME], cat_file[MAX_FILENAME];

	//¿¿¿¿ESTO SOLO LO DEBERIA HACER EL TASK 0?????????

	snprintf(infile, MAX_FILENAME, "%s.%d.hdf5", data_prefix, 0); // Gadget snapshot file
	snprintf(cat_file, MAX_FILENAME, "%s.%d.hdf5", catalog_prefix, 0); // FOF catalog file

	param_file = argv[3];


    Tini = tiempoc();

	//¿¿¿¿ESTO SOLO LO DEBERIA HACER EL TASK 0?????????

    GalCos_units(param_file); // Load the units and parameters from the parameter file
    GalCos_preload_Gadget(infile); // Preload the Gadget snapshot file to get other information
	GalCos_preload_Catalog(cat_file); // Preload the FOF catalog file to get the halo information
    

    /////////////////////////////////////////////////////////////////////////////
    /*                          DOMAIN DECOMPOSITION                           */
    /////////////////////////////////////////////////////////////////////////////
    
    domain_info = (struct domain *) malloc((size_t) Number_Processes*sizeof(struct domain));
        
    if(task == 0) // If task 0, initialize the domain information
    {
	Nparts_in_node = (int) (Npart_Total/Number_Processes);
	
	printf("Le toca %d particulas a cada nodo\n",Nparts_in_node); fflush(stdout);
	printf("sobran %d particles\n",(Npart_Total%Number_Processes)); fflush(stdout);
	
	for(i=0; i<Number_Processes; i++)
	{
	    domain_info[i].istart = i*Nparts_in_node;
	    domain_info[i].iend = (i+1)*Nparts_in_node;
	    domain_info[i].Nparts_per_node = domain_info[i].iend - domain_info[i].istart;
	}
	
	i = Number_Processes-1;
	domain_info[i].iend = Npart_Total;
	domain_info[i].Nparts_per_node = domain_info[i].iend - domain_info[i].istart;
	
	i = 0;
	Oistart = domain_info[i].istart;
	Oiend = domain_info[i].iend;
	ONparts_per_node = domain_info[i].Nparts_per_node;
	
	for(i=0; i<Number_Processes; i++)
	    printf("task %d start at %d and ends at %d, nparts=%d\n",i,domain_info[i].istart,domain_info[i].iend,
		   domain_info[i].Nparts_per_node); fflush(stdout);
	
    }
    
    MPI_Bcast(&domain_info[0],Number_Processes*sizeof(struct domain),MPI_BYTE,0,MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    


    /////////////////////////////////////////////////////////////////////////////////////
    //                 STARTING WITH FOF HALO IDENTIFICATION                           // 
    /////////////////////////////////////////////////////////////////////////////////////
    

    if(task == 0) 
	printf(" >> LOADING PARTICLE DATA FROM GADGET FILE...\n"); fflush(stdout);
	
	    
    GalCos_load_Gadget(&UNCLUSTERED,&NGROUPS,infile);
    
    MPI_Bcast(&NGROUPS,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&UNCLUSTERED,1,MPI_INT,0,MPI_COMM_WORLD);
    
    if(task == 0)	
    {
	printf(" %g Memory in particles = %g %d\n",sizeof(struct part)/(1024*1024.0),
	       Npart_Total*sizeof(struct part)/(1024*1024.0),Npart_Total); fflush(stdout);
	printf(" %g Memory in halos = %g\n",sizeof(struct halo)/(1024*1024.0),
	       60000*sizeof(struct halo)/(1024*1024.0)); fflush(stdout);
	printf(" THERE ARE %d UNCLUSTERED PARTICLES (%f percent)\n",UNCLUSTERED,
	       100.0*(1.0*UNCLUSTERED)/(Npart_Total*1.0)); fflush(stdout);
	//printf(" There are %d fof halos\n",NGROUPS); fflush(stdout);
    }

	
    
    /////////////////////////////////////////////////////////////////////////////////////
    //                 FILLING HALO STRUCTURES WITH PARTICLE ID's                      //
    /////////////////////////////////////////////////////////////////////////////////////
    
	if(task == 0)
	{

    Halos = (struct halo *) malloc((size_t) NGROUPS*sizeof(struct halo));
    if(Halos == NULL)
    {
	printf("There are no memory left to allocate Halos\n");
	MPI_Finalize();
	exit(0);
    }
    
	

	int Halo_counter = 0; //Contador de halos
	int Nmembers1 = 0; //Contador de miembros de halos, para asignar los ids de las particulas del halo

	for(int fileid = 0; fileid<NumFiles; fileid++)
	{
		
		printf("\nReading catalog %d.....\n", fileid); fflush(stdout);

		char cat_file[100];    

		hid_t file_id, group_id, dataset_id, dataspace_id, header_group, attr_id;
		herr_t status;
		hsize_t dims[1];
		
		sprintf(cat_file, "%s.%d.hdf5", argv[2], fileid);

		file_id = H5Fopen(cat_file, H5F_ACC_RDONLY, H5P_DEFAULT);
		if(file_id < 0) {
			printf("Cannot open %s\n", cat_file);
			MPI_Finalize();
			exit(0);
		}

		//Read Header to get the number of groups in this file
		int NGROUPS_THISFILE;

		if((header_group = H5Gopen(file_id, "/Header", H5P_DEFAULT)) < 0) {
			printf("Error opening Header group in %s\n", cat_file);
			exit(FAILURE);
		}

		if((attr_id = H5Aopen(header_group, "Ngroups_ThisFile", H5P_DEFAULT)) < 0) {
			printf("Error opening NgroupsThisFile attribute in %s\n", cat_file);
			exit(FAILURE);
		}
		
		status = H5Aread(attr_id, H5T_NATIVE_INT, &NGROUPS_THISFILE);
		
		H5Aclose(attr_id);
		H5Gclose(header_group);

		printf("Number of groups in this file: %d\n", NGROUPS_THISFILE); fflush(stdout);
	
		group_id = H5Gopen(file_id, "Group", H5P_DEFAULT);
		if(group_id < 0) {
			printf("Cannot open Group in %s\n", cat_file);
			MPI_Finalize();
			exit(0);
		}

		//Open the GroupLen dataset to get the number of members in each group
		dataset_id = H5Dopen(group_id, "GroupLen", H5P_DEFAULT);
		if(dataset_id < 0) {
			printf("Cannot open GroupLen in %s\n", cat_file);
			MPI_Finalize();
			exit(0);
		}

		int *group_members = (int *) malloc(NGROUPS_THISFILE * sizeof(int));

		if(group_members == NULL) {
			printf("No memory available for group members\n");
			MPI_Finalize();
			exit(0);
		}
		status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, group_members);
		if(status < 0) {
			printf("Cannot read GroupLen in %s\n", cat_file);
			MPI_Finalize();
			exit(0);
		}
		H5Dclose(dataset_id);

		//Open the GroupMass dataset to get the mass of each group
		dataset_id = H5Dopen(group_id, "GroupMass", H5P_DEFAULT);
		if(dataset_id < 0) {
			printf("Cannot open GroupMass in %s\n", cat_file);
			MPI_Finalize();
			exit(0);
		}
		double *group_mass = (double *) malloc(NGROUPS_THISFILE * sizeof(double));
		if(group_mass == NULL) {
			printf("No memory available for group mass\n");
			MPI_Finalize();
			exit(0);
		}
		status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, group_mass);
		if(status < 0) {
			printf("Cannot read GroupMass in %s\n", cat_file);
			MPI_Finalize();
			exit(0);
		}
		H5Dclose(dataset_id);

		//Open the GroupPos dataset to get the positions of each group
		dataset_id = H5Dopen(group_id, "GroupPos", H5P_DEFAULT);
		if(dataset_id < 0) {
			printf("Cannot open GroupPos in %s\n", cat_file);
			MPI_Finalize();
			exit(0);
		}
		float (*group_pos)[3] = (float (*)[3]) malloc(NGROUPS_THISFILE * 3 * sizeof(float));
		if(group_pos == NULL) {
			printf("No memory available for group positions\n");
			MPI_Finalize();
			exit(0);
		}
		status = H5Dread(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, group_pos);
		if(status < 0) {
			printf("Cannot read GroupPos in %s\n", cat_file);
			MPI_Finalize();
			exit(0);
		}
		H5Dclose(dataset_id);

		//Open the Group_R_TopHat200 dataset to get the radius of each group
		dataset_id = H5Dopen(group_id, "Group_R_TopHat200", H5P_DEFAULT);
		if(dataset_id < 0) {
			printf("Cannot open Group_R_TopHat200 in %s\n", cat_file);
			MPI_Finalize();
			exit(0);
		}
		float *group_radius = (float *) malloc(NGROUPS_THISFILE * sizeof(float));
		if(group_radius == NULL) {
			printf("No memory available for group radius\n");
			MPI_Finalize();
			exit(0);
		}
		status = H5Dread(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, group_radius);
		if(status < 0) {
			printf("Cannot read Group_R_TopHat200 in %s\n", cat_file);
			MPI_Finalize();
			exit(0);
		}
		H5Dclose(dataset_id);

		//Open the Group_M_TopHat200 dataset to get the mass of each group
		dataset_id = H5Dopen(group_id, "Group_M_TopHat200", H5P_DEFAULT);
		if(dataset_id < 0) {
			printf("Cannot open Group_M_TopHat200 in %s\n", cat_file);
			MPI_Finalize();
			exit(0);
		}
		double *group_mass_top_hat = (double *) malloc(NGROUPS_THISFILE * sizeof(double));
		if(group_mass_top_hat == NULL) {
			printf("No memory available for group mass top hat\n");
			MPI_Finalize();
			exit(0);
		}
		status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, group_mass_top_hat);
		if(status < 0) {
			printf("Cannot read Group_M_TopHat200 in %s\n", cat_file);
			MPI_Finalize();
			exit(0);
		}
		H5Dclose(dataset_id);

		H5Gclose(group_id);
		H5Fclose(file_id);
		

		for (int i = 0; i < NGROUPS_THISFILE; i++) {
			Halos[Halo_counter].Nmembers = group_members[i];
			Halos[Halo_counter].mass = group_mass[i];
			Halos[i].NDomain_particles = 0;	
			Halos[Halo_counter].pos[0] = group_pos[i][0];
			Halos[Halo_counter].pos[1] = group_pos[i][1];
			Halos[Halo_counter].pos[2] = group_pos[i][2];
			Halos[Halo_counter].IDcluster = Halo_counter; // Assign the cluster ID
			Halos[Halo_counter].Rvir = group_radius[i]; // Assign the radius of the group
			Halos[Halo_counter].Mvir = group_mass_top_hat[i]; // Assign the mass of the group
			
			//Esto es para las particulas del halo, las cuales no necesitamos
			// Halos[Halo_counter].NDomain_particles = 0; //
			// Halos[Halo_counter].Halo_particles = realloc(Halos[Halo_counter].Halo_particles, (size_t) group_members[i] * sizeof(int));

			// for (int j = 0; j < group_members[i]; j++) {
			// 	Halos[Halo_counter].Halo_particles[j] = Nmembers1;
			// 	Nmembers1++;
			// }
			
			Halo_counter++;

			// //imprimir el id de la ultima particula del ultimo halo
			// if (Halo_counter == NGROUPS) {
			// 	printf("\nLast particle ID of halo %d: %d\n", Halo_counter - 1, Halos[Halo_counter - 1].Halo_particles[Halos[Halo_counter - 1].Nmembers - 1]);
			// }

			if (i==0) {
				printf("Group %d has %d members\n", Halo_counter-1, group_members[i]); fflush(stdout);
				printf("Group %d has mass %g\n", Halo_counter-1, group_mass[i]); fflush(stdout);
				printf("Group %d position: (%g, %g, %g)\n", Halo_counter-1, group_pos[i][0], group_pos[i][1], group_pos[i][2]); fflush(stdout);
				printf("Group %d radius (TopHat200): %g kpc\n", Halo_counter-1, group_radius[i]); fflush(stdout);
				printf("Group %d mass (TopHat200): %ge10 Msun\n", Halo_counter-1, group_mass_top_hat[i]); fflush(stdout);
			}
			
		}
		free(group_members);
		free(group_mass);
		free(group_pos);
		
	}
	}
	
	
	//Esto ya no, está en el catalogo de halos, puedo obtener la info de ahí
    // for(i=0; i<NGROUPS; i++)
    // {
	// Halos[i].Nmembers = 0;
	// Halos[i].NDomain_particles = 0;
	// Halos[i].Halo_particles = NULL;
	// Halos[i].Domain_particles = NULL;
    // }
    
    // if(task == 0)
    // {
	
	// iadvance = 1000000;
	
	// for(i=0; i<domain_info[task].Nparts_per_node; i++)
	// {
	    
	//     if((i%iadvance) == 0)
	// 	printf(" *Evaluated %d of %d at %f min\n",i,Npart_Total,(tiempoc()-Tini)/60.0); fflush(stdout);
	    
	//     if(Particle[i].Cluster_ID != EMPTY_FLAG)
	//     {
	// 	ihalo = Particle[i].Cluster_ID;
	// 	Halos[ihalo].Nmembers++;
	// 	Halos[ihalo].Halo_particles = realloc(Halos[ihalo].Halo_particles, (size_t) Halos[ihalo].Nmembers*sizeof(int));
	// 	Halos[ihalo].Halo_particles[Halos[ihalo].Nmembers-1]= Particle[i].id;
	//     }
	    
	// }

	if(task == 0)
    {
	
	aux_pf=fopen("Mass_function.dat","w");
	for(i=0; i<NGROUPS; i++)
	    fprintf(aux_pf,"%g\n",Halos[i].mass);

	printf("\nMass function written\n"); fflush(stdout);
	fclose(aux_pf);

	printf("There are %d halos\n",NGROUPS); fflush(stdout);
	NCLUSTERS = NGROUPS;
	}
	
	if (task == 0)
	{ 
	//old
	// for(i=0; i<NGROUPS; i++)
	// {
	//     Halos[i].ID_CenterHalo = Halos[i].Halo_particles[0];
	//     Halos[i].IDcluster = i;
	//     Halos[i].Domain_particles=NULL;
	//     Halos[i].NDomain_particles=0;
	// }
	
	/////////////////////////////////////////////////////////////////////////////////////
	//                     CENTERING AND PERIODIC BOUNDARY CONDITIONS                  //
	/////////////////////////////////////////////////////////////////////////////////////
	
	iadvance = 10000;
	auxHalos = (struct halo *) malloc((size_t) NCLUSTERS*sizeof(struct halo));
	aux_pf=fopen("Halo_Catalog.dat","w");
	
	for(i=0; i<NCLUSTERS; i++)
	{
		//periodic_boundary_corrections(i);
		//Halos[i].ID_CenterHalo = GalCos_clusterBound(Halos[i].Halo_particles,Halos[i].Nmembers,i);
		//GalCos_Halo_prop(&(Halos[i]));

		fprintf(aux_pf,"%g %g %g %g %g\n",Halos[i].pos[0],Halos[i].pos[1],Halos[i].pos[2],Halos[i].Mvir,Halos[i].Rvir);
	    fflush(stdout);
		//auxHalos[i].Halo_particles = (int *) malloc((size_t) Halos[i].Nmembers*sizeof(int));
		//for(j=0; j<Halos[i].Nmembers; j++)
			//auxHalos[i].Halo_particles[j] = Halos[i].Halo_particles[j];//define auxHalos without Particle
			//auxHalos[i].Halo_particles[j] = Particle[Halos[i].Halo_particles[j]].Oid; //old way, using Particle
	    
		// free(Halos[i].Halo_particles);
	    // Halos[i].Halo_particles=NULL;
	    
	    
	}

	fclose(aux_pf);
	free(Particle);
	
	domain_info[task].istart = Oistart;
	domain_info[task].iend = Oiend;
	domain_info[task].Nparts_per_node = ONparts_per_node;
	
      //reload_task0(infile);
      
      printf("Done halo properties\n"); fflush(stdout);
      
    }

	if (task == 0)
	{
		GalCos_PartitionParticles(4400000, 5400000);
	}
	/////////////////////////////
	MPI_Finalize();
	exit(0);
	/////////////////////////////
	
  
  MPI_Bcast(&NCLUSTERS,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&Halos[0],NCLUSTERS*sizeof(struct halo),MPI_BYTE,0,MPI_COMM_WORLD);
  
  MPI_Barrier(MPI_COMM_WORLD);
  
  /////////////////////////////////////////////////////////////////////////////////
  /*                       Evaluating the domain of every halo                   */
  /////////////////////////////////////////////////////////////////////////////////
  
  counter=0;
  for(i=0; i<domain_info[task].Nparts_per_node; i++)
  {
      
      if((counter%1000000) == 0)
	  printf(" (%d) *Computed domains for %d at %f min\n",task,counter,
		 (tiempoc()-Tini)/60.0); fflush(stdout);
      
      if(Particle[i].Cluster_ID == EMPTY_FLAG)
	  GalCos_domain(i);
      
      counter++;
  }
  
  MPI_Barrier(MPI_COMM_WORLD);
  
  /////////////////////////////////////////////////////////////////////////////////
  /*      Sending Particle information from the other machines to the 0 task     */
  /////////////////////////////////////////////////////////////////////////////////
  
  /*
    sprintf(buf,"%s%s%d",infile,".parts.rescue.",task);
    pf=fopen(buf,"w");
    for(i=0; i<Nparts_per_node; i++)
    fprintf(pf,"%d %g %g %g\n",Particle[i].Cluster_ID,Particle[i].pos[0],Particle[i].pos[1],
    Particle[i].pos[2]);
    fclose(pf);
  */
  
  for(l=0; l<Number_Processes; l++)
    {
      
      if(task == l)
	{
	  sprintf(buf,"%s%s",infile,".parts.rescue");
	  
	  if(task == 0)
	    pf=fopen(buf,"w");
	  else
	    pf=fopen(buf,"a");
	  
	  for(i=0; i<domain_info[task].Nparts_per_node; i++)
	      fprintf(pf,"%d\n",Particle[i].Cluster_ID);
	  
	  fclose(pf);
	}
      
      MPI_Barrier(MPI_COMM_WORLD);
      
    }
  
  MPI_Barrier(MPI_COMM_WORLD);
  
  
  for(i=1; i<Number_Processes; i++)
    {
      
      if(task == i)
	{
	  
	  /* Sending information of Halos */
	  
	  Parts_in_Domain_per_node=0;
	  for(j=0; j<NCLUSTERS; j++)
	    Parts_in_Domain_per_node += Halos[j].NDomain_particles; 
	  
	  Info_domains = (int *) malloc((size_t) (Parts_in_Domain_per_node + 2*NCLUSTERS)*sizeof(int));
	  if(Info_domains == NULL)
	    {
	      printf("No memory available for allocation\n");
	      exit(0);
	    }
	  
	  k=0;
	  for(j=0; j<NCLUSTERS; j++)
	    {
	      
	      Info_domains[k] = Halos[j].IDcluster;
	      k++;
	      Info_domains[k] = Halos[j].NDomain_particles;
	      k++;
	      
	      for(l=0; l<Halos[j].NDomain_particles; l++)
		{
		  Info_domains[k] = Halos[j].Domain_particles[l];
		  k++;
		}
	      
	    }
	  
	  sendbuff = Parts_in_Domain_per_node + 2*NCLUSTERS;
	  MPI_Ssend(&sendbuff,1,MPI_INT,0,COLLECT_TAG,MPI_COMM_WORLD);
	  MPI_Ssend(&Info_domains[0],sendbuff,MPI_INT,0,COLLECT_TAG,MPI_COMM_WORLD);
	  
	  free(Info_domains);
	  
	}
      
    }
  
  /* Now I have to collect the data of particles from the other tasks */
  
  if(task == 0)
    {
      
      for(i=1; i<Number_Processes; i++)
	{
	  
	  /* Receiving data of halos */
	  
	  MPI_Recv(&sendbuff,1,MPI_INT,i,COLLECT_TAG,MPI_COMM_WORLD,&status);
	  
	  printf("Receiving info of halos from %d %d\n",i,sendbuff); fflush(stdout);
	  
	  Info_domains = (int *) malloc((size_t) sendbuff*sizeof(int));
	  if(Info_domains == NULL)
	    {
	      printf("No memory available for allocation\n");
	      exit(0);
	    }
	  
	  MPI_Recv(&Info_domains[0],sendbuff,MPI_INT,i,COLLECT_TAG,MPI_COMM_WORLD,&status);
	  printf("recibidos\n"); fflush(stdout);
	  
	  k=0;
	  for(j=0; j<NCLUSTERS; j++)
	    {
	      ihalo = Info_domains[k];
	      k++;
	      
	      Old_NDomain_particles = Halos[ihalo].NDomain_particles;
	      NDomain_parts_per_proc = Info_domains[k];
	      
	      Halos[ihalo].NDomain_particles += Info_domains[k];
	      k++;
	      
	      Halos[ihalo].Domain_particles = realloc(Halos[ihalo].Domain_particles,(size_t) Halos[ihalo].NDomain_particles*sizeof(int));
	      
	      m = Old_NDomain_particles;
	      for(l=0; l<NDomain_parts_per_proc; l++)
		{
		  Halos[ihalo].Domain_particles[m] = Info_domains[k];
		  k++;
		  m++;
		}
	    }
	  
	  free(Info_domains);
	  
	}
      
    }
  
  MPI_Barrier(MPI_COMM_WORLD);
  
  if(task == 0) 
    {
      
      for(i=0; i<NCLUSTERS; i++)
	{
	  Halos[i].Halo_particles = (int *) malloc((size_t) Halos[i].Nmembers*sizeof(int));
	  
	  for(j=0; j<Halos[i].Nmembers; j++)
	    Halos[i].Halo_particles[j] = auxHalos[i].Halo_particles[j];
	}

      printf("Writing info files\n");
      write_rescue(infile);
      
      printf("all is done!\n");
      fflush(stdout);
    }
  
  free(Particle);
  free(Halos);
  
  MPI_Finalize();
  
  return 0;
  
}
