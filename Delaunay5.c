/*
  Densfield is a code to compute the radial density distrubution of halos and
  its domains.

  The inputs of the code are the simulation snapshot, the file parameters (and
  imoplicit the data file with the info of the FOF groups).

  The outputs are a set of files density_#bin.dat with the radial density
  distribution of every input poplation in the simulation.

  For reference about terminology, see Wang et al. 2008.

*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>

int Number_Processes;
char arg1[30];

#include "GalCos_variables.h"
#include "GalCos_variables.c"

void usage_main(void)
{
  printf("Densfield = Density field estimator\n");
  printf("Juan Carlos Muñoz C.\n");
  printf("Andrés López Echeverri.\n");
  printf("Usage:./Delaunay5.x <initial_of_snapshots>\n");
  printf("For example, for snapshots called snap_000.d.hdf5 with d different number file just write:./Delaunay5.x snap_000\n");
  exit(SUCCES);
}

struct domain
{
  /* coords of the border of the complete box */
  float zmin; 
  float zmax;
  float ymin;
  float ymax;
  float xmin;
  float xmax;
  
  /* coords of the border of the region enclosing only the domain
     region. Not the border*/
  float zmind;
  float zmaxd;
  float ymind;
  float ymaxd;
  float xmind;
  float xmaxd;
  
  float pos[3]; // center of the box
  int Nparticles; // Number of particles in that box
}*domain_info;

#include "GalCos_routines.c"
#include "GalCos_load_Gadget3.c"
#include "Delaunay_volumes2.c"

int main(int argc, char *argv[])
{
  
  int i,Tini,counter,j,k,l,ipart;
  float cellsize,overlap,zlevel,ylevel,xlevel;
  double Vol;
  int Np1D=4;
  
  char buff[200];
  char infile[100];
  FILE *pf = NULL;
  FILE *pf2 = NULL;

  Number_Processes = Np1D*Np1D*Np1D;
  
  if (argc < 2)
    usage_main();

  // Asignar argv[1] a la variable global arg1
  strcpy(arg1, argv[1]);
  sprintf(infile, "%s.%d.hdf5", argv[1], 0);
  Tini = tiempoc();

  GalCos_preload_Gadget(infile);
  
  /////////////////////////////////////////////////////////////
  /*                     Data Decomposition                  */
  /////////////////////////////////////////////////////////////
  domain_info = (struct domain *) malloc((size_t) Number_Processes*sizeof(struct domain));
  
  cellsize = BoxSize/(1.0*Np1D);
  overlap = cellsize/10.0;
  
  printf(" \nCellsize = %g\n",cellsize);
  printf(" overlap = %g\n",overlap);
  printf(" Np1D = %d\n\n",Np1D);
  
  zlevel = 0.0;
  
  counter=0;
  for(i=0; i<Np1D; i++)
    {
      
      ylevel = 0.0;
      
      for(j=0; j<Np1D; j++)
	{
	  
	  xlevel = 0.0;
	  
	  for(k=0; k<Np1D; k++)
	    {
	      
	      domain_info[counter].zmin = zlevel - overlap;
	      domain_info[counter].zmax = zlevel + cellsize + overlap;
	      domain_info[counter].pos[2] = (domain_info[counter].zmin + domain_info[counter].zmax)*0.5;
	      
	      domain_info[counter].ymin = ylevel - overlap;
	      domain_info[counter].ymax = ylevel + cellsize + overlap;
	      domain_info[counter].pos[1] = (domain_info[counter].ymin + domain_info[counter].ymax)*0.5;
	      
	      domain_info[counter].xmin = xlevel - overlap;
	      domain_info[counter].xmax = xlevel + cellsize + overlap;
	      domain_info[counter].pos[0] = (domain_info[counter].xmin + domain_info[counter].xmax)*0.5;
	      
	      domain_info[counter].zmind = zlevel;
	      domain_info[counter].zmaxd = zlevel + cellsize;
	      
	      domain_info[counter].ymind = ylevel;
	      domain_info[counter].ymaxd = ylevel + cellsize;
	      
	      domain_info[counter].xmind = xlevel;
	      domain_info[counter].xmaxd = xlevel + cellsize;
	      
	      xlevel += cellsize;
	      counter++;
	    }
	  
	  ylevel += cellsize;
      }
      
      zlevel += cellsize;
      
  }

  //Let's read all the particles and assign a BoxID, ID, FileID to each one of them corresponding to the domain they belong to
  //We will use the domain_info structure to store the BoxID of each particle 
  //The idea is to have each particle well identified (ID) with the domain it belongs to (BoxID) or if the particle enters into
  //another domain, in that case it will appear in more than one Box
  int ID = 0;
  printf("Npart_Total = %d\n", Npart_Total);

  int maxb = 8;
  float Plength,Npos[3],posa[3];

  Rpart = (struct part *) malloc((size_t) maxb*Npart_Total*sizeof(struct part));
  
  Plength = pow(0.5*BoxSize,2);
  printf("\nPlength: %f\n", Plength);

  int N_Rpart_Total = 0;
  
  //abrir el archivo pf2 para meter las posiciones de las particulas
  pf2 = fopen("posiciones.dat","w");  

  for (int fileID=0; fileID<NumFiles; fileID++)
  {  
    
    printf("Reading file %d\n", fileID);
    char infile[100];    

    hid_t file_id, group_id, dataset_id, dataspace_id, header_group, attr_id;
    herr_t status;
    hsize_t dims[2];
    
    float (*positions)[3];  // Matriz para almacenar las posiciones
    sprintf(infile, "%s.%d.hdf5", argv[1], fileID);
  
    file_id = H5Fopen(infile, H5F_ACC_RDONLY, H5P_DEFAULT);
    
    if (file_id < 0) {
      fprintf(stderr, "Error: no se pudo abrir el archivo %s\n", infile);
      continue;
    }
      
    // Open the coordinate dataset
    dataset_id = H5Dopen2(file_id, "/PartType1/Coordinates", H5P_DEFAULT);
    dataspace_id = H5Dget_space(dataset_id);
    H5Sget_simple_extent_dims(dataspace_id, dims, NULL);
    

    // Read the positions of the particles one by one and assign the BoxID depending on the domain_info structure
    positions = malloc(3 * dims[0] * sizeof(float));  
    status = H5Dread(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, positions);
    if (status < 0) {
        fprintf(stderr, "Error al leer datos de HDF5\n");
        exit(1);
    }


  printf(">>>>>>Assigning particles of file %d to domains>>>>>>>>>>>>>\n", fileID);
  printf("Number of particles: %d\n", (int)dims[0]);
  struct part dp;
    for (int i = 0; i < (int)dims[0]; i++) {

      if (i%1000000 == 0) {
        printf("i = %d\n", i);
      }

      dp.pos[0] = positions[i][0];
      dp.pos[1] = positions[i][1];
      dp.pos[2] = positions[i][2];
      
      int IDx = (int)(dp.pos[0] / cellsize);
      int IDy = (int)(dp.pos[1] / cellsize);
      int IDz = (int)(dp.pos[2] / cellsize);
   
      
      for (int dx = -1; dx <= 1; dx++) {
	int ix = (IDx + dx + Np1D)%Np1D;
        for (int dy = -1; dy <= 1; dy++) {
	  int iy = (IDy + dy + Np1D)%Np1D;
          for (int dz = -1; dz <= 1; dz++) {
	    int iz = (IDz + dz + Np1D)%Np1D;
            
            //linear index
            int j = ix + iy*Np1D + iz*Np1D*Np1D;
            
            posa[0] = domain_info[j].pos[0];
            posa[1] = domain_info[j].pos[1];
            posa[2] = domain_info[j].pos[2];

            Npos[0] = dp.pos[0];
            Npos[1] = dp.pos[1];
            Npos[2] = dp.pos[2];

            //Periodic adjustment
            if(pow(dp.pos[0]-posa[0],2) > Plength) {
              if(posa[0] > dp.pos[0])
                Npos[0] = dp.pos[0] + BoxSize;
              else
                Npos[0] = dp.pos[0] - BoxSize;
              }
          
            if(pow(dp.pos[1]-posa[1],2) > Plength) {
              if(posa[1] > dp.pos[1])
                Npos[1] = dp.pos[1] + BoxSize;
              else
                Npos[1] = dp.pos[1] - BoxSize;
              }
          
            if(pow(dp.pos[2]-posa[2],2) > Plength) {
              if(posa[2] > dp.pos[2])
                Npos[2] = dp.pos[2] + BoxSize;
              else
                Npos[2] = dp.pos[2] - BoxSize;
              }
              
            if((Npos[0] >= domain_info[j].xmin) && (Npos[0] <= domain_info[j].xmax)  &&
               (Npos[1] >= domain_info[j].ymin) && (Npos[1] <= domain_info[j].ymax) && 
               (Npos[2] >= domain_info[j].zmin) && (Npos[2] <= domain_info[j].zmax)) {

              Rpart[ID].BoxID = j;
              Rpart[ID].FileID = fileID;
              Rpart[ID].ID = i; //this id the ID of the particle in the file
              Rpart[ID].tag = 0;
	            Rpart[ID].Volume = 0.0;

              if((dp.pos[0] >= domain_info[j].xmind) && (dp.pos[0] <= domain_info[j].xmaxd)){
                if((dp.pos[1] >= domain_info[j].ymind) && (dp.pos[1] <= domain_info[j].ymaxd)){
                  if((dp.pos[2] >= domain_info[j].zmind) && (dp.pos[2] <= domain_info[j].zmaxd)){
                    Rpart[ID].tag = 1;
                    // Store the position in the file
                    fprintf(pf2, "%lf %lf %lf\n", dp.pos[0], dp.pos[1], dp.pos[2]);
                  }
                }
              }
              domain_info[j].Nparticles++;
              N_Rpart_Total++;
              ID = ID + 1;

          }
        }
      }
    }
  }

    // Close the file
  free(positions);
  H5Dclose(dataset_id);
  H5Sclose(dataspace_id);
  H5Fclose(file_id);
  printf(">>>>>>Particles of file %d assigned to domains>>>>>>>>>>>>>\n", fileID);

  }
  fclose(pf2);

  for(l=0; l<Number_Processes; l++)
      {
        printf(" >> LOADING PARTICLE DATA FROM GADGET FILE...in task %d\n",l);
        int Nparts = GalCos_load_Gadget(infile,overlap,l,N_Rpart_Total);
        
        if(Nparts>0)
    {
      
      Delaunay_volumes(l);         
      printf("Writing in task %d\n",l);
      
      for(i=0; i<domain_info[l].Nparticles; i++)
        {
      
          ipart = Particle[i].Oid;
          if(Particle[i].tag == 1) // copying just the info of particles in the data domain
          Rpart[ipart].Volume = Particle[i].Volume;
          
        }
      
      printf(" %d >> DONE %f min\n",l,(tiempoc()-Tini)/60.0); fflush(stdout);
      free(Particle);
      printf("\n");
      
    }
        
      }
      
    printf("Writing...\n");
    sprintf(buff,"%s","volumenes.dat");
    
    pf=fopen(buff,"w");

    for(i=0; i<N_Rpart_Total; i++)
      {
	      if(Rpart[i].tag==0) continue;
        Vol = Rpart[i].Volume;
        fwrite(&Vol, sizeof(double), 1, pf);
      }
    fclose(pf);
    free(Rpart);

    
  return 0;
  
}
