#include "GalCos_variables.h"

#define MAXPATHLEN 256

int GalCos_preload_Gadget(char *infile)
{
  //Leer hdf5
  //hid_t dataset_id, dataspace_id;
  hid_t file_id, header_group, attr_id;
  herr_t status;
  //hsize_t dims[2];

  // Abre el archivo HDF5
  file_id = H5Fopen(infile, H5F_ACC_RDONLY, H5P_DEFAULT);
  header_group = H5Gopen(file_id, "Header", H5P_DEFAULT);

  //Redshift
  attr_id = H5Aopen(header_group, "Redshift", H5P_DEFAULT);
  status = H5Aread(attr_id, H5T_NATIVE_FLOAT, &REDSHIFT);
  H5Aclose(attr_id);
  //Omega0
  attr_id = H5Aopen(header_group, "Omega0", H5P_DEFAULT);
  status = H5Aread(attr_id, H5T_NATIVE_FLOAT, &OMEGA_MATTER);
  H5Aclose(attr_id);
  //OmegaLambda
  attr_id = H5Aopen(header_group, "OmegaLambda", H5P_DEFAULT);
  status = H5Aread(attr_id, H5T_NATIVE_FLOAT, &OMEGALAMBDA);
  H5Aclose(attr_id);
  //HubbleParam
  attr_id = H5Aopen(header_group, "HubbleParam", H5P_DEFAULT);
  status = H5Aread(attr_id, H5T_NATIVE_FLOAT, &HUBBLEPARAM);
  H5Aclose(attr_id);
  //Time???? (Cosmic Time)
  attr_id = H5Aopen(header_group, "Time", H5P_DEFAULT);
  status = H5Aread(attr_id, H5T_NATIVE_FLOAT, &COSMIC_TIME);
  H5Aclose(attr_id);
  //FlagSfr
  attr_id = H5Aopen(header_group, "Flag_Sfr", H5P_DEFAULT);
  status = H5Aread(attr_id, H5T_NATIVE_INT, &FlagSfr);
  H5Aclose(attr_id);
  //FlagFeedback
  attr_id = H5Aopen(header_group, "Flag_Feedback", H5P_DEFAULT);
  status = H5Aread(attr_id, H5T_NATIVE_INT, &FlagFeedback);
  H5Aclose(attr_id);
  //FlagCooling
  attr_id = H5Aopen(header_group, "Flag_Cooling", H5P_DEFAULT);
  status = H5Aread(attr_id, H5T_NATIVE_INT, &FlagCooling);
  H5Aclose(attr_id);
  //NumFiles
  attr_id = H5Aopen(header_group, "NumFilesPerSnapshot", H5P_DEFAULT);
  status = H5Aread(attr_id, H5T_NATIVE_INT, &NumFiles);
  H5Aclose(attr_id);
  //BoxSize
  attr_id = H5Aopen(header_group, "BoxSize", H5P_DEFAULT);
  status = H5Aread(attr_id, H5T_NATIVE_FLOAT, &BoxSize);
  H5Aclose(attr_id);
  //NpartTotal
  int particles[6];
  // attr_id = H5Aopen(header_group, "NumPart_ThisFile", H5P_DEFAULT);
  attr_id = H5Aopen(header_group, "NumPart_Total", H5P_DEFAULT);
  status = H5Aread(attr_id, H5T_NATIVE_INT, &particles[0]);
  H5Aclose(attr_id);
  Npart_Total=particles[1];
  printf("Npart_Total: %d\n", Npart_Total);
  
  printf("\n");
  printf(" * Frame's Time... %g\n",COSMIC_TIME); 
  printf(" * Redshift... %g\n",REDSHIFT);
  printf(" * Flagsfr... %d\n",FlagSfr);
  printf(" * Flagfed... %d\n",FlagFeedback);
  printf(" * Flagcool... %d\n",FlagCooling);
  printf(" * Numfiles... %d\n",NumFiles);
  printf(" * Boxsize... %g\n",BoxSize);
  printf(" * Omega0... %g\n",OMEGA_MATTER);
  printf(" * OmageLambda... %g\n",OMEGALAMBDA);
  printf(" * Hubbleparam... %g\n",HUBBLEPARAM);
      
  H5Gclose(header_group);
  H5Fclose(file_id);

  return 0;
  
}


/// @brief 
/// @param infile 
/// @param overlap 
/// @param l 
/// @return 
int GalCos_load_Gadget(char *infile,float overlap,int l,int N_Rpart_Total)
{
  
  int i,NPH,j,oNPH;
  float Plength,Npos[3],posa[3];
  struct part dp;

  int task = l;
  ///////////////////////////////////Leer HDF5///////////////////
  Plength = pow(0.5*BoxSize,2);
  //printf("Number of particles: %d\n", Npart_Total);
  int Npartbox = domain_info[task].Nparticles;
  printf("Number of particles in box %d: %d\n", task, Npartbox);

  Particle = (struct Vpart *) malloc(Npartbox * sizeof(struct Vpart));
  if (Particle == NULL) {
      fprintf(stderr, "Error al asignar memoria para Particle.\n");
      exit(EXIT_FAILURE);
  }
    
  //////////////////// positions
   
  float coord[3];

  posa[0] = domain_info[task].pos[0];
  posa[1] = domain_info[task].pos[1];
  posa[2] = domain_info[task].pos[2];


  //Let's check all over the Rpart array and get the positions of the particles
  //in the BoxID=task
  int counter = 0;
  for(int FILEid=0; FILEid<NumFiles; FILEid++){
    
    ////leer archivo
    printf("Reading file %d....\n", FILEid);
    char infile[100];

    hid_t file_id, dataset_id, dataspace_id, memspace_id;
    herr_t status;
    hsize_t offset[2], count[2];
    double pos[3];  // Buffer para almacenar la fila leída (x,y,z)
    
    
    sprintf(infile, "%s.%d.hdf5", arg1, FILEid);

    file_id = H5Fopen(infile, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_id < 0) {
      fprintf(stderr, "Error al abrir el archivo HDF5: %s\n", infile);
      return -1;
    }

    dataset_id = H5Dopen2(file_id, "/PartType1/Coordinates", H5P_DEFAULT);
    dataspace_id = H5Dget_space(dataset_id);
    // Read the positions of the particles one by one and assign the BoxID
    // depending on the domain_info structure

    for(int i=0; i<N_Rpart_Total; i++){
      if(Rpart[i].FileID==FILEid){
        if(Rpart[i].BoxID==task){
          //obtener las posiciones para almacenarlas aqui
	        Particle[counter].Oid = i;
	        Particle[counter].tag = Rpart[i].tag;
          offset[0] = Rpart[i].ID; // Fila en el dataset
          offset[1] = 0;             // Inicia en la primera columna (x)
          count[0] = 1;              // Una sola fila
          count[1] = 3;              // Tres columnas (x, y, z)
          status = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
          
          if (status < 0) {
            fprintf(stderr, "Error al seleccionar la hiperslab para la fila %llu.\n", offset[0]);
            continue;
          }
          
          /// Crear un dataspace de memoria para leer la hiperslab (dimensiones 1x3)
          hsize_t dimsm[2] = {1, 3};
          memspace_id = H5Screate_simple(2, dimsm, NULL);
          status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, pos);

          if (status < 0) {
            fprintf(stderr, "Error al leer los datos de la fila %llu.\n", offset[0]);
          } else {
            // Almacenar las coordenadas en el arreglo Particle
	    if (counter==0){
	      printf("/////Posiciones primera particula/////\n");
	      printf("x = %f\n",pos[0]);
	      printf("y = %f\n",pos[1]);
	      printf("z = %f\n",pos[2]);
	    }
            Npos[0] = pos[0];
            Npos[1] = pos[1];
            Npos[2] = pos[2];

            //periodic adjustment
            if(pow(pos[0]-posa[0],2) > Plength) {
              if(posa[0] > pos[0])
                Npos[0] = pos[0] + BoxSize;
              else
                Npos[0] = pos[0] - BoxSize;
              }
          
            if(pow(pos[1]-posa[1],2) > Plength) {
              if(posa[1] > pos[1])
                Npos[1] = pos[1] + BoxSize;
              else
                Npos[1] = pos[1] - BoxSize;
              }
          
            if(pow(pos[2]-posa[2],2) > Plength) {
              if(posa[2] > pos[2])
                Npos[2] = pos[2] + BoxSize;
              else
                Npos[2] = pos[2] - BoxSize;
              }


	    Particle[counter].pos[0] = Npos[0];
            Particle[counter].pos[1] = Npos[1];
            Particle[counter].pos[2] = Npos[2];
            counter++;
	    
          }

          // Cerrar el memspace creado para esta lectura
          H5Sclose(memspace_id);
        }
      
      
      }
    }
  
    // Cerrar los identificadores de HDF5
    H5Sclose(dataspace_id);
    H5Dclose(dataset_id);
    H5Fclose(file_id);
  
  }
  
  printf("Task %d has %d particles\n",task,domain_info[task].Nparticles);
  
  Particle = (struct Vpart *) realloc(Particle,(size_t) domain_info[task].Nparticles*sizeof(struct Vpart));
  
  return domain_info[task].Nparticles;
 
  
}
