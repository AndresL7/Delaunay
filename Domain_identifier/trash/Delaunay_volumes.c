#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include <gsl/gsl_linalg.h>

#include "GalCos_SPH_density.c"

int compute_spherical_volume_domain(int IDhalo)
{
  
  int i,*SPH_NGB,ipart,j,jpart,counter;
  float *SPH_NGB_DISTANCES,Radius;
  TREENODEPTR root=NULL;
  FILE *pf;
  
  pf=fopen("plota_spheres.dat","w");
  fprintf(pf,"set parametric\n");
  
  if(Halos[IDhalo].NDomain_particles > NGB_MAX)
    {
      
      /* Building the BH for NGB search */
      root = Build_tree(Halos[IDhalo].Domain_particles,Halos[IDhalo].NDomain_particles);
      
      SPH_NGB = (int *) malloc(NGB_MAX*sizeof(int));
      SPH_NGB_DISTANCES = (float *) malloc(NGB_MAX*sizeof(float));
      
      for(i=0; i<Halos[IDhalo].NDomain_particles; i++)
	{
	  ipart = Halos[IDhalo].Domain_particles[i];
	  GalCos_SPH_NGB(ipart,SPH_NGB,SPH_NGB_DISTANCES,Halos[IDhalo].NDomain_particles,&root,Halos[IDhalo].Domain_particles);
	  
	  Radius = SPH_NGB_DISTANCES[0]/4.0;
	  
	  fprintf(pf,"replot (%f-%f*cos(u)*cos(v)),(%f-%f*cos(u)*sin(v)),(%f-%f*sin(u)) not w l 3\n",Particle[ipart].pos[0],
		  Radius,Particle[ipart].pos[1],Radius,Particle[ipart].pos[2],Radius);
	  
	  Particle[ipart].Volume = (4.0/3.0)*M_PI*Radius*Radius*Radius;
	  
	  if(Particle[ipart].Volume < 0)
	    {
	      Particle[ipart].Volume=0;
	      printf("que putas esta pasando?\n");
	      getchar();
	    }
	  
	}
      
      GalCos_Free_BHtree(&root);
      
      free(SPH_NGB);
      free(SPH_NGB_DISTANCES);
      
    }
  else
    {
      
      SPH_NGB = (int *) malloc((Halos[IDhalo].NDomain_particles-1)*sizeof(int));
      SPH_NGB_DISTANCES = (float *) malloc((Halos[IDhalo].NDomain_particles-1)*sizeof(float));
      
      for(i=0; i<Halos[IDhalo].NDomain_particles; i++)
	{
	  
	  ipart = Halos[IDhalo].Domain_particles[i];
	  counter=0;
	  
	  for(j=0; j<Halos[IDhalo].NDomain_particles; j++)
	    {
	      
	      jpart = Halos[IDhalo].Domain_particles[j];
	      
	      if(ipart != jpart)
		{
		  SPH_NGB_DISTANCES[counter] = distance(Particle[ipart].pos[0],Particle[ipart].pos[1],Particle[ipart].pos[2],
							Particle[jpart].pos[0],Particle[jpart].pos[1],Particle[jpart].pos[2]);
		  SPH_NGB[counter]=jpart;
		  counter++;
		}
	    }
	  
	  gsl_fisort(counter,SPH_NGB_DISTANCES,SPH_NGB);
	  
	  Radius = SPH_NGB_DISTANCES[0]/4.0;
	  
	  fprintf(pf,"replot (%f-%f*cos(u)*cos(v)),(%f-%f*cos(u)*sin(v)),(%f-%f*sin(u)) not w l 3\n",Particle[ipart].pos[0],
		  Radius,Particle[ipart].pos[1],Radius,Particle[ipart].pos[2],Radius);
	  
	  Particle[ipart].Volume = (4.0/3.0)*M_PI*Radius*Radius*Radius;

	  if(Particle[ipart].Volume < 0)
	    {
	      printf("que putas esta pasando?\n");
	      getchar();
	      Particle[ipart].Volume=0;
	    }
	  
	}

      free(SPH_NGB);
      free(SPH_NGB_DISTANCES);
      
    }
  
  fclose(pf);
  
  return 0;
  
}

int compute_spherical_volume_halo(int IDhalo)
{
  
  int i,*SPH_NGB,ipart;
  float *SPH_NGB_DISTANCES,Radius;
  TREENODEPTR root=NULL;
  
  /* Building the BH for NGB search */
  root = Build_tree(Halos[IDhalo].Halo_particles,Halos[IDhalo].Nmembers);
  
  SPH_NGB = (int *) malloc(NGB_MAX*sizeof(int));
  SPH_NGB_DISTANCES = (float *) malloc(NGB_MAX*sizeof(float));
  
  for(i=0; i<Halos[IDhalo].Nmembers; i++)
    {
      
      ipart = Halos[IDhalo].Halo_particles[i];
      GalCos_SPH_NGB(ipart,SPH_NGB,SPH_NGB_DISTANCES,Halos[IDhalo].Nmembers,&root,Halos[IDhalo].Halo_particles);
      
      Radius = SPH_NGB_DISTANCES[0]/2.0;
      Particle[ipart].Volume = (4.0/3.0)*M_PI*Radius*Radius*Radius;

      if(Particle[ipart].Volume < 0)
	{
	  printf("que putas esta pasando?\n");
	  //getchar();
	  Particle[ipart].Volume=0;
	}
      
    }
  
  GalCos_Free_BHtree(&root);
  
  free(SPH_NGB);
  free(SPH_NGB_DISTANCES);
  
  return 0;
  
}



int fill_matrix(double *distances,int **facets,double **coords,int ifacet)
{
  ////////////////////////////
  
  distances[0] = 0.0;
  distances[1] = 1.0;    
  distances[2] = 1.0;    
  distances[3] = 1.0;    
  distances[4] = 1.0;    
  
  ///////////////////////////
  
  distances[5] = 1.0;    
  distances[6] = 0.0;
  //1,2//
  distances[7] = distance(coords[facets[ifacet][0]][0],coords[facets[ifacet][0]][1],coords[facets[ifacet][0]][2],
  coords[facets[ifacet][1]][0],coords[facets[ifacet][1]][1],coords[facets[ifacet][1]][2]);
  distances[7] = distances[7]*distances[7];
  
  //1,3//
  distances[8] = distance(coords[facets[ifacet][0]][0],coords[facets[ifacet][0]][1],coords[facets[ifacet][0]][2],
  coords[facets[ifacet][2]][0],coords[facets[ifacet][2]][1],coords[facets[ifacet][2]][2]);
  distances[8] = distances[8]*distances[8];
  
  //1,4//
  distances[9] = distance(coords[facets[ifacet][0]][0],coords[facets[ifacet][0]][1],coords[facets[ifacet][0]][2],
  coords[facets[ifacet][3]][0],coords[facets[ifacet][3]][1],coords[facets[ifacet][3]][2]);
  distances[9] = distances[9]*distances[9];
  
  //////////////////////////
  
  distances[10] = 1.0;
  //2,1//
  distances[11] = distances[7];
    
  distances[12] = 0.0;
  
  //2,3//
  distances[13] = distance(coords[facets[ifacet][1]][0],coords[facets[ifacet][1]][1],coords[facets[ifacet][1]][2],
			   coords[facets[ifacet][2]][0],coords[facets[ifacet][2]][1],coords[facets[ifacet][2]][2]);
  distances[13] = distances[13]*distances[13];
  
  //2,4//
  distances[14] = distance(coords[facets[ifacet][1]][0],coords[facets[ifacet][1]][1],coords[facets[ifacet][1]][2],
			   coords[facets[ifacet][3]][0],coords[facets[ifacet][3]][1],coords[facets[ifacet][3]][2]);
  distances[14] = distances[14]*distances[14];
  
  ///////////////////////////
    
  distances[15] = 1.0;
  distances[16] = distances[8];
  distances[17] = distances[13]; 
  distances[18] = 0.0;
  
  //3,4//
  distances[19] = distance(coords[facets[ifacet][2]][0],coords[facets[ifacet][2]][1],coords[facets[ifacet][2]][2],
			   coords[facets[ifacet][3]][0],coords[facets[ifacet][3]][1],coords[facets[ifacet][3]][2]);
  distances[19] = distances[19]*distances[19];
  
  ///////////////////////////    
  
  distances[20] = 1.0;
  distances[21] = distances[9]; 
  distances[22] = distances[14]; 
  distances[23] = distances[19]; 
  distances[24] = 0.0;
  
  return 0;
  
}

double compute_det(double *data)
{
  int s;
  double Det;
  
  gsl_matrix_view m = gsl_matrix_view_array(data,5,5);
  gsl_permutation * p = gsl_permutation_alloc (5);
  gsl_linalg_LU_decomp (&m.matrix, p, &s);
  
  Det=gsl_linalg_LU_det (&m.matrix,s);
  gsl_permutation_free (p);
  
  return Det;
  
}


int print_delaunuay(int **facets, double **coords, int i, int j)
{
  
  FILE *pf;
  
  pf=fopen("plotador.gpl","w");
  
  fprintf(pf,"set arrow from %f,%f,%f to %f,%f,%f nohead\n",coords[facets[i][0]][0],coords[facets[i][0]][1],
	  coords[facets[i][0]][2],coords[facets[i][1]][0],coords[facets[i][1]][1],coords[facets[i][1]][2]);
  
  fprintf(pf,"set arrow from %f,%f,%f to %f,%f,%f nohead\n",coords[facets[i][0]][0],coords[facets[i][0]][1],
	  coords[facets[i][0]][2],coords[facets[i][2]][0],coords[facets[i][2]][1],coords[facets[i][2]][2]);
  
  fprintf(pf,"set arrow from %f,%f,%f to %f,%f,%f nohead\n",coords[facets[i][0]][0],coords[facets[i][0]][1],
	  coords[facets[i][0]][2],coords[facets[i][3]][0],coords[facets[i][3]][1],coords[facets[i][3]][2]);
    
  fprintf(pf,"set arrow from %f,%f,%f to %f,%f,%f nohead\n",coords[facets[i][1]][0],coords[facets[i][1]][1],
	  coords[facets[i][1]][2],coords[facets[i][3]][0],coords[facets[i][3]][1],coords[facets[i][3]][2]);
  
  fprintf(pf,"set arrow from %f,%f,%f to %f,%f,%f nohead\n",coords[facets[i][1]][0],coords[facets[i][1]][1],
	  coords[facets[i][1]][2],coords[facets[i][2]][0],coords[facets[i][2]][1],coords[facets[i][2]][2]);
  
  fprintf(pf,"set arrow from %f,%f,%f to %f,%f,%f nohead\n",coords[facets[i][2]][0],coords[facets[i][2]][1],
	  coords[facets[i][2]][2],coords[facets[i][3]][0],coords[facets[i][3]][1],coords[facets[i][3]][2]);
  
  fprintf(pf,"set label 'Here' at %f,%f,%f\n",coords[facets[i][j]][0],coords[facets[i][j]][1],coords[facets[i][j]][2]);
  
  fclose(pf);
  
  return 0;
}


int Delaunay_volumes(int IDhalo)
{
    
  int i,j,**facets,aux,k,SIZEMATRIX,NDATS,NFACETS,ipart,counter,NDIM=3;
  double **coords,*distances,*Volumenes,Vi;
  char buff[80];
  FILE *pf;
  
  NDATS=Halos[IDhalo].NDomain_particles;
  
  //Printing file with coords for the qhull task//
  
  system("rm data");
  if((pf=fopen("data","w")) == NULL)
    {
      printf("I can not open file data\n");
      exit(0);
    }
  
  fprintf(pf,"%d\n",NDIM);
  fprintf(pf,"%d\n",Halos[IDhalo].NDomain_particles);
  
  for(i=0; i<NDATS; i++)
    {
      ipart=Halos[IDhalo].Domain_particles[i];
      fprintf(pf,"%g %g %g\n",Particle[ipart].pos[0],Particle[ipart].pos[1],
	      Particle[ipart].pos[2]);
    }
  fclose(pf);
  
  //system("rm delaunay_coors");
  sprintf(buff,"%s","qdelaunay Fv QJ < data > delaunay_coors");
  system(buff);
  printf("%s\n",buff);
  printf("qdelaunay done...\n");
  
  // Allocating space for facets
  
  if((pf=fopen("delaunay_coors","r")) == NULL)
  {
      printf("I can not open delaunay_coors file\n");
      exit(0);
  }
  
  fscanf(pf,"%d",&NFACETS);
  
  facets=(int **) malloc((size_t) NFACETS*sizeof(int *));
  
  for(i=0; i<NFACETS; i++)
    facets[i]=(int *) malloc((size_t) 4*sizeof(int));
  
  // Allocating space for coords
  coords=(double **) malloc((size_t) NDATS*sizeof(double *));
  
  for(i=0; i<NDATS; i++)
    coords[i]=(double *) malloc((size_t) 3*sizeof(double));
  
  Volumenes=(double *) malloc((size_t) NFACETS*sizeof(double));
  
  // Allocating space for matrix distances
  SIZEMATRIX=5;
  distances=(double *) malloc((size_t) SIZEMATRIX*SIZEMATRIX*sizeof(double));
  
  for(i=0; i<NFACETS; i++)
    {
      fscanf(pf,"%d %d %d %d %d",&aux,&facets[i][0],&facets[i][1],
	     &facets[i][2],&facets[i][3]);
    }
  fclose(pf);
    
  for(i=0; i<NDATS; i++)
    {
      ipart = Halos[IDhalo].Domain_particles[i];
      
      coords[i][0] = Particle[ipart].pos[0];
      coords[i][1] = Particle[ipart].pos[1];
      coords[i][2] = Particle[ipart].pos[2];
      
      Particle[ipart].Volume=0;
    }
    
    
  for(i=0; i<NFACETS; i++)
  {
      
      fill_matrix(distances,facets,coords,i);
      
      Vi=compute_det(distances);
      
      if(Vi < 0)
      {
	Vi=0.0;
      }
      
      Volumenes[i] = sqrt(Vi/288.0);
      
  }
  
  //#pragma omp parallel for
  for(k=0; k<NDATS; k++)
    {
      
      ipart=Halos[IDhalo].Domain_particles[k];
      Particle[ipart].Volume=0.0;
      counter=0;
            
      for(i=0; i<NFACETS; i++)
	{
	  
	  for(j=0; j<4; j++)
	    {
	      
	      //print_delaunuay(facets,coords,i,j);
	      if(facets[i][j] == k)
		{
		  Particle[ipart].Volume = Particle[ipart].Volume + Volumenes[i];
		  counter++;
		}
	      
	  }
	  
	  if(counter > 0)
	    Particle[ipart].Volume=Particle[ipart].Volume/4.0;
	  else
	    {
	      Particle[ipart].Volume=0;
	      printf("Domain ahhh... no jodas!! %d\n",ipart);
	      exit(0);
	    }
	  
	}
      
    }

  for(i=0; i<NFACETS; i++)
    free(facets[i]);

  free(facets);
  
  for(i=0; i<NDATS; i++)
    free(coords[i]);

  free(coords);
  free(distances);
  free(Volumenes);

  return 0;
  
}


int Delaunay_volumes_halo(int IDhalo)
{
  
    int i,j,**facets,aux,k,SIZEMATRIX,NDATS,NFACETS,ipart,counter;
    double **coords,*distances,*Volumenes, Vi;
    char buff[80];
    FILE *pf;
    
    NDATS=Halos[IDhalo].Nmembers;
    
    //Printing file with coords for the qhull task
    
    if((pf=fopen("data","w")) == NULL)
    {
	printf("I can not open file data\n");
	exit(0);
    }
    
    fprintf(pf,"3\n");
    fprintf(pf,"%d\n",Halos[IDhalo].Nmembers);
    
    for(i=0; i<NDATS; i++)
    {
	ipart=Halos[IDhalo].Halo_particles[i];
	fprintf(pf,"%g %g %g\n",Particle[ipart].pos[0],Particle[ipart].pos[1],
		Particle[ipart].pos[2]);
    }
    fclose(pf);
    
    system("rm delaunay_coors");
    sprintf(buff,"%s","qdelaunay Fv QJ < data > delaunay_coors");
    system(buff);
    
    // Allocating space for facets
    
    if((pf=fopen("delaunay_coors","r")) == NULL)
    {
	printf("I can not open delaunay_coors file\n");
	exit(0);
    }
    
    fscanf(pf,"%d",&NFACETS);
    
    facets=(int **) malloc((size_t) NFACETS*sizeof(int *));
    
    for(i=0; i<NFACETS; i++)
	facets[i]=(int *) malloc((size_t) 4*sizeof(int));
    
    // Allocating space for coords
    coords=(double **) malloc((size_t) NDATS*sizeof(double *));
    
    for(i=0; i<NDATS; i++)
	coords[i]=(double *) malloc((size_t) 3*sizeof(double));
    
    Volumenes=(double *) malloc((size_t) NFACETS*sizeof(double));
    
    // Allocating space for matrix distances
    SIZEMATRIX=5;
    distances=(double *) malloc((size_t) SIZEMATRIX*SIZEMATRIX*sizeof(double));
    
    for(i=0; i<NFACETS; i++)
    {
	fscanf(pf,"%d %d %d %d %d",&aux,&facets[i][0],&facets[i][1],
	     &facets[i][2],&facets[i][3]);
    }
    fclose(pf);
    
    for(i=0; i<NDATS; i++)
    {
	ipart=Halos[IDhalo].Halo_particles[i];
	
	coords[i][0] = Particle[ipart].pos[0];
	coords[i][1] = Particle[ipart].pos[1];
	coords[i][2] = Particle[ipart].pos[2];
	Particle[ipart].Volume=0;
    }
    
  for(i=0; i<NFACETS; i++)
  {
      
      fill_matrix(distances,facets,coords,i);
      
      Vi=compute_det(distances);
                  
      if(Vi < 0)
	{
	  Vi=0;
	}
      
      Volumenes[i] = sqrt(Vi/288.0);
      
  }
  
  for(k=0; k<NDATS; k++)
    {
      
      ipart=Halos[IDhalo].Halo_particles[k];
      Particle[ipart].Volume=0;
      counter=0;
      for(i=0; i<NFACETS; i++)
	{
	  
	  for(j=0; j<4; j++)
	    {
	      //print_delaunuay(facets,coords,i,j);
              if(facets[i][j] == k)
		{
		  Particle[ipart].Volume = Particle[ipart].Volume + Volumenes[i];
		  counter++;
		}
	      
	  }
	  
	  if(counter > 0)
	    Particle[ipart].Volume=Particle[ipart].Volume/(4.0);
	  else
	    {
	      Particle[ipart].Volume=0;
	      printf("Halo ahhh... no jodas!! %d \n",ipart);
	      exit(0);
	    }
	  
	}
      
    }
  
  for(i=0; i<NFACETS; i++)
    free(facets[i]);
  
  free(facets);
  
  for(i=0; i<NDATS; i++)
    free(coords[i]);
  
  free(coords);
  free(distances);
  free(Volumenes);
    
  return 0;
  
}



