#include "GalCos_variables.h"

#include <gsl/gsl_sort.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sort_double.h>
#include <gsl/gsl_sort_float.h>

int GLOBAL_COUNTER=0;
float GLOBAL_MEMORY=0;

int get_node_center(TREENODEPTR father,float *pos,int i)
{
  
  /*
    Here is defined the numberin of nodes in the tree. starting with the
    lower left (i=0) growing in counter clockwise and from bottom to top
  */
  
  if(i == 0)
    {
      pos[0] = father->pos[0] - father->size/4.0;
      pos[1] = father->pos[1] - father->size/4.0;
      pos[2] = father->pos[2] - father->size/4.0;
    }
  
  if(i == 1)
    {
      pos[0] = father->pos[0] + father->size/4.0;
      pos[1] = father->pos[1] - father->size/4.0;
      pos[2] = father->pos[2] - father->size/4.0;
    }
  
  if(i == 2)
    {
      pos[0] = father->pos[0] + father->size/4.0;
      pos[1] = father->pos[1] + father->size/4.0;
      pos[2] = father->pos[2] - father->size/4.0;
    }
  
  if(i == 3)
    {
      pos[0] = father->pos[0] - father->size/4.0;
      pos[1] = father->pos[1] + father->size/4.0;
      pos[2] = father->pos[2] - father->size/4.0;
    }
  
  if(i == 4)
    {
      pos[0] = father->pos[0] - father->size/4.0;
      pos[1] = father->pos[1] - father->size/4.0;
      pos[2] = father->pos[2] + father->size/4.0;
    }
  
  if(i == 5)
    {
      pos[0] = father->pos[0] + father->size/4.0;
      pos[1] = father->pos[1] - father->size/4.0;
      pos[2] = father->pos[2] + father->size/4.0;
    }
  
  if(i == 6)
    {
      pos[0] = father->pos[0] + father->size/4.0;
      pos[1] = father->pos[1] + father->size/4.0;
      pos[2] = father->pos[2] + father->size/4.0;
    }
  
  if(i == 7)
    {
      pos[0] = father->pos[0] - father->size/4.0;
      pos[1] = father->pos[1] + father->size/4.0;
      pos[2] = father->pos[2] + father->size/4.0;
    }
  
  return 0;
}

int *auxpointer=NULL,*copy_inputs=NULL,TREE_BUILD_NPRTICLES;

//int get_number_of_particles(float *pos,float size,int *IDpart,TREENODEPTR *p,int *inputs,int NpartsInBox)
int get_number_of_particles(float *pos,float size,int *IDpart,int *inputs,int NpartsInBox)
{
  
  float xmin,xmax,ymin,ymax,zmin,zmax,Hbox;
  int counter,i,ipart,*cola_get_particles=NULL;
  //struct COLA cola_get_particles;
  //load_cola(&cola_get_particles);
  
  Hbox=size/2.0;
  
  xmin=pos[0] - Hbox;
  xmax=pos[0] + Hbox;
  
  ymin=pos[1] - Hbox;
  ymax=pos[1] + Hbox;
  
  zmin=pos[2] - Hbox;
  zmax=pos[2] + Hbox;
  
  counter=0;
  for(i=0; i<NpartsInBox; i++)
    {
      
      ipart=inputs[i];
      
      if( (Particle[ipart].pos[0] >= xmin) && (Particle[ipart].pos[0] < xmax) )
	{
	  if( (Particle[ipart].pos[1] >= ymin) && (Particle[ipart].pos[1] < ymax) )
	    {
	      if( (Particle[ipart].pos[2] >= zmin) && (Particle[ipart].pos[2] < zmax) )
		{
		  counter++;
		  *IDpart=Particle[ipart].id;
		  
		  //encola(&cola_get_particles,ipart);
		  
		  cola_get_particles = (int *) realloc(cola_get_particles,(size_t) counter*sizeof(int));
		  
		  cola_get_particles[counter-1] = ipart;

		  if((ipart != Particle[ipart].id) || (ipart != *IDpart))
		    {
		      printf("mmm un problema!!\n");
		      exit(0);
		    }
		}
	    }
	}
      
    }
  
  auxpointer = (int *) malloc((size_t) counter*sizeof(int));
  
  for(i=0; i<counter; i++)
    auxpointer[i] = cola_get_particles[i];
  
  //reinitialize_cola(&cola_get_particles);
  if(cola_get_particles != NULL)
    free(cola_get_particles);
  
  return counter;
  
}


void alloca_node(TREENODEPTR *p,struct tnode *father,float size,float *pos)
{
  
  int i,iparticle,Nparticles;
  float center[3],sizenode;
  
  if((*p) == NULL)
    {
      
      (*p)=(struct tnode *) malloc(sizeof(struct tnode));
      
      if((*p) != NULL)
	{
	  	  
	  (*p)->particlesInNode = NULL;
	  (*p)->sons = NULL;
	  
	  (*p)->index = GLOBAL_COUNTER;
	  
	  if((*p)->index == 0)
	    {
	      
	      (*p)->particlesInNode = (int *) malloc((size_t) TREE_BUILD_NPRTICLES*sizeof(int));
	      
	      (*p)->NumParts = TREE_BUILD_NPRTICLES;
	      
	      for(i=0; i<TREE_BUILD_NPRTICLES; i++)
		(*p)->particlesInNode[i] = copy_inputs[i];

	      free(copy_inputs);
	      
	    }
	  else
	    {
	      //get_number_of_particles(pos,size,&iparticle,p,(*father)->particlesInNode,(*father)->NumParts);
	      
	      (*p)->NumParts = get_number_of_particles(pos,size,&iparticle,father->particlesInNode,father->NumParts);
	      
	      (*p)->particlesInNode = (int *) malloc((size_t) (*p)->NumParts*sizeof(int));
	      
	      for(i=0; i<(*p)->NumParts; i++)
		(*p)->particlesInNode[i] = auxpointer[i];
	      
	      free(auxpointer);
	      	      
	    }
	  
	  GLOBAL_COUNTER++;
	  
	  (*p)->pos[0] = pos[0];
	  (*p)->pos[1] = pos[1];
	  (*p)->pos[2] = pos[2];
	  
	  (*p)->poscm[0] = 0.0;
	  (*p)->poscm[1] = 0.0;
	  (*p)->poscm[2] = 0.0;
	  
	  (*p)->Q = 0.0;
	  (*p)->D = 0.0;
	  
	  (*p)->size = size;
	  
	  Nparticles = (*p)->NumParts;
	  
	  /* Empty node */
	  if(Nparticles == 0)
	    {
	      (*p)->tag = EMPTY_FLAG;
	      (*p)->IDparticle = EMPTY_FLAG;
	      (*p)->mass = 0.0;
	      (*p)->Nsons = 0;
	    }
	  
	  /* leaf node... PARTICLE */
	  if(Nparticles == 1)
	    {
	      (*p)->tag = 0;
	      (*p)->IDparticle = iparticle;
	      (*p)->mass = Particle[iparticle].mass;
	      (*p)->Nsons = 0;
	      
	      (*p)->poscm[0] = Particle[iparticle].pos[0];
	      (*p)->poscm[1] = Particle[iparticle].pos[1];
	      (*p)->poscm[2] = Particle[iparticle].pos[2];
	    }
	  
	  /* Twing node */
	  if(Nparticles > 1)
	    {
	      
	      (*p)->tag=1;
	      (*p)->IDparticle=EMPTY_FLAG;
	      (*p)->mass=0.0;
	      (*p)->Nsons=8;
	      
	      (*p)->sons=(struct tnode **) malloc((*p)->Nsons*sizeof(struct tnode*));
	      
	      for(i=0; i<(*p)->Nsons; i++)
		(*p)->sons[i]=NULL;
	      
	      for(i=0; i<(*p)->Nsons; i++)
		{
		  get_node_center(*p,center,i);
		  sizenode=(*p)->size/2.0;
		  alloca_node(&((*p)->sons[i]),*p,sizenode,center);
		}
	      
	    }
	  
	}
      else
	{
	  printf("No memory available to allocate pointer\n");
	  exit(0);	    
	}
	
    }
 
}

void GalCos_tree_mass(TREENODEPTR *p)
{
  
  int i;
  
  if((*p) != NULL)
    {
      
      for(i=0; i<(*p)->Nsons; i++)
	{
	  GalCos_tree_mass(&((*p)->sons[i]));
	}
      
      if((*p)->tag == 1)
	{
	  for(i=0; i<(*p)->Nsons; i++)
	    {
	      (*p)->mass = (*p)->mass + (*p)->sons[i]->mass;
	    }
	}
    }
  
}

void GalCos_tree_CM(TREENODEPTR *p)
{
  
  int i;
  
  if((*p) != NULL)
    {
      
      for(i=0; i<(*p)->Nsons; i++)
	{
	  GalCos_tree_CM(&((*p)->sons[i]));
	}
      
      
      if((*p)->tag == 1)
	{
	  (*p)->poscm[0] = 0.0;
	  (*p)->poscm[1] = 0.0;
	  (*p)->poscm[2] = 0.0;
	  
	  for(i=0; i<(*p)->Nsons; i++)
	    {
	      (*p)->poscm[0] = (*p)->poscm[0] + ((*p)->sons[i]->poscm[0])*((*p)->sons[i]->mass);
	      (*p)->poscm[1] = (*p)->poscm[1] + ((*p)->sons[i]->poscm[1])*((*p)->sons[i]->mass);
	      (*p)->poscm[2] = (*p)->poscm[2] + ((*p)->sons[i]->poscm[2])*((*p)->sons[i]->mass);
	    }
	  
	  (*p)->poscm[0] = (*p)->poscm[0]/(*p)->mass;
	  (*p)->poscm[1] = (*p)->poscm[1]/(*p)->mass;
	  (*p)->poscm[2] = (*p)->poscm[2]/(*p)->mass;
	  
	}
      
    }
  
}

void GalCos_tree_walk(TREENODEPTR *p)
{
  
  int i;
  
  GalCos_tree_Multipole(p,*p);
  
  for(i=0; i<(*p)->Nsons; i++)
    {
      GalCos_tree_Multipole(&((*p)->sons[i]),(*p)->sons[i]);
    }
  
}

void GalCos_Free_BHtree(TREENODEPTR *p)
{
  
  int i;
  
  if(*p != NULL)
    {
      
      if((*p)->tag == 1)
	{
	  
	  for(i=0; i<(*p)->Nsons; i++)
	    GalCos_Free_BHtree(&((*p)->sons[i]));
	  
	}
      
      if((*p)->particlesInNode != NULL)
	free((*p)->particlesInNode);
            
      if((*p)->sons != NULL)
	free((*p)->sons);
      
      free((*p));
      
    }
  
}


void GalCos_tree_Multipole(TREENODEPTR *p,TREENODEPTR walk)
{
  
  int i;
  float Rcm,R3,R5,xi,yi,zi;
  
  if(p != NULL)
    {
      
      for(i=0; i<walk->Nsons; i++)
	{
	  GalCos_tree_Multipole(p,walk->sons[i]);
	}
      
      
      if(walk->tag == 0)
	{
	  
	  Rcm=distance((*p)->poscm[0],(*p)->poscm[1],(*p)->poscm[2],0,0,0);
	  R3=Rcm*Rcm*Rcm;
	  R5=R3*Rcm*Rcm;
	  
	  xi = walk->poscm[0] - (*p)->poscm[0];
	  yi = walk->poscm[1] - (*p)->poscm[1];
	  zi = walk->poscm[2] - (*p)->poscm[2];
	  
	  (*p)->D = (*p)->D - walk->mass*Dipole_summ(xi,yi,zi,(*p)->poscm[0],(*p)->poscm[1],(*p)->poscm[2],R3);
	  (*p)->Q = (*p)->Q - walk->mass*Quadrupole_summ(xi,yi,zi,(*p)->poscm[0],(*p)->poscm[1],(*p)->poscm[2],Rcm,R5);
	  
	}
      
    }
  
}


 /*
   void GalCos_tree_Multipole(TREENODEPTR *p,TREENODEPTR walk)
   {
   
   int i,ipart;
   float Rcm,R3,R5,xi,yi,zi;
  
   if(p != NULL)
   {
   
   if((*p)->index != 0)
   {
   
   Rcm=distance((*p)->poscm[0],(*p)->poscm[1],(*p)->poscm[2],0,0,0);
   R3=Rcm*Rcm*Rcm;
   R5=R3*Rcm*Rcm;
   
   for(i=0; i<(*p)->NumParts; i++)
   {
   ipart = (*p)->particlesInNode[i];
   
   xi = Particle[ipart].pos[0] - (*p)->poscm[0];
   yi = Particle[ipart].pos[1] - (*p)->poscm[1];
   zi = Particle[ipart].pos[2] - (*p)->poscm[2];
   
   (*p)->D = (*p)->D - Particle[ipart].mass*Dipole_summ(xi,yi,zi,(*p)->poscm[0],(*p)->poscm[1],(*p)->poscm[2],R3);
   (*p)->Q = (*p)->Q - Particle[ipart].mass*Quadrupole_summ(xi,yi,zi,(*p)->poscm[0],(*p)->poscm[1],(*p)->poscm[2],Rcm,R5);
   }
   
   }
   
   for(i=0; i<(*p)->Nsons; i++)
   {
   GalCos_tree_Multipole(&(*p),walk);
   }
   
   }
  
   }
 */

void GalCos_tree_force(int ipart,TREENODEPTR p)
{
  
  float dist;
  int i;
  
  if(p != NULL)
    {
      
      if(p->tag != EMPTY_FLAG)
	{
	  
	  dist=distance(Particle[ipart].pos[0],Particle[ipart].pos[1],Particle[ipart].pos[2],
			p->poscm[0],p->poscm[1],p->poscm[2]);
	  
	  if((p->size/dist) <= BH_OPENING ) // Particle-Node
	    {
	      
	      /*
	       * This line is just including MONOPOLE contributions. I
	       * have to check the computaion of dipole and quadrupole
	       * because I think is not fine !! ^o^
	       */
	      
	      //Particle[ipart].EP = Particle[ipart].EP - G_INTERNAL_UNITS*Particle[ipart].mass*(p->mass/GRAV_SOFT)*grav_soft_spline(dist,GRAV_SOFT);
	      //Particle[ipart].EP = Particle[ipart].EP - G_INTERNAL_UNITS*Particle[ipart].mass*p->mass/dist;
	      
	      /* FULL POTENTIAL!! */
	      //Particle[ipart].EP = Particle[ipart].EP - G_INTERNAL_UNITS*Particle[ipart].mass*((p->mass/GRAV_SOFT)*grav_soft_spline(dist,GRAV_SOFT) + p->D + p->Q );
	      Particle[ipart].EP = Particle[ipart].EP - G_INTERNAL_UNITS*((p->mass/GRAV_SOFT)*grav_soft_spline(dist,GRAV_SOFT) + p->D + p->Q );
	      
	    }
	  else /* Opening the node -- walk across subnodes inside this node */
	    {
	      
	      for(i=0; i<p->Nsons; i++)
		{
		  GalCos_tree_force(ipart,p->sons[i]);
		}
	      
	    }
	  
	}
       
    }

}

void GalCos_tree_print(TREENODEPTR p)
{
  
  int i;
  //float xmin,xmax,ymin,ymax,zmin,zmax,Hbox;
  //FILE *pf;
  
  if(p != NULL)
    {
      
      if(p->tag == 1)
	{
	  //if(p->father != NULL)
	    printf("%f %f\n",p->NumParts*Particle[0].mass,p->mass);
	    //printf("%d %d %d %f %f %f %d\n",p->tag,p->index,p->father->index,p->pos[0],p->pos[1],p->pos[2],p->NumParts);
	  
	  //printf("%d %f %f %f %f %f %f %d %f %g %g\n",p->father->index,
	  //p->pos[0],p->pos[1],p->pos[2],p->poscm[0],p->poscm[1],p->poscm[2],p->NumParts,p->mass,p->Q,p->D);
	  //else
	  //printf("%d %d -1 %f %f %f %d\n",p->tag,p->index,p->pos[0],p->pos[1],p->pos[2],p->NumParts);
	  //printf("-1 %f %f %f %f %f %f %d %f %g %g\n",p->index,p->father->index,p->pos[0],
	  //p->pos[1],p->pos[2],p->poscm[0],p->poscm[1],p->poscm[2],p->NumParts,p->mass,p->Q,p->D);
	}
      
      
      if(p->NumParts > 1)
	{
	  
	  for(i=0; i<p->Nsons; i++)
	    {
	      GalCos_tree_print(p->sons[i]);
	    }
	}
      
      /*
      pf=fopen("plotar.gpl","a");
      
      Hbox=p->size/2.0;
      
      xmin=p->pos[0] - Hbox;
      xmax=p->pos[0] + Hbox;
      
      ymin=p->pos[1] - Hbox;
      ymax=p->pos[1] + Hbox;
      
      zmin=p->pos[2] - Hbox;
      zmax=p->pos[2] + Hbox;
      
      //fprintf(pf,"set ticslevel 0\n");
      fprintf(pf,"set arrow from %.3f,%.3f,%.3f to %.3f,%.3f,%.3f nohead\n",xmin,ymin,zmin,xmax,ymin,zmin);
      fprintf(pf,"set arrow from %.3f,%.3f,%.3f to %.3f,%.3f,%.3f nohead\n",xmin,ymin,zmin,xmin,ymax,zmin);
      fprintf(pf,"set arrow from %.3f,%.3f,%.3f to %.3f,%.3f,%.3f nohead\n",xmin,ymin,zmin,xmin,ymin,zmax);
      fprintf(pf,"set arrow from %.3f,%.3f,%.3f to %.3f,%.3f,%.3f nohead\n",xmax,ymin,zmin,xmax,ymax,zmin);
      fprintf(pf,"set arrow from %.3f,%.3f,%.3f to %.3f,%.3f,%.3f nohead\n",xmax,ymin,zmin,xmax,ymin,zmax);
      fprintf(pf,"set arrow from %.3f,%.3f,%.3f to %.3f,%.3f,%.3f nohead\n",xmax,ymax,zmin,xmax,ymax,zmax);
      fprintf(pf,"set arrow from %.3f,%.3f,%.3f to %.3f,%.3f,%.3f nohead\n",xmin,ymax,zmin,xmin,ymax,zmax);
      fprintf(pf,"set arrow from %.3f,%.3f,%.3f to %.3f,%.3f,%.3f nohead\n",xmin,ymax,zmin,xmax,ymax,zmin);
      fprintf(pf,"set arrow from %.3f,%.3f,%.3f to %.3f,%.3f,%.3f nohead\n",xmin,ymin,zmax,xmax,ymin,zmax);
      fprintf(pf,"set arrow from %.3f,%.3f,%.3f to %.3f,%.3f,%.3f nohead\n",xmin,ymin,zmax,xmin,ymax,zmax);
      fprintf(pf,"set arrow from %.3f,%.3f,%.3f to %.3f,%.3f,%.3f nohead\n",xmin,ymax,zmax,xmax,ymax,zmax);
      fprintf(pf,"set arrow from %.3f,%.3f,%.3f to %.3f,%.3f,%.3f nohead\n",xmax,ymin,zmax,xmax,ymax,zmax);
      //fprintf(pf,"pause -1\n");
      
      fclose(pf);
      */
    }
    
}

TREENODEPTR Build_tree(int *inputs, int NPARTICLES)
{
  
  float xmin,xmax,ymin,ymax,zmin,zmax,BOXSIZE,pos[3]; 
  float xsize,ysize,zsize,BoxXCenter,BoxYCenter,BoxZCenter;
  int i,ipart;
  TREENODEPTR root=NULL;
  
  copy_inputs = (int *) malloc((size_t) NPARTICLES*sizeof(int));

  xmin=0; ymin=0; zmin=0;
  xmax=0; xmax=0; ymax=0;
  
  for(i=0; i<NPARTICLES; i++)
    {
      ipart = inputs[i];
      
      // min
      if(Particle[ipart].pos[0] < xmin)
        xmin = Particle[ipart].pos[0];

      if(Particle[ipart].pos[1] < ymin)
        ymin = Particle[ipart].pos[1];
      
      if(Particle[ipart].pos[2] < zmin)
        zmin = Particle[ipart].pos[2];
      
      // max
      
      if(Particle[ipart].pos[0] > xmax)
        xmax = Particle[ipart].pos[0];
      
      if(Particle[ipart].pos[1] > ymax)
        ymax = Particle[ipart].pos[1];
      
      if(Particle[ipart].pos[2] > zmax)
        zmax = Particle[ipart].pos[2];
      
      copy_inputs[i] = inputs[i];
      
    }
  
  xsize = (xmax - xmin);
  ysize = (ymax - ymin);
  zsize = (zmax - zmin);
  
  BoxXCenter = xmin + xsize/2.0;
  BoxYCenter = ymin + ysize/2.0;
  BoxZCenter = zmin + zsize/2.0;
  
  if((xsize >= ysize) && (xsize >= zsize))
    BOXSIZE=xsize;
  else if((ysize >= xsize) && (ysize >= zsize))
    BOXSIZE=ysize;
  else if((zsize >= xsize) && (zsize >= ysize))
    BOXSIZE=zsize;
  else
    {
      printf("problem computing the root cell of the tree %f %f %f %f\n",xsize,ysize,zsize,BOXSIZE);
      exit(0);
    }

  
  BOXSIZE = BOXSIZE + 2.0*BOXSIZE/pow(1.0*NPARTICLES,1.0/3.0);
  
  pos[0] = BoxXCenter;
  pos[1] = BoxYCenter;
  pos[2] = BoxZCenter;
  
  TREE_BUILD_NPRTICLES = NPARTICLES;

  alloca_node(&root,NULL,BOXSIZE,pos);
  //printf("Tree is OK!\n");
    
  GalCos_tree_mass(&root);
  GalCos_tree_CM(&root);
    
  GalCos_tree_walk(&root);
  
  //printf("computing multipoles\n");
  //GalCos_tree_Multipole(&root,root);
  //printf("computed multipoles\n");
  //GalCos_tree_print(root);
  
  GLOBAL_COUNTER=0;
  //GLOBAL_MEMORY=0;

  return root;
  
}

