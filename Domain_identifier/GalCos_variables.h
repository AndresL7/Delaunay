#ifndef INIT_GALCOS_H
#define INIT_GALCOS_H

#include<stdlib.h>
#include<stdio.h>
#include<math.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_sort_double.h>

#define EMPTY_FLAG -1
#define SUCCES 0
#define FAILURE 1
#define BH_OPENING 0.7

/* Some physical and astrophysical constants (in cgs units) */

#define G_GRAVITY         6.672e-8   
#define HUBBLE            3.2407789e-18	/* in h/sec */
#define HUBBLE_TIME       3.09e+17   // in sec/h

extern int Npart_Total,MINIMUM_MEMBERS,NGB_MAX,FOF_PROCEDURE;
extern int FLAG_SUBFIND,FLAG_WRITE_INDIVIDUAL_CLUSTERS;
extern int MINIMUM_NSUBSTRUCT,NCLUSTERS; 

extern double BoxSize,REDSHIFT,OMEGA_MATTER,OMEGALAMBDA,HUBBLEPARAM,OMEGABARYON,COSMIC_TIME,PARTMASS;
extern float Llenght,b_Link;
extern double GRAV_SOFT;

extern double G_INTERNAL_UNITS,LENGHT_INTERNAL_UNITS;
extern double VELOCITY_INTERNAL_UNITS,MASS_INTERNAL_UNITS,TIME_INTERNAL_UNITS;
extern double ENERGY_INTERNAL_UNITS,DENSITY_INTERNAL_UNITS,HUBBLE_INTERNAL_UNITS;


struct tnode
{
  int NumParts;             // Number of particles in this node
  int tag;                  // tag=1 if twig, tag=0 if leaf
  int PartsInSon[8];        // Number of particles in every son node
  int index;
  int IDparticle;           // IDs of the particles in the node
  int Nsons;
  int *particlesInNode;
  int IDson;

  struct tnode **sons;      // Pointers to the 8 possible sons of the node
  float size;               // Side size of the node
  float Q;                  // Quadrupole momentum of the node
  float D;                  // Dipole momentum of the node
  float pos[3];
  float poscm[3];
  float mass;  
  struct tnode *father;
  //TREENODEPTR *father;
};

typedef struct tnode TREENODE;
typedef TREENODE *TREENODEPTR;

extern struct COLA
{
  int head;
  int tail;
  int cola_Nmembers;
  int *inputs;
}temp_cola;

extern struct COLA_float
{
  int head;
  int tail;
  int cola_Nmembers;
  float *inputs;
}NGB_DISTcand;

extern struct halo
{
  float mass;
  float pos[3];
  float vel[3];
  float Hmr;
  float Rvir,Mvir;
  int ID_CenterHalo;
  int Nmembers;
  int NDomain_particles;
  int IDcluster;
  int *Halo_particles;
  int *Domain_particles;
  //Vvir,Tvir,TotalEnergy,TotalAngularMom,Lambda_spin;
  //float Mass_fraction;
  //float TotalEk;
  //float TotalEP;
  float Radius;
  //int NumberOfSubhalos;
  //int *ID_subhalos;
}*Halos;

extern struct part
{
  float pos[3];
  float vel[3];
  float mass;
  float EP;
  int id;
  int Oid;
  int Cluster_ID;
}*Particle;


//int GalCos_load_Gadget(char *infile);
void GalCos_InitHalo(int size_init);
float distance(float xi, float yi, float zi, float xj, float yj, float zj);
void GalCos_BuildCluster(struct COLA cola,int Ncluster,int HaloIDCenter);
void GalCos_Halo_prop(struct halo *Haloinp);

int get_node_center(TREENODEPTR father,float pos[],int i);
//int get_number_of_particles(float *pos,float size,int *IDpart,TREENODEPTR *p,int *inputs,int NpartsInBox);
int get_number_of_particles(float *pos,float size,int *IDpart,int *inputs,int NpartsInBox);

int gsl_fisort(int dimension,float *fvector,int *ivector);
TREENODEPTR Build_tree(int *inputs, int NPARTICLES);

//void alloca_node(TREENODEPTR *p,struct tnode *father,float size,float *pos,int *inputs,int Nparticles, int IDson);
//void alloca_node(TREENODEPTR *p,TREENODEPTR *father,float size,float *pos,int *inputs,int Nparticles, int IDson);
void alloca_node(TREENODEPTR *p,struct tnode *father,float size,float *pos);

void GalCos_tree_walk(TREENODEPTR *p);
void GalCos_tree_Multipole(TREENODEPTR *p,TREENODEPTR walk);


#endif /* INIT_GALCOS_H */
