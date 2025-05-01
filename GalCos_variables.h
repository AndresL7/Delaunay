#ifndef INIT_GALCOS_H
#define INIT_GALCOS_H

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include <hdf5.h>

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

extern long int *idum,seed;
extern int Npart_Total;
extern int MINIMUM_MEMBERS,NGB_MAX,FOF_PROCEDURE;
extern int FLAG_SUBFIND,FLAG_WRITE_INDIVIDUAL_CLUSTERS;
extern int MINIMUM_NSUBSTRUCT,NCLUSTERS; 
extern int *Mbins,NBINS,NITER;

extern float BoxSize,REDSHIFT,OMEGA_MATTER,OMEGALAMBDA,HUBBLEPARAM,OMEGABARYON;
extern float Llenght,COSMIC_TIME,b_Link;
extern float GRAV_SOFT,RMIN,RMAX;

extern double G_INTERNAL_UNITS,LENGHT_INTERNAL_UNITS;
extern double VELOCITY_INTERNAL_UNITS,MASS_INTERNAL_UNITS,TIME_INTERNAL_UNITS;
extern double ENERGY_INTERNAL_UNITS,DENSITY_INTERNAL_UNITS,HUBBLE_INTERNAL_UNITS;

extern FILE *fp_inp;

struct tnode
{
  int NumParts;             // Number of particles in this node
  int tag;                  // tag=1 if twig, tag=0 if leaf
  int index;
  int IDparticle;           // IDs of the particles in the node
  int *particlesInNode;
  int Nsons;
  
  struct tnode **sons;      // Pointers to the 8 possible sons of the node
  float size;               // Side size of the node
  float Q;                  // Quadrupole momentum of the node
  float D;                  // Dipole momentum of the node
  float pos[3];
  float poscm[3];
  float mass;  
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
}NGB_DISTcand;;

extern struct halo
{
  float mass;
  float pos[3];
  float Bpos[3];
  float vel[3];
  float Mass_fraction;
  float Hmr;
  float TotalEk;
  float TotalEP;
  float Radius;
  int ID_CenterHalo;
  int Nmembers;
  int NDomain_particles;
  int IDcluster;
  int *Halo_particles;
  int *Domain_particles;
  int NumberOfSubhalos;
  int *ID_subhalos;
  float Rvir,Mvir,Vvir,Tvir,TotalEnergy,TotalAngularMom,Lambda_spin;
}*Halos;

/*
extern struct part
{
  float pos[3];
  double Volume;
  int Oid;
  short int tag;
}*Particle;
*/

extern struct part
{
  float pos[3];
  double Volume;
  int Oid;
  short int tag;
  int BoxID;
  int ID;
  short int FileID;
}*Rpart;


extern struct Vpart
{
  double pos[3];
  double Volume;
  int Oid;
  short int tag;
  int BoxID;
  int ID;
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
void alloca_node(TREENODEPTR *p,struct tnode *father,float size,float *pos);

void GalCos_tree_walk(TREENODEPTR *p);
void GalCos_tree_Multipole(TREENODEPTR *p,TREENODEPTR walk);


#endif /* INIT_GALCOS_H */
