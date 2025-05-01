#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include "GalCos_variables.h"

long int *idum,seed;
int Npart_Total;
int MINIMUM_MEMBERS,NGB_MAX,FOF_PROCEDURE;
int FLAG_SUBFIND,FLAG_WRITE_INDIVIDUAL_CLUSTERS;
int MINIMUM_NSUBSTRUCT,NCLUSTERS;
int *Mbins,NBINS,NITER;

int FlagSfr, FlagFeedback, FlagCooling, NumFiles;


float BoxSize,REDSHIFT,OMEGA_MATTER,OMEGALAMBDA,HUBBLEPARAM,OMEGABARYON;
float Llenght,COSMIC_TIME,b_Link;
float GRAV_SOFT,RMIN,RMAX;

double G_INTERNAL_UNITS,LENGHT_INTERNAL_UNITS;
double VELOCITY_INTERNAL_UNITS,MASS_INTERNAL_UNITS,TIME_INTERNAL_UNITS;
double ENERGY_INTERNAL_UNITS,DENSITY_INTERNAL_UNITS,HUBBLE_INTERNAL_UNITS;

FILE *fp_inp;

struct grid *GridPoint;
//struct part *Particle;
struct Vpart *Particle;
struct part *Rpart;
struct halo *Halos;
struct COLA temp_cola;
struct COLA Cola_subhalos;
struct COLA Cola_SUBFIND;
struct COLA_float NGB_DISTcand;
struct COLA NGB_cand;
