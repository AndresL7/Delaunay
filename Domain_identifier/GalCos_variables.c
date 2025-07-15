#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include "GalCos_variables.h"



char *data_prefix, *catalog_prefix;

int Npart_Total,MINIMUM_MEMBERS,NGB_MAX,FOF_PROCEDURE;
int FLAG_SUBFIND,FLAG_WRITE_INDIVIDUAL_CLUSTERS;
int MINIMUM_NSUBSTRUCT,NCLUSTERS, Npart_clustered, Npart_unclustered;

double BoxSize,REDSHIFT,OMEGA_MATTER,OMEGALAMBDA,HUBBLEPARAM,OMEGABARYON;
//flags
int FlagSfr,FlagFeedback,FlagCooling;
int NumFiles;
double COSMIC_TIME;
float Llenght,b_Link;
double GRAV_SOFT,PARTMASS;

double G_INTERNAL_UNITS,LENGHT_INTERNAL_UNITS;
double VELOCITY_INTERNAL_UNITS,MASS_INTERNAL_UNITS,TIME_INTERNAL_UNITS;
double ENERGY_INTERNAL_UNITS,DENSITY_INTERNAL_UNITS,HUBBLE_INTERNAL_UNITS;

//struct grid *GridPoint;
struct part *Particle=NULL;
struct halo *Halos=NULL;
struct COLA temp_cola;
struct COLA Cola_subhalos;
struct COLA Cola_SUBFIND;
struct COLA_float NGB_DISTcand;
struct COLA NGB_cand;
