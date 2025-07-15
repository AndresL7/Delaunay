#include <gsl/gsl_sort.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sort_double.h>
#include <gsl/gsl_sort_float.h>

#include "GalCos_variables.h"

float distance_eps(float xi, float yi, float zi, float xj, float yj, float zj, float eps)
{
  
  float dist,x,y,z;
  
  x=(xi-xj)*(xi-xj);
  y=(yi-yj)*(yi-yj);
  z=(zi-zj)*(zi-zj);
  
  dist=sqrt( x + y + z + eps*eps);
  
  return dist;
  
}


float distance(float xi, float yi, float zi, float xj, float yj, float zj)
{
  
  float dist,x,y,z;
  
  x=(xi-xj)*(xi-xj);
  y=(yi-yj)*(yi-yj);
  z=(zi-zj)*(zi-zj);
  
  dist=sqrt( x + y + z );
  
  return dist;
  
}


int GalCos_interp_dist(void)
{

  Llenght=b_Link*BoxSize/pow(1.0*Npart_Total,1.0/3.0);
  printf("The linking lenght is %f\n",Llenght);
  
  return 0;
    
}


// organizo en orden creciente de fvector, indexando en ivector
int gsl_fisort(int dimension,float *fvector,int *ivector)
{
  
  int *aux_ivector,i;
  float *aux_fvector;
  size_t *P;

  P= (size_t *) malloc((size_t) dimension*sizeof(size_t));
  aux_fvector=(float *) malloc(dimension*sizeof(float));
  aux_ivector=(int *) malloc(dimension*sizeof(int));
  
  
  //#pragma omp parallel for
  for(i=0; i<dimension; i++)
    {
      aux_ivector[i]=ivector[i];
      aux_fvector[i]=fvector[i];
    }
  
  gsl_sort_float_index(P,fvector,1,dimension);
  
  //#pragma omp parallel for
  for(i=0; i<dimension; i++)
    {
      fvector[i]=aux_fvector[P[i]];
      ivector[i]=aux_ivector[P[i]];
    }

  free(P);
  free(aux_fvector);
  free(aux_ivector);

  return 0;

}


// organizo en orden creciente de fvector, indexando en ivector
int gsl_int_int_sort(int dimension,int *fvector,int *ivector)
{
  
  int *aux_ivector,i;
  int *aux_fvector;
  size_t *P;

  P=(size_t *) malloc((size_t) dimension*sizeof(size_t));
  if(P == NULL){
    printf("Allocation error routines:81\n");
    exit(0);
  }
    
  aux_fvector=(int *) malloc((size_t) dimension*sizeof(int));
  if(aux_fvector == NULL){
    printf("Allocation error routines:81\n");
    exit(0);
  }
  
  aux_ivector=(int *) malloc((size_t) dimension*sizeof(int));
  if(aux_ivector == NULL){
    printf("Allocation error routines:81\n");
    exit(0);
  }
 
  for(i=0; i<dimension; i++)
    {
      aux_ivector[i]=ivector[i];
      aux_fvector[i]=fvector[i];
    }

  gsl_sort_int_index(P,fvector,1,dimension);
  
  for(i=0; i<dimension; i++)
    {
      fvector[i]=aux_fvector[P[i]];
      ivector[i]=aux_ivector[P[i]];
    }

  free(P);
  free(aux_fvector);
  free(aux_ivector);

  return 0;

}

double compute_mean_mol_weight(float X, float Y, float Z)
{

  double mu;
  mu=2.0*X + (3.0/4.0)*Y + 0.5*Z;
  mu=1.0/mu;
  return mu;
  
}


double Cosmic_Time(double redshift)
{
  
  double cosmictime;
  
  //cosmictime=(2.0/3.0)*(HUBBLE_TIME/TIME_INTERNAL_UNITS)/0.73;
  cosmictime=(2.0/3.0)*(HUBBLE_TIME/31536000)/0.73;
  cosmictime=cosmictime/pow((1.0 + redshift),1.5);
  
  return cosmictime;// regresa el tiempo en años
  
}

float Virial_criterion(void)
{

  float EZsquare,OMEGA_Z,X,deltaD,rhocrit_z,value;
  float HUBBLE_PARAMETER=0.73;
  
  EZsquare=OMEGA_MATTER*pow((1.0 + REDSHIFT),3) + (1.0-OMEGA_MATTER);
  
  OMEGA_Z = OMEGA_MATTER*pow((1.0+REDSHIFT),3)/EZsquare;
  X=OMEGA_Z - 1.0;
  
  deltaD=18.0*M_PI*M_PI + 82.0*X - 39.0*X*X;
  
  rhocrit_z=(3.0/(8.0*M_PI*G_INTERNAL_UNITS))*10000*pow(HUBBLE_PARAMETER,2)*EZsquare;
  
  value=deltaD*rhocrit_z;

  return value;
  
}


long int tiempoc(void)
{
  struct tm *hora;
  time_t lt;
  long int segundos;
  
  lt = time(NULL);
  hora=localtime(&lt);
  segundos = hora->tm_sec + hora->tm_min * 60 + hora->tm_hour*3600;
  segundos = segundos + hora->tm_yday*86400;
  
  return segundos;
}


float Quadrupole_summ(float xi,float yi,float zi,float x,float y,float z,float R,float R5)
{
  float suma,ri;
  float Qxx,Qyy,Qzz,Qxy,Qxz,Qyz,Rxx,Ryy,Rzz,Rxy,Rxz,Ryz;
  
  ri=sqrt(xi*xi + yi*yi + zi*zi);
  
  Qxx = (3.0*xi*xi - ri*ri);
  Qyy = (3.0*yi*yi - ri*ri);
  Qzz = (3.0*zi*zi - ri*ri);

  Qxy = 3.0*xi*yi;
  Qxz = 3.0*xi*zi;
  Qyz = 3.0*yi*zi;
  
  Rxx = (3.0*x*x -R*R);
  Ryy = (3.0*y*y -R*R);
  Rzz = (3.0*z*z -R*R);

  Rxy = 3.0*x*y;
  Rxz = 3.0*x*z;
  Ryz = 3.0*y*z;

  suma= (Qxx*Rxx + Qxy*Rxy + Qxz*Rxz + Qxy*Rxy + Qyy*Ryy + Qyz*Ryz + Qxz*Rxz + Qyz*Ryz + Qzz*Rzz)/(6.0*R5);
    
  return suma;
    
}

float Dipole_summ(float xi,float yi,float zi,float x,float y,float z,float R3)
{
  
  float suma;

  suma= (xi*x + yi*y + zi*z)/R3;
  return suma;

}


float grav_soft_spline(float x, float h)
{
  
  float kern,u;

  u=x/h;

  if( (u >= 0) && (u < 0.5) )
    kern = (16.0/3.0)*u*u - (48.0/5.0)*pow(u,4) + (32.0/5.0)*pow(u,5) - (14.0/5.0); 
  else if( (u >= 0.5) && (u < 1.0) )
    kern= (1.0/(15.0*u)) + (32.0/3.0)*u*u - 16.0*pow(u,3) + (48.0/5.0)*pow(u,4) - (32.0/15.0)*pow(u,5) - (16.0/5.0);
  else if(u >= 1.0)
    kern=-1.0/u;

  return (-kern);
      
}

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

double ran1(long int *idum)
{
    int j;
    long int k;
    
    static long int iy=0;
    static long int iv[NTAB];
    double temp;
    
    if(*idum <= 0 || !iy){
	
    if(-(*idum) < 1) *idum=1;
    else *idum = -(*idum);
    
    for(j=NTAB+7;j>=0;j--)
    {
	k=(*idum)/IQ;
	*idum=IA*(*idum-k*IQ)-IR*k;
	if(*idum < 0) *idum += IM;
	if(j < NTAB) iv[j] = *idum;
    }
    iy=iv[0];
    }
    
    k=(*idum)/IQ;
    *idum=IA*(*idum-k*IQ)-IR*k;
    
    if(*idum < 0) *idum += IM;
    j=iy/NDIV;
    iy=iv[j];
    iv[j] = *idum;
    temp= AM*iy;
    if(temp > RNMX) return RNMX;
    else return temp;
    
}

#undef IA 
#undef IM 
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX
