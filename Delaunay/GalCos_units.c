#include "GalCos_variables.h"

#define EXIT_ERROR printf("Error in parameter %s in parameter file\n",buf1); exit(0);
#define SKIP   fgets(buf,200,par_pf);


int GalCos_units(char *param_file)
{

  int aux_int,i;
  char buf[200],buf1[200],buf2[200];
  FILE *par_pf;
  
  
  if(NULL==(par_pf=fopen(param_file,"r")))
    {
      printf("Parameterfile %s not found...\n",param_file);
      exit(0);
    }
  
  /*
    Reading internal units of the simulation. Each of them are the
    equivalence in cgs units, example, LENGHT_INTERNAL_UNITS is the
    internal unit of lenght in cm.
  */
  
  printf("\n====================================================\n");
  printf(" >> UNITS AND PARAMETERS\n");

  /* Basica data */

  fgets(buf,200,par_pf);
  if( sscanf(buf,"%s%s",buf1,buf2) < 2 ){
    EXIT_ERROR;
  }
  else{
    aux_int=atoi(buf2);
    printf("%s %d\n",buf1,aux_int);
  }
  
  if(aux_int != 0)
    MINIMUM_MEMBERS=aux_int;

  fgets(buf,200,par_pf);
  if( sscanf(buf,"%s%s",buf1,buf2) < 2 ){
    EXIT_ERROR;
  }
  else{
    MINIMUM_NSUBSTRUCT=atoi(buf2);
    printf("%s %d\n",buf1,MINIMUM_NSUBSTRUCT);
  }
    
  fgets(buf,200,par_pf);
  if( sscanf(buf,"%s%s",buf1,buf2) < 2 ){
    EXIT_ERROR;
  }
  else{
    b_Link=atof(buf2);
    printf("%s %g\n",buf1,b_Link);
  }
  
  fgets(buf,200,par_pf);
  if( sscanf(buf,"%s%s",buf1,buf2) < 2 ){
    EXIT_ERROR;
  }
  else{
    NGB_MAX=atoi(buf2);
    printf("%s %d\n",buf1,NGB_MAX);
  }
  
  fgets(buf,200,par_pf);
  if( sscanf(buf,"%s%s",buf1,buf2) < 2 ){
    EXIT_ERROR;
  }
  else{
    GRAV_SOFT=atof(buf2);
    printf("%s %g\n",buf1,GRAV_SOFT);
  }
  
  fgets(buf,200,par_pf);
  if( sscanf(buf,"%s%s",buf1,buf2) < 2 ){
    EXIT_ERROR;
  }
  else{
    FLAG_SUBFIND=atoi(buf2);
    printf("%s %d\n",buf1,FLAG_SUBFIND);
  }
  SKIP;
  printf("\n");
  
  fgets(buf,200,par_pf);
  if( sscanf(buf,"%s%s",buf1,buf2) < 2 ){
    EXIT_ERROR;
  }
  else{
    NBINS=atoi(buf2);
    printf("%s %d\n",buf1,NBINS);
  }
  
  Mbins = (int *) malloc((size_t) (NBINS+1)*sizeof(int));
  
  for(i=0; i<NBINS+1; i++)
    {
      fgets(buf,200,par_pf);
      if( sscanf(buf,"%s%s",buf1,buf2) < 2 ){
	EXIT_ERROR;
      }
      else{
	Mbins[i]=atoi(buf2);
	printf("%s %d\n",buf1,Mbins[i]);
      }
      
    }
  SKIP;
  printf("\n");
  
  fgets(buf,200,par_pf);
  if( sscanf(buf,"%s%s",buf1,buf2) < 2 ){
    EXIT_ERROR;
  }
  else{
    RMIN=atof(buf2);
    printf("%s %g\n",buf1,RMIN);
  }
  
  fgets(buf,200,par_pf);
  if( sscanf(buf,"%s%s",buf1,buf2) < 2 ){
    EXIT_ERROR;
  }
  else{
    RMAX=atof(buf2);
    printf("%s %g\n",buf1,RMAX);
  }

  fgets(buf,200,par_pf);
  if( sscanf(buf,"%s%s",buf1,buf2) < 2 ){
    EXIT_ERROR;
  }
  else{
    NITER=atoi(buf2);
    printf("%s %d\n",buf1,NITER);
  }
  SKIP;
  printf("\n");
  

  fgets(buf,200,par_pf);
  if( sscanf(buf,"%s%s",buf1,buf2) < 2 ){
    EXIT_ERROR;
  }
  else{
    OMEGABARYON=atof(buf2);
    printf("%s %g\n",buf1,OMEGABARYON);
  }
  
  fgets(buf,200,par_pf);
  if( sscanf(buf,"%s%s",buf1,buf2) < 2 ){
    EXIT_ERROR;
  }
  else{
    OMEGA_MATTER=atof(buf2);
    printf("%s %g\n",buf1,OMEGA_MATTER);
  }
  SKIP;
  printf("\n");
  
  fgets(buf,200,par_pf);
  if( sscanf(buf,"%s%s",buf1,buf2) < 2 ){
    EXIT_ERROR;
  }
  else{
    G_INTERNAL_UNITS=atof(buf2);
    printf("%s %g\n",buf1,G_INTERNAL_UNITS);
  }

  fgets(buf,200,par_pf);
  if( sscanf(buf,"%s%s",buf1,buf2) < 2 ){
    EXIT_ERROR;
  }
  else{
    LENGHT_INTERNAL_UNITS=atof(buf2);
    printf("%s %g\n",buf1,LENGHT_INTERNAL_UNITS);
  }

  fgets(buf,200,par_pf);
  if( sscanf(buf,"%s%s",buf1,buf2) < 2 ){
    EXIT_ERROR;
  }
  else{
    VELOCITY_INTERNAL_UNITS=atof(buf2);
    printf("%s %g\n",buf1,VELOCITY_INTERNAL_UNITS);
  }
  
  fgets(buf,200,par_pf);
  if( sscanf(buf,"%s%s",buf1,buf2) < 2 ){
    EXIT_ERROR;
  }
  else{
    MASS_INTERNAL_UNITS=atof(buf2);
    printf("%s %g\n",buf1,MASS_INTERNAL_UNITS);
  }
  
  fgets(buf,200,par_pf);
  if( sscanf(buf,"%s%s",buf1,buf2) < 2 ){
    EXIT_ERROR;
  }
  else{
    TIME_INTERNAL_UNITS=atof(buf2);
    printf("%s %g\n",buf1,TIME_INTERNAL_UNITS);
  }

  fgets(buf,200,par_pf);
  if( sscanf(buf,"%s%s",buf1,buf2) < 2 ){
    EXIT_ERROR;
  }
  else{
    ENERGY_INTERNAL_UNITS=atof(buf2);
    printf("%s %g\n",buf1,ENERGY_INTERNAL_UNITS);
  }
  
  fgets(buf,200,par_pf);
  if( sscanf(buf,"%s%s",buf1,buf2) < 2 ){
    EXIT_ERROR;
  }
  else{
    DENSITY_INTERNAL_UNITS=atof(buf2);
    printf("%s %g\n",buf1,DENSITY_INTERNAL_UNITS);
  }
  
  fgets(buf,200,par_pf);
  if( sscanf(buf,"%s%s",buf1,buf2) < 2 ){
    EXIT_ERROR;
  }
  else{
    HUBBLE_INTERNAL_UNITS=atof(buf2);
    printf("%s %g\n",buf1,HUBBLE_INTERNAL_UNITS);
  }
  SKIP;
  printf("\n");
  
  fclose(par_pf);

  printf("\n====================================================\n");
  
  return (SUCCES);


}
