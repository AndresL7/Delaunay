#include "GalCos_variables.h"

#define EXIT_ERROR printf("Error in parameter %s in parameter file\n",buf1); exit(0);
#define SKIP   fgets(buf,200,par_pf);


int GalCos_units(char *param_file)
{

  int aux_int;
  char buf[200],buf1[200],buf2[200];
  FILE *par_pf=NULL;
  
  
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
  
  if(task == 0) 
    {
      printf("\n====================================================\n"); fflush(stdout);
      printf(" >> UNITS AND PARAMETERS\n"); fflush(stdout);
    }
  
  fgets(buf,200,par_pf);
  if( sscanf(buf,"%s%s",buf1,buf2) < 2 ){
    EXIT_ERROR;
  }
  else{
    aux_int=atoi(buf2);
    if(task == 0) printf("%s %d\n",buf1,aux_int); fflush(stdout);
  }
  
  if(aux_int != 0)
    MINIMUM_MEMBERS=aux_int;
  
  fgets(buf,200,par_pf);
  if( sscanf(buf,"%s%s",buf1,buf2) < 2 ){
    EXIT_ERROR;
  }
  else{
    MINIMUM_NSUBSTRUCT=atoi(buf2);
    if(task == 0) printf("%s %d\n",buf1,MINIMUM_NSUBSTRUCT); fflush(stdout);
  }
  
  fgets(buf,200,par_pf);
  if( sscanf(buf,"%s%s",buf1,buf2) < 2 ){
    EXIT_ERROR;
  }
  else{
    b_Link=atof(buf2);
    if(task == 0) printf("%s %g\n",buf1,b_Link); fflush(stdout);
  }
  
  fgets(buf,200,par_pf);
  if( sscanf(buf,"%s%s",buf1,buf2) < 2 ){
    EXIT_ERROR;
  }
  else{
    NGB_MAX=atoi(buf2);
    if(task == 0) printf("%s %d\n",buf1,NGB_MAX); fflush(stdout);
  }
  
  fgets(buf,200,par_pf);
  if( sscanf(buf,"%s%s",buf1,buf2) < 2 ){
    EXIT_ERROR;
  }
  else{
    GRAV_SOFT=atof(buf2);
    if(task == 0) printf("%s %g\n",buf1,GRAV_SOFT); fflush(stdout);
  }
  SKIP;
  if(task == 0) printf("\n"); fflush(stdout);
  
  fgets(buf,200,par_pf);
  if( sscanf(buf,"%s%s",buf1,buf2) < 2 ){
    EXIT_ERROR;
  }
  else{
    FLAG_SUBFIND=atoi(buf2);
    if(task == 0) printf("%s %d\n",buf1,FLAG_SUBFIND); fflush(stdout);
  }
  SKIP;
  if(task == 0) printf("\n"); fflush(stdout);
  
  fgets(buf,200,par_pf);
  if( sscanf(buf,"%s%s",buf1,buf2) < 2 ){
    EXIT_ERROR;
  }
  else{
    OMEGABARYON=atof(buf2);
    if(task == 0) printf("%s %g\n",buf1,OMEGABARYON); fflush(stdout);
  }
  
  fgets(buf,200,par_pf);
  if( sscanf(buf,"%s%s",buf1,buf2) < 2 ){
    EXIT_ERROR;
  }
  else{
    OMEGA_MATTER=atof(buf2);
    if(task == 0) printf("%s %g\n",buf1,OMEGA_MATTER); fflush(stdout);
  }
  SKIP;
  if(task == 0) printf("\n"); fflush(stdout);

  fgets(buf,200,par_pf);
  if( sscanf(buf,"%s%s",buf1,buf2) < 2 ){
    EXIT_ERROR;
  }
  else{
    G_INTERNAL_UNITS=atof(buf2);
    if(task == 0) printf("%s %g\n",buf1,G_INTERNAL_UNITS); fflush(stdout);
  }

  fgets(buf,200,par_pf);
  if( sscanf(buf,"%s%s",buf1,buf2) < 2 ){
    EXIT_ERROR;
  }
  else{
    LENGHT_INTERNAL_UNITS=atof(buf2);
    if(task == 0) printf("%s %g\n",buf1,LENGHT_INTERNAL_UNITS); fflush(stdout);
  }

  fgets(buf,200,par_pf);
  if( sscanf(buf,"%s%s",buf1,buf2) < 2 ){
    EXIT_ERROR;
  }
  else{
    VELOCITY_INTERNAL_UNITS=atof(buf2);
    if(task == 0) printf("%s %g\n",buf1,VELOCITY_INTERNAL_UNITS); fflush(stdout);
  }

  fgets(buf,200,par_pf);
  if( sscanf(buf,"%s%s",buf1,buf2) < 2 ){
    EXIT_ERROR;
  }
  else{
    MASS_INTERNAL_UNITS=atof(buf2);
    if(task == 0) printf("%s %g\n",buf1,MASS_INTERNAL_UNITS); fflush(stdout);
  }

  fgets(buf,200,par_pf);
  if( sscanf(buf,"%s%s",buf1,buf2) < 2 ){
    EXIT_ERROR;
  }
  else{
    TIME_INTERNAL_UNITS=atof(buf2);
    if(task == 0) printf("%s %g\n",buf1,TIME_INTERNAL_UNITS); fflush(stdout);
  }

  fgets(buf,200,par_pf);
  if( sscanf(buf,"%s%s",buf1,buf2) < 2 ){
    EXIT_ERROR;
  }
  else{
    ENERGY_INTERNAL_UNITS=atof(buf2);
    if(task == 0) printf("%s %g\n",buf1,ENERGY_INTERNAL_UNITS); fflush(stdout);
  }
  
  fgets(buf,200,par_pf);
  if( sscanf(buf,"%s%s",buf1,buf2) < 2 ){
    EXIT_ERROR;
  }
  else{
    DENSITY_INTERNAL_UNITS=atof(buf2);
    if(task == 0) printf("%s %g\n",buf1,DENSITY_INTERNAL_UNITS); fflush(stdout);
  }

  fgets(buf,200,par_pf);
  if( sscanf(buf,"%s%s",buf1,buf2) < 2 ){
    EXIT_ERROR;
  }
  else{
    HUBBLE_INTERNAL_UNITS=atof(buf2);
    if(task == 0) printf("%s %g\n",buf1,HUBBLE_INTERNAL_UNITS); fflush(stdout);
  }
  SKIP;
  if(task == 0) printf("\n"); fflush(stdout);
    
  fclose(par_pf);

  if(task == 0) 
    printf("\n====================================================\n"); fflush(stdout);
  
  return (0);


}
