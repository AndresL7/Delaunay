#include "GalCos_variables.h"
#include "GalCos_cola.h"

void load_cola(struct COLA *cola)
{
  
  cola->inputs=(int *) malloc((size_t) sizeof(int));
  
  if(cola->inputs == NULL){
    printf("An error in memory allocation 1\n");
    exit(FAILURE);
  }
  
  cola->cola_Nmembers=1;
  cola->head=0; // index number for the first element
  cola->tail=0; // index number for the last element
  
}


void reinitialize_cola(struct COLA *cola)
{
  
  if(cola->tail != 0) 
    free(cola->inputs);
  
  cola->inputs=(int *) malloc((size_t) sizeof(int));
  
  if(cola->inputs == NULL){
    printf("NO memory available cola:30\n");
    exit(0);
  }
  
  cola->cola_Nmembers=1;
  cola->head=0;
  cola->tail=0;
  
}


void encola(struct COLA *cola,int entrance)
{  
  int *temp=NULL;
  
  temp=(int *) realloc(cola->inputs,(size_t) (cola->cola_Nmembers+1)*sizeof(int));
  cola->inputs=temp;
  cola->cola_Nmembers=cola->cola_Nmembers+1;
  
  if(cola->tail == 0)
    {
      cola->inputs[0]=entrance;
      cola->tail=cola->tail+1;
    }
  else
    { 
      cola->inputs[cola->tail]=entrance;  
      cola->tail=cola->tail+1;
    }
  
}


int desencola(struct COLA *cola)
{
  
  int i,head_cola,tail_cola;
  int *temp=NULL;
  
  head_cola=cola->inputs[0];
  tail_cola=cola->tail;
    
  for(i=1; i<tail_cola; i++)
    {
      cola->inputs[i-1]=cola->inputs[i];      
    }
  
  cola->tail=cola->tail-1;
  cola->cola_Nmembers=cola->cola_Nmembers-1;
  
  temp=(int *) realloc(cola->inputs,(size_t) (cola->cola_Nmembers-1)*sizeof(int));
  
  cola->inputs=temp;
  
  return head_cola;
  
}

/*
  COLA FLOAT PARA NGBs
 */


void load_cola_float(struct COLA_float *colaf)
{
  
  colaf->inputs=(float *) malloc((size_t) sizeof(float));
  
  if(colaf->inputs == NULL){
    printf("An error in memory allocation 1\n");
    exit(FAILURE);
  }
  
  colaf->cola_Nmembers=1;
  colaf->head=0; // index number for the first element
  colaf->tail=0; // index number for the last element
  
}


void reinitialize_cola_float(struct COLA_float *colaf)
{
  
  if(colaf->tail != 0) 
    free(colaf->inputs);
  
  colaf->inputs=(float *) malloc((size_t) sizeof(float));
  
  if(colaf->inputs == NULL){
    printf("NO memory available cola:30\n");
    exit(0);
  }
  
  colaf->cola_Nmembers=1;
  colaf->head=0;
  colaf->tail=0;
  
}


void encola_float(struct COLA_float *colaf,float entrance)
{  
  
  float *temp=NULL;
  
  temp = (float *) realloc(colaf->inputs,(size_t) (colaf->cola_Nmembers+1)*sizeof(float));
  
  colaf->inputs = temp;
  colaf->cola_Nmembers = colaf->cola_Nmembers+1;
  
  if(colaf->tail == 0)
    {
      colaf->inputs[0] = entrance;
      colaf->tail = colaf->tail+1;
    }
  else
    { 
      colaf->inputs[colaf->tail] = entrance;  
      colaf->tail = colaf->tail+1;
    }
  
}


float desencola_float(struct COLA_float *colaf)
{
    
  int i;
  float head_cola,tail_cola;
  float *temp=NULL;
  
  head_cola = colaf->inputs[0];
  tail_cola = colaf->tail;
  
  for(i=1; i<tail_cola; i++)
    {
      colaf->inputs[i-1] = colaf->inputs[i];      
    }
  
  colaf->tail = colaf->tail-1;
  colaf->cola_Nmembers = colaf->cola_Nmembers-1;
  
  temp = (float *) realloc(colaf->inputs,(size_t) (colaf->cola_Nmembers-1)*sizeof(float));
  
  colaf->inputs=temp;
  
  return head_cola;
  
}

