int GalCos_write_cluster(struct halo Haloinp)
{
  
  int k,ipart,ID_cluster;
  char buf[100];  
  
  FILE *pf_out;
  
  ID_cluster=Haloinp.IDcluster;
  
  sprintf(buf,"%s%d%s","clusterH_",ID_cluster,".glc");
  
  pf_out=fopen(buf,"w");
  
  for(k=0; k<Haloinp.Nmembers; k++)
    {
      ipart=Haloinp.Halo_particles[k];
      
      fprintf(pf_out,"%d  %d  %f  %f  %f  %f  %f  %f  %f %f %f\n",Particle[ipart].id,Particle[ipart].Cluster_ID,
	      Particle[ipart].mass,Particle[ipart].pos[0]+Haloinp.Bpos[0],Particle[ipart].pos[1]+Haloinp.Bpos[1],Particle[ipart].pos[2]+Haloinp.Bpos[2],
	      Particle[ipart].vel[0],Particle[ipart].vel[1],Particle[ipart].vel[2],Particle[ipart].EP,
	      Particle[ipart].Volume);
      
    }
  
  fclose(pf_out);
  
  
  /////
  
  sprintf(buf,"%s%d%s","clusterD_",ID_cluster,".glc");
  
  pf_out=fopen(buf,"w");
  
  for(k=0; k<Haloinp.NDomain_particles; k++)
    {
      ipart=Haloinp.Domain_particles[k];
      
      fprintf(pf_out,"%d  %d  %f  %f  %f  %f  %f  %f  %f %f %f\n",Particle[ipart].id,Particle[ipart].Cluster_ID,
	      Particle[ipart].mass,Particle[ipart].pos[0]+Haloinp.Bpos[0],Particle[ipart].pos[1]+Haloinp.Bpos[1],Particle[ipart].pos[2]+Haloinp.Bpos[2],
	      Particle[ipart].vel[0],Particle[ipart].vel[1],Particle[ipart].vel[2],Particle[ipart].EP,
	      Particle[ipart].Volume);
      
    }
  
  fclose(pf_out);
  
  return 0;

}

int GalCos_write_halo_Bproperties(char *infile,int NCLUSTERS, int FINAL_NHALOS)
{

  int k,i;
  char buf[100];  
  
  FILE *pf_out;

  sprintf(buf,"%s%s",infile,".Bprop");
  
  pf_out=fopen(buf,"w");
  
  
  if(NCLUSTERS > 0)
    {
      fprintf(pf_out,"%d\n",FINAL_NHALOS);
      fprintf(pf_out,"%g\n",REDSHIFT);
    }
  else
    {
      fprintf(pf_out,"%d\n",0);
      fprintf(pf_out,"%g\n",REDSHIFT);
    }
  
  if(NCLUSTERS > 0)
    {
      
      for(i=0; i<NCLUSTERS; i++)
	{
	  
	  fprintf(pf_out,"%d ",Halos[i].ID_CenterHalo);
	  fprintf(pf_out,"%g ",Halos[i].mass);
	  
	  fprintf(pf_out,"%g ",Halos[i].pos[0]); 
	  fprintf(pf_out,"%g ",Halos[i].pos[1]);    
	  fprintf(pf_out,"%g ",Halos[i].pos[2]);
	      
	  fprintf(pf_out,"%g ",Halos[i].vel[0]);    
	  fprintf(pf_out,"%g ",Halos[i].vel[1]);    
	  fprintf(pf_out,"%g ",Halos[i].vel[2]);
	  
	  fprintf(pf_out,"%d ",Halos[i].Nmembers);
	  fprintf(pf_out,"%d ",Halos[i].IDcluster);
	  fprintf(pf_out,"%d ",Halos[i].NumberOfSubhalos);
	  fprintf(pf_out,"%g ",Halos[i].Mass_fraction);
	  fprintf(pf_out,"%g\n",Halos[i].Hmr);
	  
	  for(k=0; k<Halos[i].Nmembers; k++)
	    {
	      fprintf(pf_out,"%d ",Halos[i].Halo_particles[k]);
	    }
	  fprintf(pf_out,"\n");
	  
	  fprintf(pf_out,"%g %g %g %g %g %g %g\n",Halos[i].Rvir,Halos[i].Mvir,
		  Halos[i].Vvir,Halos[i].Tvir,Halos[i].TotalEnergy,Halos[i].TotalAngularMom,Halos[i].Lambda_spin);
	  
	  fprintf(pf_out,"\n\n");
	  
	}
    }
  else
    {
      fprintf(pf_out,"%g ",0.); //Halos[i].IDcenterpos
      fprintf(pf_out,"%g ",0.); //Halos[i].mass
      
      fprintf(pf_out,"%g ",0.); //Halos[i].pos[0]
      fprintf(pf_out,"%g ",0.); //Halos[i].pos[1]   
      fprintf(pf_out,"%g ",0.); //Halos[i].pos[2]
      
      fprintf(pf_out,"%g ",0.); //Halos[i].vel[0]    
      fprintf(pf_out,"%g ",0.); //Halos[i].vel[1]   
      fprintf(pf_out,"%g ",0.); //Halos[i].vel[2]
      
      fprintf(pf_out,"%d ",0); //Halos[i].Nmembers
      fprintf(pf_out,"%d ",0); //Halos[i].IDcluster
      fprintf(pf_out,"%d ",0); //Halos[i].NumberOfSubhalos
      fprintf(pf_out,"%g ",0.0); //Halos[i].Mass_fraction
      fprintf(pf_out,"%g\n",0.0); //Halos[i].Hmr
      
      fprintf(pf_out,"%g %g %g %g %g %g %g\n",0.,0.,0.,0.,0.,0.,0.);
      
      fprintf(pf_out,"\n");
      
    }

      fclose(pf_out);      
  
      
  return 0;
  
}


int GalCos_write_histograms(int NCLUSTERS, char *infile)
{

  int i;
  char buf[100];  
  
  FILE *pf_out;
  
  ////////////////////////////////////////////////////////
  /// TOTAL MASS DISTRIBUTION

  sprintf(buf,"%s%s%s","MasesTotal_",infile,".glc");
  pf_out=fopen(buf,"w");
  
  for(i=0; i<NCLUSTERS; i++)
    {
      fprintf(pf_out,"%g\n",Halos[i].mass);
    }
  
  fclose(pf_out);
  
  ////////////////////////////////////////////////////////
  /// VIRIAL MASS DISTRIBUTION
    
  sprintf(buf,"%s%s%s","MasesVirial_",infile,".glc");
  pf_out=fopen(buf,"w");
  
  for(i=0; i<NCLUSTERS; i++)
    {
      fprintf(pf_out,"%g\n",Halos[i].Mvir);
    }
  
  fclose(pf_out);
  
  ////////////////////////////////////////////////////////
  /// VIRIAL RADIUS

  sprintf(buf,"%s%s%s","VirialRadius_",infile,".glc");
  pf_out=fopen(buf,"w");
  
  for(i=0; i<NCLUSTERS; i++)
    {
      fprintf(pf_out,"%g\n",Halos[i].Rvir);
    }
  
  fclose(pf_out);
  
  ////////////////////////////////////////////////////////
  /// LAMBDA PARAMETER

  sprintf(buf,"%s%s%s","LambdaParam_",infile,".glc");
  pf_out=fopen(buf,"w");
  
  for(i=0; i<NCLUSTERS; i++)
    {
      fprintf(pf_out,"%g\n",Halos[i].Lambda_spin);
    }
  
  fclose(pf_out);

  ////////////////////////////////////////////////////////
  /// TOTAL MEAN DENSITY

  sprintf(buf,"%s%s%s","TotalMeanDens_",infile,".glc");
  pf_out=fopen(buf,"w");
  
  for(i=0; i<NCLUSTERS; i++)
    {
      fprintf(pf_out,"%g\n",3.0*Halos[i].mass/4.0*M_PI*pow(Halos[i].Radius,3));
    }
  
  fclose(pf_out);
  
  
  ////////////////////////////////////////////////////////
  /// VIRIAL MEAN DENSITY
  
  sprintf(buf,"%s%s%s","VirialMeanDens_",infile,".glc");
  pf_out=fopen(buf,"w");
  
  for(i=0; i<NCLUSTERS; i++)
    {
      fprintf(pf_out,"%g\n",3.0*Halos[i].Mvir/4.0*M_PI*pow(Halos[i].Rvir,3));
    }
  
  fclose(pf_out);
  
  ////////////////////////////////////////////////////////
  /// TOTAL MEAN DENSITY
  
  sprintf(buf,"%s%s%s","HmrMeanDens_",infile,".glc");
  pf_out=fopen(buf,"w");
  
  for(i=0; i<NCLUSTERS; i++)
    {
      fprintf(pf_out,"%g\n",3.0*Halos[i].mass/8.0*M_PI*pow(Halos[i].Hmr,3));
    }
  
  fclose(pf_out);

  ////////////////////////////////////////////////////////
  /// VIRIAL FACTOR 
  
  sprintf(buf,"%s%s%s","VirialDist_",infile,".glc");
  pf_out=fopen(buf,"w");
  
  for(i=0; i<NCLUSTERS; i++)
    {
      fprintf(pf_out,"%g\n",2.0*Halos[i].TotalEk+Halos[i].TotalEP);
    }
  
  fclose(pf_out);

  return 0;
}  


void number_of_subhalos(int IDhalo,int *Total_subhalos)
{
  
  int i,subhaloID;
  
  for(i=0; i<Halos[IDhalo].NumberOfSubhalos; i++)
    {
      subhaloID=Halos[IDhalo].ID_subhalos[i];
      *Total_subhalos = *Total_subhalos + Halos[subhaloID].NumberOfSubhalos;
      number_of_subhalos(subhaloID,Total_subhalos);
    }
    
}

  
int GalCos_writeAll(struct halo Halo,char *infile)
{
  
  int Total_subhalos;
  char buf[100];
  FILE *pf_out;

  sprintf(buf,"%s%s%s","All_",infile,".glc");
  pf_out=fopen(buf,"a");
  
  Total_subhalos=0;
  number_of_subhalos(Halo.IDcluster,&Total_subhalos);
  
  Total_subhalos=Total_subhalos+Halo.NumberOfSubhalos;

  fprintf(pf_out,"%d\t %d\t %g\t %g\t %g\t %g\t %g\t %d\t %g\t %g\t %g\t %g\t %g\t %g %d\n",Halo.IDcluster,Halo.Nmembers,Halo.mass,Halo.Mass_fraction,
	  Halo.Hmr,Halo.TotalEk,Halo.TotalEP,Halo.NumberOfSubhalos,Halo.Radius,Halo.Rvir,Halo.Mvir,Halo.Vvir,
	  Halo.Tvir,Halo.Lambda_spin,Total_subhalos);
  
  fclose(pf_out);

  //printf("El numero total de subhalos es %d",Total_subhalos);

  return 0;
 
}

int GalCos_writeSubhaloDists(int NCLUSTERS,char *infile)
{
  
  return 0;
  
}

