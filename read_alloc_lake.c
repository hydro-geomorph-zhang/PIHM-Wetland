#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

//#include "sundialstypes.h"
#include "pihm.h" 
#include "lake.h"


void read_alloc_lake(char *filename, LakeData LD, Model_Data DS)
{
    int i, j,k;
  	char *fn;
  	char tempchar[100];
  
  	FILE *lakegeom_file;	/* Pointer to .lakegeom file */
	FILE *lakebathy_file;	/* Pointer to .lakebathy file */
  	FILE *lakeatt_file;		/* Pointer to .lakeatt file */
  	FILE *lakesoil_file;	/* Pointer to .lakesoil file*/
    FILE *lakeriverinterest_file;	/* Pointer to .lakesoil file*/
 
  
 	printf("\nStart reading in Lake input files ... \n");
  

  	/*========== open *.lakegeom file ==========*/
//  	printf("\n  1) reading %s.lakegeom ... ", filename);
  	fn = (char *)malloc((strlen(filename)+1024)*sizeof(char));
  	strcpy(fn, filename);
  	lakegeom_file = fopen(strcat(fn, ".lakegeom"), "r");
	free(fn);
  	if(lakegeom_file == NULL)
  		{
    		printf("\n  Fatal Error: %s.mesh is in use or does not exist!\n", filename);
    		exit(1);
  		}
    
  	/* start reading lakegeom_file */ 
  	fscanf(lakegeom_file,"%d", &LD->NumLake);
	fscanf(lakegeom_file, "%s %s %s %s %s %s %s" ,tempchar,tempchar,tempchar,tempchar,tempchar,tempchar,tempchar);
  
  	LD->Lake = (Lake *)malloc((LD->NumLake)*sizeof(Lake));
    //LD->LakeBoundary =(int **)malloc((LD->NumLake)*sizeof(int));

  
  	/* read in elements information */ 
  	for (i=0; i<LD->NumLake; i++)
  	{
		fscanf(lakegeom_file, "%s %d %lf %lf %lf %lf %lf" ,tempchar, &(LD->Lake[i].NumBankEle), &(LD->Lake[i].SurfArea),&(LD->Lake[i].BankElevAdjust), &(LD->Lake[i].OutElevDecline), &(LD->Lake[i].Weir_length),&(LD->Lake[i].Orifice_height));
		LD->Lake[i].BankEle = (int *)malloc(LD->Lake[i].NumBankEle*sizeof(int));
        LD->Lake[i].StreamIndex = (int *)malloc(LD->Lake[i].NumBankEle*sizeof(int));        
		fscanf(lakegeom_file, "%s" ,tempchar);
		for(k=0;k<LD->Lake[i].NumBankEle;k++)
		{
			LD->Lake[i].StreamIndex[k] = -1;
			fscanf(lakegeom_file, "%d" ,&(LD->Lake[i].BankEle[k]));
			DS->Land_lake_index[LD->Lake[i].BankEle[k]-1] = i+1;
		}
	}
  	/* finish reading lakegeom_file */  
  	fclose(lakegeom_file);
	
	
	/*========== open *.lakebathy file ==========*/
//  	printf("\n  2) reading %s.lakebathy ... ", filename);
  	fn = (char *)malloc((strlen(filename)+1024)*sizeof(char));
  	strcpy(fn, filename);
  	lakebathy_file = fopen(strcat(fn, ".lakebathy"), "r");
	free(fn);
  	if(lakebathy_file == NULL)
  		{
    		printf("\n  Fatal Error: %s.mesh is in use or does not exist!\n", filename);
    		exit(1);
  		}
    
  	/* start reading lakebathy_file */ 
  	/* read in elements information */ 
  	for (i=0; i<LD->NumLake; i++)
  	{
		fscanf(lakebathy_file, "%s %lf %lf %lf" ,tempchar, &(LD->Lake[i].SurfIC), &(LD->Lake[i].BedElev), &(LD->Lake[i].BaseElev));
	}
 
  	/* finish reading lakebathy_file */  
  	fclose(lakebathy_file);
	
		/*========== open *.lakeatt_file  ==========*/
//  	printf("\n  3) reading %s.lakeatt ... ", filename);
  	fn = (char *)malloc((strlen(filename)+1024)*sizeof(char));
  	strcpy(fn, filename);
  	lakeatt_file = fopen(strcat(fn, ".lakeatt"), "r");
	free(fn);
  	if(lakeatt_file == NULL)
  		{
    		printf("\n  Fatal Error: %s.mesh is in use or does not exist!\n", filename);
    		exit(1);
  		}
    
  	/* start reading lakeatt_file */ 
  	/* read in elements information */ 
  	for (i=0; i<LD->NumLake; i++)
  	{
		fscanf(lakeatt_file, "%s %d  %d %d" ,tempchar, &(LD->Lake[i].Soil), &(LD->Lake[i].Meteo), &(LD->Lake[i].Macropore));
	}
 
  	/* finish reading lakeatt_file */  
  	fclose(lakeatt_file);
	
	
	/*========== open *.lakesoil_file file ==========*/  
//  	printf("\n  4) reading %s.lakesoil ... ", filename);
  	fn = (char *)malloc((strlen(filename)+1024)*sizeof(char));
  	strcpy(fn, filename);
  	lakesoil_file = fopen(strcat(fn, ".lakesoil"), "r");
	free(fn);
  	if(lakesoil_file == NULL)
  		{
    		printf("\n  Fatal Error: %s.soil is in use or does not exist!\n", filename);
    		exit(1);
  		}
  
  	/* start reading soil_file */  
  	fscanf(lakesoil_file, "%d", &LD->NumLakeSoil);
  	LD->LakeSoil = (Lake_soil *)malloc(LD->NumLakeSoil*sizeof(Lake_soil));
  
  	for (i=0; i<LD->NumLakeSoil; i++)
  	{
    	fscanf(lakesoil_file, "%d", &(LD->LakeSoil[i].index));
    	fscanf(lakesoil_file, "%lf %lf %lf %lf %lf",&(LD->LakeSoil[i].Kh),&(LD->LakeSoil[i].Kv),&(LD->LakeSoil[i].ThetaS),&(LD->LakeSoil[i].ThetaR),&(LD->LakeSoil[i].macKh));
  	} 
 
  	fclose(lakesoil_file);
    
    
    /*========== open *.lakeriverInterest file ==========*/
    //  	printf("\n  1) reading %s.lakegeom ... ", filename);
    fn = (char *)malloc((strlen(filename)+1024)*sizeof(char));
    strcpy(fn, filename);
    lakeriverinterest_file = fopen(strcat(fn, ".riverInterest"), "r");
    free(fn);
    if(lakeriverinterest_file == NULL)
    {
        printf("\n  Fatal Error: %s.riverInterest is in use or does not exist!\n", filename);
        printf("\n  The river flow will sum up together!\n");
        LD->NumRiverInterest = 0;
    }
    else
    {
        /* start reading lakegeom_file */
        fscanf(lakeriverinterest_file,"%d", &(LD->NumRiverInterest));
    
        LD->RiverInterestInformation = (int ***)malloc((LD->NumLake)*sizeof(int**));
        LD->RiverInterest = (int *)malloc((LD->NumRiverInterest)*sizeof(int));
        
        //LD->LakeBoundary =(int **)malloc((LD->NumLake)*sizeof(int));
        for (i=0; i<LD->NumLake; i++)
        {
            LD->RiverInterestInformation[i] = (int **)malloc((LD->NumRiverInterest)*sizeof(int*));
            for (j=0; j<LD->NumRiverInterest; j++)
            {
                LD->RiverInterestInformation[i][j] = (int *)malloc((2)*sizeof(int));
            }
        }
        
        for (i=0; i<LD->NumRiverInterest; i++)
        {
            fscanf(lakeriverinterest_file, "%d %s" ,&(LD->RiverInterest[i]),tempchar);
        }
        fclose(lakeriverinterest_file);
    }
}

	
	
	
	
	
	
  