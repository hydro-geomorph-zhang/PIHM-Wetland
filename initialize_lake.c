#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


#include "pihm.h"
#include "lake.h"

void initialize_lake(char *filename, LakeData LD, Model_Data DS, Control_Data * CS, FILE *Lakeboundarynode,FILE *LakeReference)
{
	int           i, j, k, t, q,tmpBool, BoolR = 0,inabr, count1=0,count2=0, EleID;
	double        a_x, a_y, b_x, b_y, c_x, distX, distY;
	double        a_zmin, a_zmax, b_zmin, b_zmax, c_zmin, c_zmax;
	double        tempvalue1, tempvalue2, tempvalue3,tempvalue4,tempvalue5;
	FILE           *init_file;
	char           *fn;
	int			  LakeNodeIndex[2];
	int			  count=0;
	int 		  StartNode, EndNode, FindIndex;
	int			  find = 0;
	int 		  edge = 0;


	printf("\nInitializing Lake data structure ... ");

	/* allocate memory storage to flux terms */
  	LD->FluxSurf = (double **)malloc(LD->NumLake*sizeof(double*));
  	LD->FluxStream = (double **)malloc(LD->NumLake*sizeof(double*));
  	LD->FluxSub = (double **)malloc(LD->NumLake*sizeof(double*));
	LD->Infil = (double *)malloc(LD->NumLake*sizeof(double));
  	LD->ET = (double *)malloc(LD->NumLake*sizeof(double));
	LD->LakeSurfDepth = (double *)malloc(LD->NumLake*sizeof(double));
	LD->LakeGW = (double *)malloc(LD->NumLake*sizeof(double));
	LD->LakeDummyY = (double *)malloc(LD->NumLake*sizeof(double));
	LD->BankSurf = (double **)malloc(LD->NumLake*sizeof(double*));
	LD->BankGW = (double **)malloc(LD->NumLake*sizeof(double*));
	for(i=0;i<LD->NumLake;i++)
	{
		LD->Lake[i].NumLakeNode = 0;
		LD->Lake[i].Kh = CS->Cal.LakeKh * LD->LakeSoil[(LD->Lake[i].Soil - 1)].Kh;
		LD->Lake[i].Kv = CS->Cal.LakeKv * LD->LakeSoil[(LD->Lake[i].Soil - 1)].Kv;
		LD->Lake[i].ThetaS = CS->Cal.LakeThetaS * LD->LakeSoil[(LD->Lake[i].Soil - 1)].ThetaS;
		LD->Lake[i].ThetaR = CS->Cal.LakeThetaR * LD->LakeSoil[(LD->Lake[i].Soil - 1)].ThetaR;
		LD->Lake[i].macKh = CS->Cal.LakemacKh * LD->LakeSoil[(LD->Lake[i].Soil - 1)].macKh;
		
		for (k=0; k<3; k++)
        {
			if (DS->Ele[LD->Lake[i].BankEle[0]-1].nabr[k]==0)
            {
				edge+=1;
            }
		}
		
		for (k=0; k<3; k++)
        {
            if (edge==1)
			{
				if (DS->Ele[LD->Lake[i].BankEle[0]-1].nabr[k]!=0)
				{
					LakeNodeIndex[count++] = k;
				}
			}
			if (edge==2)
			{
				if (DS->Ele[LD->Lake[i].BankEle[0]-1].nabr[k]==0)
				{
					LakeNodeIndex[count++] = k;
				}
			}
        }
		
		for (k=0; k<3; k++)
        {
			if (edge==2)
			{
				if (DS->Ele[LD->Lake[i].BankEle[0]-1].nabr[k]!=0)
				{
					fprintf(Lakeboundarynode, "%lf\t%lf\n", DS->Node[DS->Ele[LD->Lake[i].BankEle[0]-1].node[k] - 1].x,DS->Node[DS->Ele[LD->Lake[i].BankEle[0]-1].node[k] - 1].y);
					LD->Lake[i].NumLakeNode = LD->Lake[i].NumLakeNode + 1;
				}
			}
        }
		
		StartNode = DS->Ele[LD->Lake[i].BankEle[0]-1].node[LakeNodeIndex[0]]; 
		EndNode = DS->Ele[LD->Lake[i].BankEle[0]-1].node[LakeNodeIndex[1]];
		
		fprintf(Lakeboundarynode, "%lf\t%lf\n", DS->Node[StartNode - 1].x,DS->Node[StartNode - 1].y);
		LD->Lake[i].NumLakeNode = LD->Lake[i].NumLakeNode + 1;
		count =1;
		edge = 0;
		printf("\n lake %d StartNode: %d  EndNode: %d\n",i, StartNode, EndNode);
		while(StartNode != EndNode)
		{
			for (j=count;j<LD->Lake[i].NumBankEle;j++)
			{
				for (k=0; k<3; k++)
				{
					if (DS->Ele[LD->Lake[i].BankEle[j]-1].nabr[k]==0)
					{
						edge+=1;
					}
				}
				
				for (k=0; k<3; k++)
				{
					if (edge==1)
					{
						if (DS->Ele[LD->Lake[i].BankEle[j]-1].nabr[k]!=0)
						{
							if (DS->Ele[LD->Lake[i].BankEle[j]-1].node[k]==StartNode)
							{
								FindIndex = LD->Lake[i].BankEle[j];
								
								if (k==0)
								{
									if (DS->Ele[LD->Lake[i].BankEle[j]-1].nabr[k+1]!=0) 
									{
										StartNode = DS->Ele[LD->Lake[i].BankEle[j]-1].node[k+1];
										fprintf(Lakeboundarynode, "%lf\t%lf\n", DS->Node[StartNode - 1].x,DS->Node[StartNode - 1].y);
										LD->Lake[i].NumLakeNode = LD->Lake[i].NumLakeNode + 1;
									}
										
									if (DS->Ele[LD->Lake[i].BankEle[j]-1].nabr[k+2]!=0)
									{
										StartNode = DS->Ele[LD->Lake[i].BankEle[j]-1].node[k+2];
										fprintf(Lakeboundarynode, "%lf\t%lf\n", DS->Node[StartNode - 1].x,DS->Node[StartNode - 1].y);
										LD->Lake[i].NumLakeNode = LD->Lake[i].NumLakeNode + 1;
									}
										
								}
								if (k==1)
								{
									if (DS->Ele[LD->Lake[i].BankEle[j]-1].nabr[k-1]!=0)
									{
										StartNode = DS->Ele[LD->Lake[i].BankEle[j]-1].node[k-1];
										fprintf(Lakeboundarynode, "%lf\t%lf\n", DS->Node[StartNode - 1].x,DS->Node[StartNode - 1].y);
										LD->Lake[i].NumLakeNode = LD->Lake[i].NumLakeNode + 1;
									}
										
									if (DS->Ele[LD->Lake[i].BankEle[j]-1].nabr[k+1]!=0)
									{
										StartNode = DS->Ele[LD->Lake[i].BankEle[j]-1].node[k+1];
										fprintf(Lakeboundarynode, "%lf\t%lf\n", DS->Node[StartNode - 1].x,DS->Node[StartNode - 1].y);
										LD->Lake[i].NumLakeNode = LD->Lake[i].NumLakeNode + 1;
									}
										
								}
								if (k==2)
								{
									if (DS->Ele[LD->Lake[i].BankEle[j]-1].nabr[k-2]!=0)
									{
										StartNode = DS->Ele[LD->Lake[i].BankEle[j]-1].node[k-2];
										fprintf(Lakeboundarynode, "%lf\t%lf\n", DS->Node[StartNode - 1].x,DS->Node[StartNode - 1].y);
										LD->Lake[i].NumLakeNode = LD->Lake[i].NumLakeNode + 1;
									}
										
									if (DS->Ele[LD->Lake[i].BankEle[j]-1].nabr[k-1]!=0)
									{
										StartNode = DS->Ele[LD->Lake[i].BankEle[j]-1].node[k-1];
										fprintf(Lakeboundarynode, "%lf\t%lf\n", DS->Node[StartNode - 1].x,DS->Node[StartNode - 1].y);
										LD->Lake[i].NumLakeNode = LD->Lake[i].NumLakeNode + 1;
									}
										
								}
								LD->Lake[i].BankEle[j] = LD->Lake[i].BankEle[count];
								LD->Lake[i].BankEle[count] = FindIndex;
								count++;
								find = 1;
								printf("	--%d",StartNode);
								break;
							}						
						}
					}
					if (edge==2)
					{
						if (DS->Ele[LD->Lake[i].BankEle[j]-1].nabr[k]==0)
						{
							if (DS->Ele[LD->Lake[i].BankEle[j]-1].node[k]==StartNode)
							{
								FindIndex = LD->Lake[i].BankEle[j];
								
								if (k==0)
								{
									if (DS->Ele[LD->Lake[i].BankEle[j]-1].nabr[k+1]==0)
									{
										fprintf(Lakeboundarynode, "%lf\t%lf\n", DS->Node[DS->Ele[LD->Lake[i].BankEle[j]-1].node[k+2] - 1].x,DS->Node[DS->Ele[LD->Lake[i].BankEle[j]-1].node[k+2] - 1].y);
										LD->Lake[i].NumLakeNode = LD->Lake[i].NumLakeNode + 1;
										StartNode = DS->Ele[LD->Lake[i].BankEle[j]-1].node[k+1];
										fprintf(Lakeboundarynode, "%lf\t%lf\n", DS->Node[StartNode - 1].x,DS->Node[StartNode - 1].y);
										LD->Lake[i].NumLakeNode = LD->Lake[i].NumLakeNode + 1;
									}
										
									if (DS->Ele[LD->Lake[i].BankEle[j]-1].nabr[k+2]==0)
									{
										fprintf(Lakeboundarynode, "%lf\t%lf\n", DS->Node[DS->Ele[LD->Lake[i].BankEle[j]-1].node[k+1] - 1].x,DS->Node[DS->Ele[LD->Lake[i].BankEle[j]-1].node[k+1] - 1].y);
										LD->Lake[i].NumLakeNode = LD->Lake[i].NumLakeNode + 1;
										StartNode = DS->Ele[LD->Lake[i].BankEle[j]-1].node[k+2];
										fprintf(Lakeboundarynode, "%lf\t%lf\n", DS->Node[StartNode - 1].x,DS->Node[StartNode - 1].y);
										LD->Lake[i].NumLakeNode = LD->Lake[i].NumLakeNode + 1;
									}										
								}
								if (k==1)
								{
									if (DS->Ele[LD->Lake[i].BankEle[j]-1].nabr[k-1]==0) 
									{
										fprintf(Lakeboundarynode, "%lf\t%lf\n", DS->Node[DS->Ele[LD->Lake[i].BankEle[j]-1].node[k+1] - 1].x,DS->Node[DS->Ele[LD->Lake[i].BankEle[j]-1].node[k+1] - 1].y);
										LD->Lake[i].NumLakeNode = LD->Lake[i].NumLakeNode + 1;
										StartNode = DS->Ele[LD->Lake[i].BankEle[j]-1].node[k-1];
										fprintf(Lakeboundarynode, "%lf\t%lf\n", DS->Node[StartNode - 1].x,DS->Node[StartNode - 1].y);
										LD->Lake[i].NumLakeNode = LD->Lake[i].NumLakeNode + 1;
									}
										
									if (DS->Ele[LD->Lake[i].BankEle[j]-1].nabr[k+1]==0)
									{
										fprintf(Lakeboundarynode, "%lf\t%lf\n", DS->Node[DS->Ele[LD->Lake[i].BankEle[j]-1].node[k-1] - 1].x,DS->Node[DS->Ele[LD->Lake[i].BankEle[j]-1].node[k-1] - 1].y);
										LD->Lake[i].NumLakeNode = LD->Lake[i].NumLakeNode + 1;
										StartNode = DS->Ele[LD->Lake[i].BankEle[j]-1].node[k+1];
										fprintf(Lakeboundarynode, "%lf\t%lf\n", DS->Node[StartNode - 1].x,DS->Node[StartNode - 1].y);
										LD->Lake[i].NumLakeNode = LD->Lake[i].NumLakeNode + 1;
									}										
								}
								if (k==2)
								{
									if (DS->Ele[LD->Lake[i].BankEle[j]-1].nabr[k-2]==0)
									{
										fprintf(Lakeboundarynode, "%lf\t%lf\n", DS->Node[DS->Ele[LD->Lake[i].BankEle[j]-1].node[k-1] - 1].x,DS->Node[DS->Ele[LD->Lake[i].BankEle[j]-1].node[k-1] - 1].y);
										LD->Lake[i].NumLakeNode = LD->Lake[i].NumLakeNode + 1;
										StartNode = DS->Ele[LD->Lake[i].BankEle[j]-1].node[k-2];
										fprintf(Lakeboundarynode, "%lf\t%lf\n", DS->Node[StartNode - 1].x,DS->Node[StartNode - 1].y);
										LD->Lake[i].NumLakeNode = LD->Lake[i].NumLakeNode + 1;
									}									
										
									if (DS->Ele[LD->Lake[i].BankEle[j]-1].nabr[k-1]==0)
									{
										fprintf(Lakeboundarynode, "%lf\t%lf\n", DS->Node[DS->Ele[LD->Lake[i].BankEle[j]-1].node[k-2] - 1].x,DS->Node[DS->Ele[LD->Lake[i].BankEle[j]-1].node[k-2] - 1].y);
										LD->Lake[i].NumLakeNode = LD->Lake[i].NumLakeNode + 1;
										StartNode = DS->Ele[LD->Lake[i].BankEle[j]-1].node[k-1];
										fprintf(Lakeboundarynode, "%lf\t%lf\n", DS->Node[StartNode - 1].x,DS->Node[StartNode - 1].y);
										LD->Lake[i].NumLakeNode = LD->Lake[i].NumLakeNode + 1;
									}										
								}
								LD->Lake[i].BankEle[j] = LD->Lake[i].BankEle[count];
								LD->Lake[i].BankEle[count] = FindIndex;
								count++;
								find = 1;
								printf("	--%d",StartNode);
								break;
							}						
						}
					}
				}
				if (find==1)
				{
					find =0;
					edge = 0;
					break;
				}
				edge = 0;
			}
		}
		
		count = 0;
        count1 = 0;

        
		for (j=0;j<LD->Lake[i].NumBankEle;j++)
		{
			for (k=0;k<DS->NumRiv;k++)
			{
                for (t=0;t<3;t++)
				{
					if (DS->Riv[k].FromNode == DS->Ele[LD->Lake[i].BankEle[j]-1].node[t])
					{
                        if (LD->Lake[i].StreamIndex[j]!=1)
						{
							count1 +=1;
							LD->Lake[i].StreamIndex[j] = 1;
						}						
					}
 
					if (DS->Riv[k].ToNode == DS->Ele[LD->Lake[i].BankEle[j]-1].node[t])
					{
                        if (LD->Lake[i].StreamIndex[j]!=1)
						{
							count1 +=1;
							LD->Lake[i].StreamIndex[j] = 1;
						}
					}
                }
			}
		}
		
		LD->Lake[i].NumStreamEle = count1;
		count1=0;

		LD->FluxSurf[i] = (double *)malloc((LD->Lake[i].NumBankEle)*sizeof(double));
		LD->FluxStream[i] = (double *)malloc(LD->Lake[i].NumBankEle*sizeof(double));
		LD->FluxSub[i] = (double *)malloc((LD->Lake[i].NumBankEle)*sizeof(double));
		LD->BankSurf[i] = (double *)malloc((LD->Lake[i].NumBankEle)*sizeof(double));
		LD->BankGW[i] = (double *)malloc((LD->Lake[i].NumBankEle)*sizeof(double));
	}
	
    for(i=0;i<LD->NumLake;i++)
	{
		fprintf(Lakeboundarynode, "%d\n", LD->Lake[i].NumLakeNode);
		
	}
	fflush(Lakeboundarynode);
	fclose(Lakeboundarynode);
	
    
    count2 = 0;
    for(i=0;i<LD->NumLake;i++)
    {
        for (q=0; q<LD->NumRiverInterest; q++)
        {
            for (j=0;j<LD->Lake[i].NumBankEle;j++)
            {
                for (t=0;t<3;t++) 
                {
                    if (DS->Riv[LD->RiverInterest[q]-1].FromNode == DS->Ele[LD->Lake[i].BankEle[j]-1].node[t])
                    {
                        if (count2==1||count2==2)
                        {
                            if (LD->Lake[i].BankEle[j] == LD->RiverInterestInformation[i][q][count2-1])
                            {
                                count2 = count2-1;
                            }
                        }
                        LD->RiverInterestInformation[i][q][count2]=LD->Lake[i].BankEle[j];
                        count2+=1;
                    }
                    
                    if (DS->Riv[LD->RiverInterest[q]-1].ToNode == DS->Ele[LD->Lake[i].BankEle[j]-1].node[t])
                    {
                        if (count2==1||count2==2)
                        {
                            if (LD->Lake[i].BankEle[j] == LD->RiverInterestInformation[i][q][count2-1])
                            {
                                count2 = count2-1;
                            }
                        }
                        LD->RiverInterestInformation[i][q][count2]=LD->Lake[i].BankEle[j];
                        count2+=1;
                    }
                }
            }
            count2=0;
        }
    }
	
    for (i = 0; i <2; i++) {
        
        LD->PrintVar[i] = (double *) calloc(LD->NumLake, sizeof(double));
    }
    
    for (i = 2; i <5; i++) {
        
        LD->PrintVar[i] = (double *) calloc(2*(LD->NumLake), sizeof(double));
    }
    
    for (i = 5; i <9; i++) {
        
        LD->PrintVar[i] = (double *) calloc(LD->NumLake, sizeof(double));
    }
    LD->PrintVar[9] = (double *) calloc(((LD->NumRiverInterest+4)*LD->NumLake), sizeof(double));
    
	
	for(i=0;i<LD->NumLake;i++)
	{
		LD->Lake[i].SurfIC = DS->Ele[LD->Lake[i].BankEle[0]-1].zmax;
		fprintf(LakeReference,"%d\t %lf\t\n",i+1, LD->Lake[i].SurfIC);
	}
	fflush(LakeReference);
	fclose(LakeReference);
	
	
	if(CS->init_type==0)
	{
		for(i=0;i<LD->NumLake;i++)
		{
			LD->LakeSurfDepth[i] = LD->Lake[i].SurfIC-LD->Lake[i].BedElev;
			LD->LakeGW[i] = LD->Lake[i].SurfIC-LD->Lake[i].BaseElev;
		}
	}
	else if(CS->init_type==2)
	{
		fn = (char *) malloc(1024 * sizeof(char));
		strcpy(fn, filename);
		init_file = fopen(strcat(fn, ".lakeinit"), "r");
		free(fn);
		if (init_file == NULL) 
		{
			printf("\n  Fatal Error: %s.init is in use or does not exist!\n", filename);
			exit(1);
		} 
		else 
		{
			for (i = 0; i < LD->NumLake; i++) 
			{
				fscanf(init_file, "%lf %lf", &LD->LakeSurfDepth[i], &LD->LakeGW[i]);
			}
		}
		fclose(init_file);	
	}
    
	for (i=0; i<LD->NumLake; i++)
    {
		DS->LakeInitialElev[i] = LD->Lake[i].SurfIC;
        DS->LakeOutletDecline[i] = LD->Lake[i].OutElevDecline;
        DS->BankElevAdjust[i] = LD->Lake[i].BankElevAdjust;
		DS->LakeBed[i] = LD->Lake[i].BedElev;
		DS->Orifice_height[i] = LD->Lake[i].Orifice_height;
		DS->Weir_length[i] = LD->Lake[i].Weir_length;
	}
    
    for (i=0; i<LD->NumLake; i++)
    {
        for (j=0; j<LD->NumRiverInterest; j++)
        {
            printf("\n There are %d Rivers with Users' interests", LD->NumRiverInterest);
            printf("\n The No  %d River segment with the elements %d  and %d  ", LD->RiverInterest[j],LD->RiverInterestInformation[i][j][0], LD->RiverInterestInformation[i][j][1]);
        }
        
        for (j=0;j<LD->Lake[i].NumBankEle;j++)
        {
            if (LD->Lake[i].StreamIndex[j]==1)
            {
                printf("\n The stream elements are %d", LD->Lake[i].BankEle[j]);
            }
        }
    }
}








