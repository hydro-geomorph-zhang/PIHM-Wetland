/*********************************************************************************
 * File        : read_alloc.c                                                    *
 * Function    : read parameter files for PIHM 2.0                     	         *
 * Version     : Nov, 2007 (2.0)                                                 *
 * Developer of PIHM2.0:	Mukesh Kumar (muk139@psu.edu)		         * 
 * Developer of PIHM1.0:	Yizhong Qu   (quyizhong@gmail.com)	         * 
 *-------------------------------------------------------------------------------*
 *                                                                               *
 *..............MODIFICATIONS/ADDITIONS in PIHM 2.0..............................*
 * a) Addition of three new input files: file.calib, file.lc and file.geol       *
 * b) Declaration and allocation  of new variables for new process, shape 	 *
 *    representations  and calibration (e.g. ET, Infiltration, Macropore, 	 *
 *    Stormflow, Element beneath river, river shapes, river bed property, 	 *
 *    thresholds for root zone, infiltration and macropore depths, land cover    * 
 *    attributes etc)                                                            *
 *--------------------------------------------------------------------------------*
 * For questions or comments, please contact                                      *
 *      --> Mukesh Kumar (muk139@psu.edu)                                         *
 *      --> Prof. Chris Duffy (cxd11@psu.edu)                                     *
 * This code is free for research purpose only.                                   *
 * Please provide relevant references if you use this code in your research work  *
 *--------------------------------------------------------------------------------*
 **********************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

//#include "sundialstypes.h"
#include "pihm.h"  
#include "lake.h"

void read_alloc(char *filename, Model_Data DS, Control_Data *CS)
{
    int i, j;
    int tempindex;

    int NumTout;
    char *fn[13];
    char tempchar[100];

    FILE *bid_file; /* Pointer to .bid file */
    FILE *mesh_file; /* Pointer to .mesh file */
    FILE *att_file; /* Pointer to .att file */
    FILE *forc_file; /* Pointer to .forc file*/
    FILE *ibc_file; /* Pointer to .ibc file*/
    FILE *soil_file; /* Pointer to .soil file */
    FILE *geol_file; /* Pointer to .geol file */
    FILE *lc_file; /* Pointer to .lc file */
    FILE *para_file; /* Pointer to .para file*/
    FILE *riv_file; /* Pointer to .riv file */
    FILE *global_calib; /* Pointer to .calib file */
    FILE *tide_file; /* Pointer to .calib file */
    FILE *centroid_elev; /* Pointer to .calib file */


    printf("\nStart reading in input files ... \n");

    /*========== open *.riv file ==========*/
    //  	printf("\n  1) reading %s.riv  ... ", filename);
    fn[0] = (char *) malloc((strlen(filename) + 1024) * sizeof (char));
    strcpy(fn[0], filename);
    riv_file = fopen(strcat(fn[0], ".riv"), "r");
    free(fn[0]);
    if (riv_file == NULL)
    {
        //   		printf("\n  Fatal Error: %s.riv is in use or does not exist!\n", filename);
        exit(1);
    }

    /* start reading riv_file */
    fscanf(riv_file, "%d", &(DS->NumRiv));

    DS->Riv = (river_segment *) malloc(DS->NumRiv * sizeof (river_segment));

    for (i = 0; i < DS->NumRiv; i++)
    {
        fscanf(riv_file, "%s", tempchar);
        fscanf(riv_file, "%d %d", &(DS->Riv[i].FromNode), &(DS->Riv[i].ToNode));
        fscanf(riv_file, "%d", &(DS->Riv[i].down));
        fscanf(riv_file, "%d %d", &(DS->Riv[i].LeftEle), &(DS->Riv[i].RightEle));
        fscanf(riv_file, "%s %s", tempchar, tempchar);
        fscanf(riv_file, "%s %s", tempchar, tempchar);
        fscanf(riv_file, "%s", tempchar);
    }
    fclose(riv_file);
    //  	printf("done.\n");



    /*========== open *.bid file ==========*/
    //  	printf("\n  2) reading %s.bid ... ", filename);
    fn[10] = (char *) malloc((strlen(filename) + 1024) * sizeof (char));
    strcpy(fn[10], filename);
    bid_file = fopen(strcat(fn[10], ".bid"), "r");
    free(fn[10]);
    if (bid_file == NULL)
    {
        printf("\n  Fatal Error: %s.mesh is in use or does not exist!\n", filename);
        exit(1);
    }

    /* start reading bid_file */
    fscanf(bid_file, "%d", &DS->NumOutlet);

    DS->Outlet_location = (int **) malloc(DS->NumOutlet * sizeof (int*));
    for (i = 0; i < DS->NumOutlet; i++)
    {
        DS->Outlet_location[i] = (int *) malloc(4 * sizeof (int));
    }

    for (i = 0; i < DS->NumOutlet; i++)
    {
        fscanf(bid_file, "%d", &DS->Outlet_location[i][0]);
        DS->Outlet_location[i][1] = -1;
        DS->Outlet_location[i][2] = -1;
        DS->Outlet_location[i][3] = -1;
    }
    //  	printf("done.\n");

    /* finish reading bid_files */
    fclose(bid_file);

    /*========== open *.mesh file ==========*/
    //  	printf("\n  2) reading %s.mesh ... ", filename);
    fn[1] = (char *) malloc((strlen(filename) + 1024) * sizeof (char));
    strcpy(fn[1], filename);
    mesh_file = fopen(strcat(fn[1], ".mesh"), "r");
    free(fn[1]);
    if (mesh_file == NULL)
    {
        printf("\n  Fatal Error: %s.mesh is in use or does not exist!\n", filename);
        exit(1);
    }

    /* start reading mesh_file */
    fscanf(mesh_file, "%d %d", &DS->NumEle, &DS->NumNode);

    DS->Ele = (element *) malloc((DS->NumEle) * sizeof (element));
    DS->Land_lake_index = (int *) malloc((DS->NumEle) * sizeof (int));
    DS->Node = (nodes *) malloc(DS->NumNode * sizeof (nodes));

    /* read in elements information */
    for (i = 0; i < DS->NumEle; i++)
    {
        fscanf(mesh_file, "%d", &(DS->Ele[i].index));
        fscanf(mesh_file, "%d %d %d", &(DS->Ele[i].node[0]), &(DS->Ele[i].node[1]), &(DS->Ele[i].node[2]));
        fscanf(mesh_file, "%d %d %d", &(DS->Ele[i].nabr[0]), &(DS->Ele[i].nabr[1]), &(DS->Ele[i].nabr[2]));
        DS->Land_lake_index[i] = 0;
    }

    /* read in nodes information */
    for (i = 0; i < DS->NumNode; i++)
    {
        fscanf(mesh_file, "%d", &(DS->Node[i].index));
        fscanf(mesh_file, "%lf %lf", &(DS->Node[i].x), &(DS->Node[i].y));
        fscanf(mesh_file, "%lf %lf", &(DS->Node[i].zmin), &(DS->Node[i].zmax));
    }

    //  	printf("done.\n");

    /* finish reading mesh_files */
    fclose(mesh_file);

    /*========== open *.att file ==========*/
    //  	printf("\n  3) reading %s.att  ... ", filename);
    fn[2] = (char *) malloc((strlen(filename) + 1024) * sizeof (char));
    strcpy(fn[2], filename);
    att_file = fopen(strcat(fn[2], ".att"), "r");
    free(fn[2]);
    if (att_file == NULL)
    {
        printf("\n  Fatal Error: %s.att is in use or does not exist!\n", filename);
        exit(1);
    }

    /* start reading att_file */
    DS->Ele_IC = (element_IC *) malloc(DS->NumEle * sizeof (element_IC));
    for (i = 0; i < DS->NumEle; i++)
    {
        fscanf(att_file, "%d", &(tempindex));
        fscanf(att_file, "%d %d %d", &(DS->Ele[i].soil), &(DS->Ele[i].geol), &(DS->Ele[i].LC));
        fscanf(att_file, "%lf %lf %lf %lf %lf", &(DS->Ele_IC[i].interception), &(DS->Ele_IC[i].snow), &(DS->Ele_IC[i].surf), &(DS->Ele_IC[i].unsat), &(DS->Ele_IC[i].sat));
        fscanf(att_file, "%d %d", &(DS->Ele[i].prep), &(DS->Ele[i].temp));
        fscanf(att_file, "%d %d", &(DS->Ele[i].humidity), &(DS->Ele[i].WindVel));
        fscanf(att_file, "%d %d", &(DS->Ele[i].Rn), &(DS->Ele[i].G));
        fscanf(att_file, "%d %d %d", &(DS->Ele[i].pressure), &(DS->Ele[i].source), &(DS->Ele[i].meltF));
        for (j = 0; j < 3; j++)
        {
            fscanf(att_file, "%d", &(DS->Ele[i].BC[j]));
        }
        fscanf(att_file, "%d", &(DS->Ele[i].Macropore));
    }

    //  	printf("done.\n");

    /* finish reading att_files */
    fclose(att_file);

    /*========== open *.soil file ==========*/
    //  	printf("\n  4) reading %s.soil ... ", filename);
    fn[3] = (char *) malloc((strlen(filename) + 1024) * sizeof (char));
    strcpy(fn[3], filename);
    soil_file = fopen(strcat(fn[3], ".soil"), "r");
    free(fn[3]);
    if (soil_file == NULL)
    {
        printf("\n  Fatal Error: %s.soil is in use or does not exist!\n", filename);
        exit(1);
    }

    /* start reading soil_file */
    fscanf(soil_file, "%d", &DS->NumSoil);
    DS->Soil = (soils *) malloc(DS->NumSoil * sizeof (soils));

    for (i = 0; i < DS->NumSoil; i++)
    {
        fscanf(soil_file, "%d", &(DS->Soil[i].index));
        /* Note: Soil KsatH and macKsatH is not used in model calculation anywhere */
        fscanf(soil_file, "%lf", &(DS->Soil[i].KsatV));
        fscanf(soil_file, "%lf %lf %lf", &(DS->Soil[i].ThetaS), &(DS->Soil[i].ThetaR), &(DS->Soil[i].infD));
        fscanf(soil_file, "%lf %lf", &(DS->Soil[i].Alpha), &(DS->Soil[i].Beta));
        fscanf(soil_file, "%lf %lf", &(DS->Soil[i].hAreaF), &(DS->Soil[i].macKsatV));
    }

    fclose(soil_file);
    //  	printf("done.\n");



    /*========== open *.geol file ==========*/
    //        printf("\n  5) reading %s.geol ... ", filename);
    fn[4] = (char *) malloc((strlen(filename) + 1024) * sizeof (char));
    strcpy(fn[4], filename);
    geol_file = fopen(strcat(fn[4], ".geol"), "r");
    free(fn[4]);
    if (geol_file == NULL)
    {
        printf("\n  Fatal Error: %s.geol is in use or does not exist!\n", filename);
        exit(1);
    }

    /* start reading*/
    fscanf(geol_file, "%d", &DS->NumGeol);
    DS->Geol = (geol *) malloc(DS->NumGeol * sizeof (geol));

    for (i = 0; i < DS->NumGeol; i++)
    {
        fscanf(geol_file, "%d", &(DS->Geol[i].index));
        /* Geol macKsatV is not used in model calculation anywhere */
        fscanf(geol_file, "%lf %lf", &(DS->Geol[i].KsatH), &(DS->Geol[i].KsatV));
        fscanf(geol_file, "%lf %lf", &(DS->Geol[i].ThetaS), &(DS->Geol[i].ThetaR));
        fscanf(geol_file, "%lf %lf", &(DS->Geol[i].Alpha), &(DS->Geol[i].Beta));
        fscanf(geol_file, "%lf %lf %lf", &(DS->Geol[i].vAreaF), &(DS->Geol[i].macKsatH), &(DS->Geol[i].macD));
    }
    fclose(geol_file);
    //        printf("done.\n");

    /*========== open *.lc file ==========*/
    //  	printf("\n  6) reading %s.lc ... ", filename);
    fn[5] = (char *) malloc((strlen(filename) + 1024) * sizeof (char));
    strcpy(fn[5], filename);
    lc_file = fopen(strcat(fn[5], ".lc"), "r");
    free(fn[5]);
    if (lc_file == NULL)
    {
        printf("\n  Fatal Error: %s.land cover is in use or does not exist!\n", filename);
        exit(1);
    }

    /* start reading land cover file */
    fscanf(lc_file, "%d", &DS->NumLC);

    DS->LandC = (LC *) malloc(DS->NumLC * sizeof (LC));

    for (i = 0; i < DS->NumLC; i++)
    {
        fscanf(lc_file, "%d", &(DS->LandC[i].index));
        fscanf(lc_file, "%lf", &(DS->LandC[i].LAImax));
        fscanf(lc_file, "%lf %lf", &(DS->LandC[i].Rmin), &(DS->LandC[i].Rs_ref));
        fscanf(lc_file, "%lf %lf", &(DS->LandC[i].Albedo), &(DS->LandC[i].VegFrac));
        fscanf(lc_file, "%lf %lf", &(DS->LandC[i].Rough), &(DS->LandC[i].RzD));
    }

    fclose(lc_file);
    // 	printf("done.\n");


    /*========== open *.forc file ==========*/
    //	printf("\n  7) reading %s.forc ... ", filename);
    fn[6] = (char *) malloc((strlen(filename) + 1024) * sizeof (char));
    strcpy(fn[6], filename);
    forc_file = fopen(strcat(fn[6], ".forc"), "r");
    free(fn[6]);
    if (forc_file == NULL)
    {
        printf("\n  Fatal Error: %s.forc is in use or does not exist!\n", filename);
        exit(1);
    }

    /* start reading forc_file */
    fscanf(forc_file, "%d %d", &DS->NumPrep, &DS->NumTemp);
    fscanf(forc_file, "%d %d", &DS->NumHumidity, &DS->NumWindVel);
    fscanf(forc_file, "%d %d", &DS->NumRn, &DS->NumG);
    fscanf(forc_file, "%d %d", &DS->NumP, &DS->NumLC);
    fscanf(forc_file, "%d", &DS->NumMeltF);
    fscanf(forc_file, "%d", &DS->NumSource);

    DS->TSD_Prep = (TSD *) malloc(DS->NumPrep * sizeof (TSD));
    DS->TSD_Temp = (TSD *) malloc(DS->NumTemp * sizeof (TSD));
    DS->TSD_Humidity = (TSD *) malloc(DS->NumHumidity * sizeof (TSD));
    DS->TSD_WindVel = (TSD *) malloc(DS->NumWindVel * sizeof (TSD));
    DS->TSD_Rn = (TSD *) malloc(DS->NumRn * sizeof (TSD));
    DS->TSD_G = (TSD *) malloc(DS->NumG * sizeof (TSD));
    DS->TSD_Pressure = (TSD *) malloc(DS->NumP * sizeof (TSD));
    DS->TSD_LAI = (TSD *) malloc(DS->NumLC * sizeof (TSD));
    DS->TSD_RL = (TSD *) malloc(DS->NumLC * sizeof (TSD));
    DS->TSD_MeltF = (TSD *) malloc(DS->NumMeltF * sizeof (TSD));
    DS->TSD_Source = (TSD *) malloc(DS->NumSource * sizeof (TSD));

    DS->ISFactor = (double *) malloc(DS->NumLC * sizeof (double));
    DS->windH = (double *) malloc(DS->NumWindVel * sizeof (double));

    for (i = 0; i < DS->NumPrep; i++)
    {
        fscanf(forc_file, "%s %d %d", DS->TSD_Prep[i].name, &DS->TSD_Prep[i].index, &DS->TSD_Prep[i].length);

        DS->TSD_Prep[i].TS = (double **) malloc((DS->TSD_Prep[i].length) * sizeof (double*));

        for (j = 0; j < DS->TSD_Prep[i].length; j++)
        {
            DS->TSD_Prep[i].TS[j] = (double *) malloc(2 * sizeof (double));
        }
        DS->TSD_Prep[i].iCounter = 0;
        for (j = 0; j < DS->TSD_Prep[i].length; j++)
        {
            fscanf(forc_file, "%lf %lf", &DS->TSD_Prep[i].TS[j][0], &DS->TSD_Prep[i].TS[j][1]);
        }
    }

    for (i = 0; i < DS->NumTemp; i++)
    {
        fscanf(forc_file, "%s %d %d", DS->TSD_Temp[i].name, &DS->TSD_Temp[i].index, &DS->TSD_Temp[i].length);

        DS->TSD_Temp[i].TS = (double **) malloc((DS->TSD_Temp[i].length) * sizeof (double*));

        for (j = 0; j < DS->TSD_Temp[i].length; j++)
        {
            DS->TSD_Temp[i].TS[j] = (double *) malloc(2 * sizeof (double));
        }
        DS->TSD_Temp[i].iCounter = 0;
        for (j = 0; j < DS->TSD_Temp[i].length; j++)
        {
            fscanf(forc_file, "%lf %lf", &DS->TSD_Temp[i].TS[j][0], &DS->TSD_Temp[i].TS[j][1]);
        }
    }

    for (i = 0; i < DS->NumHumidity; i++)
    {
        fscanf(forc_file, "%s %d %d", DS->TSD_Humidity[i].name, &DS->TSD_Humidity[i].index, &DS->TSD_Humidity[i].length);

        DS->TSD_Humidity[i].TS = (double **) malloc((DS->TSD_Humidity[i].length) * sizeof (double*));

        for (j = 0; j < DS->TSD_Humidity[i].length; j++)
        {
            DS->TSD_Humidity[i].TS[j] = (double *) malloc(2 * sizeof (double));
        }
        DS->TSD_Humidity[i].iCounter = 0;
        for (j = 0; j < DS->TSD_Humidity[i].length; j++)
        {
            fscanf(forc_file, "%lf %lf", &DS->TSD_Humidity[i].TS[j][0], &DS->TSD_Humidity[i].TS[j][1]);
        }
    }

    for (i = 0; i < DS->NumWindVel; i++)
    {
        fscanf(forc_file, "%s %d %d %lf", DS->TSD_WindVel[i].name, &DS->TSD_WindVel[i].index, &DS->TSD_WindVel[i].length, &DS->windH[i]);
        DS->TSD_WindVel[i].TS = (double **) malloc((DS->TSD_WindVel[i].length) * sizeof (double*));
        DS->TSD_WindVel[i].iCounter = 0;
        for (j = 0; j < DS->TSD_WindVel[i].length; j++)
        {
            DS->TSD_WindVel[i].TS[j] = (double *) malloc(2 * sizeof (double));
        }

        for (j = 0; j < DS->TSD_WindVel[i].length; j++)
        {
            fscanf(forc_file, "%lf %lf", &DS->TSD_WindVel[i].TS[j][0], &DS->TSD_WindVel[i].TS[j][1]);
        }
    }

    for (i = 0; i < DS->NumRn; i++)
    {
        fscanf(forc_file, "%s %d %d", DS->TSD_Rn[i].name, &DS->TSD_Rn[i].index, &DS->TSD_Rn[i].length);
        DS->TSD_Rn[i].iCounter = 0;
        DS->TSD_Rn[i].TS = (double **) malloc((DS->TSD_Rn[i].length) * sizeof (double*));

        for (j = 0; j < DS->TSD_Rn[i].length; j++)
        {
            DS->TSD_Rn[i].TS[j] = (double *) malloc(2 * sizeof (double));
        }

        for (j = 0; j < DS->TSD_Rn[i].length; j++)
        {
            fscanf(forc_file, "%lf %lf", &DS->TSD_Rn[i].TS[j][0], &DS->TSD_Rn[i].TS[j][1]);
        }
    }

    for (i = 0; i < DS->NumG; i++)
    {
        fscanf(forc_file, "%s %d %d", DS->TSD_G[i].name, &DS->TSD_G[i].index, &DS->TSD_G[i].length);
        DS->TSD_G[i].iCounter = 0;
        DS->TSD_G[i].TS = (double **) malloc((DS->TSD_G[i].length) * sizeof (double*));

        for (j = 0; j < DS->TSD_G[i].length; j++)
        {
            DS->TSD_G[i].TS[j] = (double *) malloc(2 * sizeof (double));
        }

        for (j = 0; j < DS->TSD_G[i].length; j++)
        {
            fscanf(forc_file, "%lf %lf", &DS->TSD_G[i].TS[j][0], &DS->TSD_G[i].TS[j][1]);
        }
    }

    for (i = 0; i < DS->NumP; i++)
    {
        fscanf(forc_file, "%s %d %d", DS->TSD_Pressure[i].name, &DS->TSD_Pressure[i].index, &DS->TSD_Pressure[i].length);
        DS->TSD_Pressure[i].iCounter = 0;
        DS->TSD_Pressure[i].TS = (double **) malloc((DS->TSD_Pressure[i].length) * sizeof (double*));

        for (j = 0; j < DS->TSD_Pressure[i].length; j++)
        {
            DS->TSD_Pressure[i].TS[j] = (double *) malloc(2 * sizeof (double));
        }

        for (j = 0; j < DS->TSD_Pressure[i].length; j++)
        {
            fscanf(forc_file, "%lf %lf", &DS->TSD_Pressure[i].TS[j][0], &DS->TSD_Pressure[i].TS[j][1]);
        }
    }

    for (i = 0; i < DS->NumLC; i++)
    {
        fscanf(forc_file, "%s %d %d %lf", DS->TSD_LAI[i].name, &DS->TSD_LAI[i].index, &DS->TSD_LAI[i].length, &DS->ISFactor[i]);
        DS->TSD_LAI[i].iCounter = 0;
        DS->TSD_LAI[i].TS = (double **) malloc((DS->TSD_LAI[i].length) * sizeof (double*));

        for (j = 0; j < DS->TSD_LAI[i].length; j++)
        {
            DS->TSD_LAI[i].TS[j] = (double *) malloc(2 * sizeof (double));
        }

        for (j = 0; j < DS->TSD_LAI[i].length; j++)
        {
            fscanf(forc_file, "%lf %lf", &DS->TSD_LAI[i].TS[j][0], &DS->TSD_LAI[i].TS[j][1]);
        }
    }

    for (i = 0; i < DS->NumLC; i++)
    {
        fscanf(forc_file, "%s %d %d", DS->TSD_RL[i].name, &DS->TSD_RL[i].index, &DS->TSD_RL[i].length);

        DS->TSD_RL[i].TS = (double **) malloc((DS->TSD_RL[i].length) * sizeof (double*));
        DS->TSD_RL[i].iCounter = 0;
        for (j = 0; j < DS->TSD_RL[i].length; j++)
        {
            DS->TSD_RL[i].TS[j] = (double *) malloc(2 * sizeof (double));
        }

        for (j = 0; j < DS->TSD_RL[i].length; j++)
        {
            fscanf(forc_file, "%lf %lf", &DS->TSD_RL[i].TS[j][0], &DS->TSD_RL[i].TS[j][1]);
        }
    }

    for (i = 0; i < DS->NumMeltF; i++)
    {
        fscanf(forc_file, "%s %d %d", DS->TSD_MeltF[i].name, &DS->TSD_MeltF[i].index, &DS->TSD_MeltF[i].length);
        DS->TSD_MeltF[i].iCounter = 0;
        DS->TSD_MeltF[i].TS = (double **) malloc((DS->TSD_MeltF[i].length) * sizeof (double*));

        for (j = 0; j < DS->TSD_MeltF[i].length; j++)
        {
            DS->TSD_MeltF[i].TS[j] = (double *) malloc(2 * sizeof (double));
        }

        for (j = 0; j < DS->TSD_MeltF[i].length; j++)
        {
            fscanf(forc_file, "%lf %lf", &DS->TSD_MeltF[i].TS[j][0], &DS->TSD_MeltF[i].TS[j][1]);
        }
    }

    for (i = 0; i < DS->NumSource; i++)
    {
        fscanf(forc_file, "%s %d %d", DS->TSD_Source[i].name, &DS->TSD_Source[i].index, &DS->TSD_Source[i].length);

        DS->TSD_Source[i].TS = (double **) malloc((DS->TSD_Source[i].length) * sizeof (double*));

        for (j = 0; j < DS->TSD_Source[i].length; j++)
        {
            DS->TSD_Source[i].TS[j] = (double *) malloc(2 * sizeof (double));
        }

        for (j = 0; j < DS->TSD_Source[i].length; j++)
        {
            fscanf(forc_file, "%lf %lf", &DS->TSD_Source[i].TS[j][0], &DS->TSD_Source[i].TS[j][1]);
        }
    }

    fclose(forc_file);
    //  	printf("done.\n");

    /*========== open *.ibc file ==========*/
    //  	printf("\n  8) reading %s.ibc  ... ", filename);  
    fn[7] = (char *) malloc((strlen(filename) + 1024) * sizeof (char));
    strcpy(fn[7], filename);
    ibc_file = fopen(strcat(fn[7], ".ibc"), "r");
    free(fn[7]);
    if (ibc_file == NULL)
    {
        printf("\n  Fatal Error: %s.ibc is in use or does not exist!\n", filename);
        exit(1);
    }

    /* start reading ibc_file */
    fscanf(ibc_file, "%d %d", &DS->Num1BC, &DS->Num2BC);

    if (DS->Num1BC + DS->Num2BC > 0)
    {
        DS->TSD_EleBC = (TSD *) malloc((DS->Num1BC + DS->Num2BC) * sizeof (TSD));
    }

    if (DS->Num1BC > 0)
    {
        /* For elements with Dirichilet Boundary Conditions */
        for (i = 0; i < DS->Num1BC; i++)
        {
            fscanf(ibc_file, "%s %d %d", DS->TSD_EleBC[i].name, &DS->TSD_EleBC[i].index, &DS->TSD_EleBC[i].length);

            DS->TSD_EleBC[i].TS = (double **) malloc((DS->TSD_EleBC[i].length) * sizeof (double*));

            for (j = 0; j < DS->TSD_EleBC[i].length; j++)
            {
                DS->TSD_EleBC[i].TS[j] = (double *) malloc(2 * sizeof (double));
            }

            for (j = 0; j < DS->TSD_EleBC[i].length; j++)
            {
                fscanf(forc_file, "%lf %lf", &DS->TSD_EleBC[i].TS[j][0], &DS->TSD_EleBC[i].TS[j][1]);
            }
        }
    }

    if (DS->Num2BC > 0)
    {
        /* For elements with Neumann (non-natural) Boundary Conditions */
        for (i = DS->Num1BC; i < DS->Num1BC + DS->Num2BC; i++)
        {
            fscanf(ibc_file, "%s %d %d", DS->TSD_EleBC[i].name, &DS->TSD_EleBC[i].index, &DS->TSD_EleBC[i].length);

            DS->TSD_EleBC[i].TS = (double **) malloc((DS->TSD_EleBC[i].length) * sizeof (double*));

            for (j = 0; j < DS->TSD_EleBC[i].length; j++)
            {
                DS->TSD_EleBC[i].TS[j] = (double *) malloc(2 * sizeof (double));
            }
            for (j = 0; j < DS->TSD_EleBC[i].length; j++)
            {
                fscanf(forc_file, "%lf %lf", &DS->TSD_EleBC[i].TS[j][0], &DS->TSD_EleBC[i].TS[j][1]);
            }
        }
    }
    fclose(ibc_file);
    // 	printf("done.\n");

    /*========== open *.para file ==========*/
    //	printf("\n  9) reading %s.para ... ", filename); 
    fn[8] = (char *) malloc((strlen(filename) + 1024) * sizeof (char));
    strcpy(fn[8], filename);
    para_file = fopen(strcat(fn[8], ".para"), "r");
    free(fn[8]);
    if (para_file == NULL)
    {
        printf("\n  Fatal Error: %s.para is in use or does not exist!\n", filename);
        exit(1);
    }

    /* start reading para_file */
    fscanf(para_file, "%s %d %s %d", tempchar, &(CS->Verbose), tempchar, &(CS->Debug));
    fscanf(para_file, "%s %d", tempchar, &(CS->init_type));
    fscanf(para_file, "%s %d", tempchar, &(CS->Lake_module));
    fscanf(para_file, "%s %d", tempchar, &(DS->Tide_mode));
    fscanf(para_file, "%s %d", tempchar, &(DS->Sea_level_rise_mode));
    fscanf(para_file, "%s %d %s %d %s %d %s %d", tempchar, &(CS->gwDInt), tempchar, &(CS->surfDInt), tempchar, &(CS->snowDInt), tempchar, &(CS->prepDInt));
    fscanf(para_file, "%s %d %s %d %s %d %s %d", tempchar, &(CS->RechInt), tempchar, &(CS->IsDInt), tempchar, &(CS->usDInt), tempchar, &(CS->etInt));

    fscanf(para_file, "%s %d %s %d", tempchar, &DS->UnsatMode, tempchar, &DS->SurfMode);
    fscanf(para_file, "%s %d", tempchar, &(CS->Solver));
    if (CS->Solver == 2)
    {
        fscanf(para_file, "%s %d %s %d %s %lf", tempchar, &CS->GSType, tempchar, &CS->MaxK, tempchar, &CS->delt);
    }
    fscanf(para_file, "%s %lf %s %lf", tempchar, &(CS->abstol), tempchar, &(CS->reltol));
    fscanf(para_file, "%s %lf %s %lf %s %lf", tempchar, &(CS->InitStep), tempchar, &(CS->MaxStep), tempchar, &(CS->ETStep));
    fscanf(para_file, "%s %lf %s %lf %s %lf %s %d", tempchar, &(CS->StartTime), tempchar, &(CS->EndTime), tempchar, &(CS->FinishTime), tempchar, &(CS->outtype));
    if (CS->outtype == 0)
    {
        fscanf(para_file, "%s %lf %s %lf", tempchar, &CS->a, tempchar, &CS->b);
    }

    if (CS->a != 1.0)
    {
        NumTout = (int) (log(1 - (CS->EndTime - CS->StartTime)*(1 - CS->a) / CS->b) / log(CS->a));
    }
    else
    {
        if ((CS->EndTime - CS->StartTime) / CS->b - ((int) (CS->EndTime - CS->StartTime) / CS->b) > 0)
        {
            NumTout = (int) ((CS->EndTime - CS->StartTime) / CS->b);
        }
        else
        {
            NumTout = (int) ((CS->EndTime - CS->StartTime) / CS->b - 1);
        }
    }

    CS->NumSteps = NumTout + 1;

    CS->Tout = (double *) malloc((CS->NumSteps + 1) * sizeof (double));

    for (i = 0; i < CS->NumSteps + 1; i++)
    {
        if (i == 0)
        {
            CS->Tout[i] = CS->StartTime;
        }
        else
        {
            CS->Tout[i] = CS->Tout[i - 1] + pow(CS->a, i) * CS->b;
        }
    }

    if (CS->Tout[CS->NumSteps] < CS->EndTime)
    {
        CS->Tout[CS->NumSteps] = CS->EndTime;
    }

    fclose(para_file);
    // 	printf("done.\n"); 

    //	printf("\nStart reading in calibration file...\n");

    /*========= open *.calib file ==========*/
    //	printf("\n  10) reading %s.calib ... ", filename);
    fn[9] = (char *) malloc((strlen(filename) + 1024) * sizeof (char));
    strcpy(fn[9], filename);
    global_calib = fopen(strcat(fn[9], ".calib"), "r");
    free(fn[9]);
    if (global_calib == NULL)
    {
        printf("\n  Fatal Error: %s.calib is in use or does not exist!\n", filename);
        exit(1);
    }

    /* start reading calib_file */
    fscanf(global_calib, "%s %lf %s %lf %s %lf %s %lf %s %lf", tempchar, &CS->Cal.KsatH, tempchar, &CS->Cal.KsatV, tempchar, &CS->Cal.infKsatV, tempchar, &CS->Cal.macKsatH, tempchar, &CS->Cal.macKsatV);
    fscanf(global_calib, "%s %lf %s %lf %s %lf", tempchar, &CS->Cal.infD, tempchar, &CS->Cal.RzD, tempchar, &CS->Cal.macD);
    fscanf(global_calib, "%s %lf %s %lf %s %lf", tempchar, &CS->Cal.Porosity, tempchar, &CS->Cal.Alpha, tempchar, &CS->Cal.Beta);
    fscanf(global_calib, "%s %lf %s %lf", tempchar, &CS->Cal.vAreaF, tempchar, &CS->Cal.hAreaF);
    fscanf(global_calib, "%s %lf %s %lf %s %lf", tempchar, &CS->Cal.VegFrac, tempchar, &CS->Cal.Albedo, tempchar, &CS->Cal.Rough);
    fscanf(global_calib, "%s %lf %s %lf", tempchar, &CS->Cal.Prep, tempchar, &CS->Cal.Temp);
    fscanf(global_calib, "%s %lf %s %lf %s %lf", tempchar, &DS->pcCal.Et0, tempchar, &DS->pcCal.Et1, tempchar, &DS->pcCal.Et2);
    fscanf(global_calib, "%s %lf %s %lf %s %lf %s %lf %s %lf", tempchar, &CS->Cal.LakeKh, tempchar, &CS->Cal.LakeKv, tempchar, &CS->Cal.LakeThetaS, tempchar, &CS->Cal.LakeThetaR, tempchar, &CS->Cal.LakemacKh);

    // 	printf("done.\n");

    /* finish reading calib file */
    fclose(global_calib);
    
    
     /*========= open *.centroid file ==========*/
    //	//printf("\n  10) reading %s.centroid ... ", filename);
    fn = (char *) malloc((strlen(filename) + 1024) * sizeof (char));
    strcpy(fn, filename);
    centroid_elev = fopen(strcat(fn, ".centroid"), "r");
    free(fn);
    if (centroid_elev == NULL)
    {
        CS->read_centroid = 0;
    }
    else
    {
        CS->read_centroid = 1;
        /* start reading centroid elevation file */
        for (i = 0; i < DS->NumEle; i++)
        {
            fscanf(centroid_elev, "%lf %lf", &DS->Ele[i].zmax, &DS->Ele[i].zmin);
        }
        fclose(centroid_elev);
    }
    
    

    /*========= open *.tide file ==========*/
    //	printf("\n  10) reading %s.tide ... ", filename);
    if (DS->Tide_mode == 1)
    {
        fn[10] = (char *) malloc((strlen(filename) + 1024) * sizeof (char));
        strcpy(fn[10], filename);
        tide_file = fopen(strcat(fn[10], ".tide"), "r");
        free(fn[10]);
        if (tide_file == NULL)
        {
            printf("\n  Fatal Error: %s.tide is in use or does not exist!\n", filename);
            exit(1);
        }

        fscanf(tide_file, "%d", &DS->NumTide);

        DS->TSD_Tide = (TSD *) malloc(DS->NumTide * sizeof (TSD));

        for (i = 0; i < DS->NumTide; i++)
        {
            fscanf(tide_file, "%s %d %d", DS->TSD_Tide[i].name, &DS->TSD_Tide[i].index, &DS->TSD_Tide[i].length);

            DS->TSD_Tide[i].TS = (double **) malloc((DS->TSD_Prep[i].length) * sizeof (double*));

            for (j = 0; j < DS->TSD_Tide[i].length; j++)
            {
                DS->TSD_Tide[i].TS[j] = (double *) malloc(2 * sizeof (double));
            }
            DS->TSD_Tide[i].iCounter = 0;
            for (j = 0; j < DS->TSD_Tide[i].length; j++)
            {
                fscanf(tide_file, "%lf %lf", &DS->TSD_Tide[i].TS[j][0], &DS->TSD_Tide[i].TS[j][1]);
            }
        }
        /* finish reading calib file */
        fclose(tide_file);
    }

}

void FreeData(Model_Data DS, Control_Data * CS, LakeData LD)
{

    /*free river*/
    int i, j;

    /*free mesh*/

    free(DS->Ele);
    free(DS->Node);
    /*free att*/
    free(DS->Ele_IC);
    /*free soil*/
    free(DS->Soil);
    /*free geol*/
    free(DS->Geol);
    /*free lc*/
    free(DS->LandC);
    /*free forc*/
    for (i = 0; i < DS->NumPrep; i++)
    {
        for (j = 0; j < DS->TSD_Prep[i].length; j++)free(DS->TSD_Prep[i].TS[j]);
        free(DS->TSD_Prep[i].TS);
    }
    free(DS->TSD_Prep);
    for (i = 0; i < DS->NumTemp; i++)
    {
        for (j = 0; j < DS->TSD_Temp[i].length; j++)free(DS->TSD_Temp[i].TS[j]);
        free(DS->TSD_Temp[i].TS);
    }
    free(DS->TSD_Temp);
    for (i = 0; i < DS->NumHumidity; i++)
    {
        for (j = 0; j < DS->TSD_Humidity[i].length; j++)free(DS->TSD_Humidity[i].TS[j]);
        free(DS->TSD_Humidity[i].TS);
    }
    free(DS->TSD_Humidity);
    for (i = 0; i < DS->NumWindVel; i++)
    {
        for (j = 0; j < DS->TSD_WindVel[i].length; j++)free(DS->TSD_WindVel[i].TS[j]);
        free(DS->TSD_WindVel[i].TS);
    }
    free(DS->TSD_WindVel);
    free(DS->windH);
    for (i = 0; i < DS->NumRn; i++)
    {
        for (j = 0; j < DS->TSD_Rn[i].length; j++)free(DS->TSD_Rn[i].TS[j]);
        free(DS->TSD_Rn[i].TS);
    }
    free(DS->TSD_Rn);
    for (i = 0; i < DS->NumG; i++)
    {
        for (j = 0; j < DS->TSD_G[i].length; j++)free(DS->TSD_G[i].TS[j]);
        free(DS->TSD_G[i].TS);
    }
    free(DS->TSD_G);
    for (i = 0; i < DS->NumP; i++)
    {
        for (j = 0; j < DS->TSD_Pressure[i].length; j++)free(DS->TSD_Pressure[i].TS[j]);
        free(DS->TSD_Pressure[i].TS);
    }
    free(DS->TSD_Pressure);
    for (i = 0; i < DS->NumLC; i++)
    {
        for (j = 0; j < DS->TSD_LAI[i].length; j++)free(DS->TSD_LAI[i].TS[j]);
        for (j = 0; j < DS->TSD_RL[i].length; j++) free(DS->TSD_RL[i].TS[j]);

        free(DS->TSD_LAI[i].TS);
        free(DS->TSD_RL[i].TS);
    }
    free(DS->TSD_LAI);
    free(DS->TSD_RL);
    for (i = 0; i < DS->NumMeltF; i++)
    {
        for (j = 0; j < DS->TSD_MeltF[i].length; j++)free(DS->TSD_MeltF[i].TS[j]);
        free(DS->TSD_MeltF[i].TS);
    }
    free(DS->TSD_MeltF);
    for (i = 0; i < DS->NumSource; i++)
    {
        for (j = 0; j < DS->TSD_Source[i].length; j++)free(DS->TSD_Source[i].TS[j]);
        free(DS->TSD_Source[i].TS);
    }
    free(DS->TSD_Source);
    free(DS->ISFactor);
    /*free ibc*/
    if (DS->Num1BC > 0)
        for (i = 0; i < DS->Num1BC; i++)
        {
            for (j = 0; j < DS->TSD_EleBC[i].length; j++)free(DS->TSD_EleBC[i].TS[j]);
            free(DS->TSD_EleBC[i].TS);
        }
    if (DS->Num2BC > 0)
        for (i = DS->Num1BC; i < DS->Num1BC + DS->Num2BC; i++)
        {
            for (j = 0; j < DS->TSD_EleBC[i].length; j++)free(DS->TSD_EleBC[i].TS[j]);
            free(DS->TSD_EleBC[i].TS);
        }


    if (DS->Num1BC + DS->Num2BC > 0) free(DS->TSD_EleBC);
    /*free para*/
    free(CS->Tout);
    /*free initialize.c*/
    for (i = 0; i < DS->NumEle; i++)free(DS->FluxSurf[i]);
    free(DS->FluxSurf);
    for (i = 0; i < DS->NumEle; i++)free(DS->FluxSub[i]);
    free(DS->FluxSub);
    for (i = 0; i < DS->NumEle; i++)free(DS->EleET[i]);
    free(DS->EleET);
    for (i = 0; i < DS->NumEle; i++)free(DS->Seepage[i]);
    free(DS->Seepage);
    for (i = 0; i < DS->NumEle; i++)free(DS->Seepage[i]);
    free(DS->Seepage);
    for (i = 0; i < DS->NumEle; i++)free(DS->Courant_surf[i]);
    free(DS->Courant_surf);
    for (i = 0; i < DS->NumEle; i++)free(DS->Courant_sub[i]);
    free(DS->Courant_sub);
    for (i = 0; i < DS->NumEle; i++)free(DS->CrossA_surf[i]);
    free(DS->CrossA_surf);
    for (i = 0; i < DS->NumEle; i++)free(DS->CrossA_sub[i]);
    free(DS->CrossA_sub);
    for (i = 0; i < DS->NumEle; i++)free(DS->Distance[i]);
    free(DS->Distance);

    free(DS->ElePrep);
    free(DS->EleViR);
    free(DS->Recharge);
    free(DS->EleIS);
    free(DS->EleISmax);
    free(DS->EleISsnowmax);
    free(DS->EleSnow);
    free(DS->EleSnowGrnd);
    free(DS->EleSnowCanopy);
    free(DS->EleTF);
    free(DS->EleETloss);
    free(DS->EleNetPrep);
    free(DS->infil_mode);

    /*free Print*/
    for (i = 0; i < 24; i++)free(DS->PrintVar[i]);
    /*free DummyY*/
    free(DS->DummyY);
    if (CS->Lake_module != 0)
    {
        free(DS->LakeSurfElev);
        free(DS->LakeGWElev);
        free(DS->LakeInitialElev);
        free(DS->LakeOutletDecline);
        free(DS->LakeBed);
        free(DS->Weir_length);
        free(DS->Orifice_height);
        free(DS->BankElevAdjust);
        free(LD->LakeSurfDepth);
        free(LD->LakeGW);
        free(LD->LakeDummyY);
        free(LD->ET);
        free(LD->Infil);
        free(LD->RiverInterestInformation);
        free(LD->LakeSoil);
        for (i = 0; i < 10; i++)free(LD->PrintVar[i]);
        for (i = 0; i < LD->NumLake; i++) free(LD->BankSurf[i]);
        free(LD->BankSurf);
        for (i = 0; i < LD->NumLake; i++) free(LD->BankGW[i]);
        free(LD->BankGW);
        for (i = 0; i < LD->NumLake; i++) free(LD->FluxStream[i]);
        free(LD->FluxStream);
        for (i = 0; i < LD->NumLake; i++) free(LD->FluxSurf[i]);
        free(LD->FluxSurf);
        for (i = 0; i < LD->NumLake; i++) free(LD->FluxSub[i]);
        free(LD->FluxSub);
        for (i = 0; i < LD->NumLake; i++) free(LD->FluxSub[i]);
        free(LD->FluxSub);
        for (i = 0; i < LD->NumLake; i++) free(LD->LakeBoundary[i]);
        free(LD->LakeBoundary);
        for (i = 0; i < LD->NumLake; i++) free(LD->LakeBoundaryNodes[i]);
        free(LD->LakeBoundaryNodes);
        for (i = 0; i < LD->NumLake; i++)
        {
            for (j = 0; j < LD->NumRiverInterest; j++)
            {
                free(LD->RiverInterestInformation[i][j]);
            }
            free(LD->RiverInterestInformation[i]);
        }
        free(LD->RiverInterestInformation);

        for (i = 0; i < LD->NumLake; i++) free(LD->Lake[i].StreamEle);
        for (i = 0; i < LD->NumLake; i++) free(LD->Lake[i].BankEle);
        for (i = 0; i < LD->NumLake; i++) free(LD->Lake[i].StreamIndex);
        free(LD->Lake);
    }


}
