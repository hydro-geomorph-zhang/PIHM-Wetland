/*********************************************************************************
 * File        : initialize.c                                                    *
 * Function    : initialization of elemental attributes using relational database*
 * Version     : Nov, 2007 (2.0)                                                 *
 * Developer of PIHM2.0:        Mukesh Kumar (muk139@psu.edu)                    *
 * Developer of PIHM1.0:        Yizhong Qu                                       *
 *-------------------------------------------------------------------------------*
 *                                                                               *
 *..............MODIFICATIONS/ADDITIONS in PIHM 2.0..............................*
 * a) Correction of edge calculation term                                        *
 * b) Initialization of new variables for new process, shape representations  and*
 *    calibration (e.g. ET, Infiltration, Macropore, Stormflow, Element beneath  *
 *    river, river shapes, river bed property, thresholds for root zone,         *
 *    infiltration and macropore depths, land cover attributes etc)              *
 *--------------------------------------------------------------------------------*
 * For questions or comments, please contact                                      *
 *      --> Mukesh Kumar (muk139@psu.edu)                                         *
 *      --> Prof. Chris Duffy (cxd11@psu.edu)                                     *
 * This code is free for research purpose only.                                   *
 * Please provide relevant references if you use this code in your research work  *
 *--------------------------------------------------------------------------------*
 *********************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "sundials_types.h"
#include "nvector_serial.h"

#include "pihm.h"
#include "lake.h"

void initialize(char *filename, Model_Data DS, Control_Data * CS, FILE *elev_output)
{
    int i, j, k, t, tmpBool, BoolBR, BoolR = 0, inabr;
    double a_x, a_y, b_x, b_y, c_x, c_y, distX, distY;
    double a_zmin, a_zmax, b_zmin, b_zmax, c_zmin, c_zmax;
    double tempvalue1, tempvalue2, tempvalue3, tempvalue4, tempvalue5;
    FILE *init_file, *newelev_file;
    char *fn;
    double *zmin_cor;

    //	zmin_cor = (double *) malloc(DS->NumEle * sizeof(double));

    printf("\nInitializing data structure ... ");

    /* allocate memory storage to flux terms */
    DS->FluxSurf = (double **) malloc(DS->NumEle * sizeof (double*));
    DS->FluxSub = (double **) malloc(DS->NumEle * sizeof (double*));
    DS->Seepage = (double **) malloc(DS->NumEle * sizeof (double*));
    DS->EleET = (double **) malloc(DS->NumEle * sizeof (double*));
    DS->ElePrep = (double *) malloc(DS->NumEle * sizeof (double));
    DS->EleViR = (double *) malloc(DS->NumEle * sizeof (double));
    DS->Recharge = (double *) malloc(DS->NumEle * sizeof (double));
    DS->EleIS = (double *) malloc(DS->NumEle * sizeof (double));
    DS->EleISmax = (double *) malloc(DS->NumEle * sizeof (double));
    DS->EleISsnowmax = (double *) malloc(DS->NumEle * sizeof (double));
    DS->EleSnow = (double *) malloc(DS->NumEle * sizeof (double));
    DS->EleSnowGrnd = (double *) malloc(DS->NumEle * sizeof (double));
    DS->EleSnowCanopy = (double *) malloc(DS->NumEle * sizeof (double));
    DS->EleTF = (double *) malloc(DS->NumEle * sizeof (double));
    DS->EleETloss = (double *) malloc(DS->NumEle * sizeof (double));
    DS->EleNetPrep = (double *) malloc(DS->NumEle * sizeof (double));
    DS->infil_mode = (double *) malloc(DS->NumEle * sizeof (double));
    DS->EleTemp = (double *) malloc(DS->NumEle * sizeof (double));
    DS->ElePET = (double *) malloc(DS->NumEle * sizeof (double));

    DS->Total_Surfwater_out = (double *) malloc(DS->NumEle * sizeof (double));
    DS->Total_Unsatwater_out = (double *) malloc(DS->NumEle * sizeof (double));
    DS->Total_Satwater_out = (double *) malloc(DS->NumEle * sizeof (double));
    DS->Total_Saltwater_out = (double *) malloc(DS->NumEle * sizeof (double));
    DS->Courant_surf = (double **) malloc(DS->NumEle * sizeof (double*));
    DS->Courant_sub = (double **) malloc(DS->NumEle * sizeof (double*));
    DS->CrossA_surf = (double **) malloc(DS->NumEle * sizeof (double*));
    DS->CrossA_sub = (double **) malloc(DS->NumEle * sizeof (double*));
    DS->Distance = (double **) malloc(DS->NumEle * sizeof (double*));
    DS->FluxSalt = (double **) malloc(DS->NumEle * sizeof (double*));
    if (CS->Lake_module != 0)
    {
        DS->LakeSurfElev = (double *) malloc(DS->NumLake * sizeof (double));
        DS->LakeGWElev = (double *) malloc(DS->NumLake * sizeof (double));
        DS->LakeInitialElev = (double *) malloc(DS->NumLake * sizeof (double));
        DS->LakeOutletDecline = (double *) malloc(DS->NumLake * sizeof (double));
        DS->LakeBed = (double *) malloc(DS->NumLake * sizeof (double));
        DS->Weir_length = (double *) malloc(DS->NumLake * sizeof (double));
        DS->Orifice_height = (double *) malloc(DS->NumLake * sizeof (double));
        DS->BankElevAdjust = (double *) malloc(DS->NumLake * sizeof (double));
    }

    for (i = 0; i < DS->NumOutlet; i++)
    {
        for (j = 0; j < 3; j++)
        {
            if (DS->Ele[DS->Outlet_location[i][0] - 1].BC[j] == 3 || DS->Ele[DS->Outlet_location[i][0] - 1].BC[j] == 4 || DS->Ele[DS->Outlet_location[i][0] - 1].BC[j] == 5 || DS->Ele[DS->Outlet_location[i][0] - 1].BC[j] == 6)
            {
                DS->Outlet_location[i][j + 1] = 1;
            }
        }
    }

    for (i = 0; i < DS->NumSoil; i++)
    {
        DS->Soil[i].ThetaS = CS->Cal.Porosity * DS->Soil[i].ThetaS;
        DS->Soil[i].ThetaR = CS->Cal.Porosity * DS->Soil[i].ThetaR;
    }


    FILE * Filename[2];
    Filename[0] = fopen("node_x.txt", "w");
    Filename[1] = fopen("node_y.txt", "w");

    for (i = 0; i < DS->NumNode; i++)
    {
        fprintf(Filename[0], "%lf\t", DS->Node[i].x);
        fprintf(Filename[0], "\n");
        fprintf(Filename[1], "%lf\t", DS->Node[i].y);
        fprintf(Filename[1], "\n");
    }

    fflush(Filename[0]);
    fflush(Filename[1]);


    for (i = 0; i < DS->NumEle; i++)
    {


        DS->Ele[i].Nabr_index = (int *) malloc(DS->NumEle * sizeof (int));

        for (t = 0; t < DS->NumEle; t++)
        {
            DS->Ele[i].Nabr_index[t] = -1;
        }


        for (j = 0; j < 3; j++)
        {
            inabr = DS->Ele[i].nabr[j] - 1;
            if (inabr >= 0)
            {
                DS->Ele[i].Nabr_index[inabr] = j;
            }
        }


        DS->FluxSurf[i] = (double *) malloc(3 * sizeof (double));
        DS->FluxSub[i] = (double *) malloc(3 * sizeof (double));
        DS->FluxSalt[i] = (double *) malloc(3 * sizeof (double));
        DS->Seepage[i] = (double *) malloc(3 * sizeof (double));
        DS->EleET[i] = (double *) malloc(6 * sizeof (double));
        DS->Courant_surf[i] = (double *) malloc(3 * sizeof (double));
        DS->Courant_sub[i] = (double *) malloc(3 * sizeof (double));
        DS->CrossA_surf[i] = (double *) malloc(3 * sizeof (double));
        DS->CrossA_sub[i] = (double *) malloc(3 * sizeof (double));
        DS->Distance[i] = (double *) malloc(3 * sizeof (double));

        a_x = DS->Node[DS->Ele[i].node[0] - 1].x;
        b_x = DS->Node[DS->Ele[i].node[1] - 1].x;
        c_x = DS->Node[DS->Ele[i].node[2] - 1].x;
        a_y = DS->Node[DS->Ele[i].node[0] - 1].y;
        b_y = DS->Node[DS->Ele[i].node[1] - 1].y;
        c_y = DS->Node[DS->Ele[i].node[2] - 1].y;

        a_zmin = DS->Node[DS->Ele[i].node[0] - 1].zmin;
        b_zmin = DS->Node[DS->Ele[i].node[1] - 1].zmin;
        c_zmin = DS->Node[DS->Ele[i].node[2] - 1].zmin;
        a_zmax = DS->Node[DS->Ele[i].node[0] - 1].zmax;
        b_zmax = DS->Node[DS->Ele[i].node[1] - 1].zmax;
        c_zmax = DS->Node[DS->Ele[i].node[2] - 1].zmax;

        DS->Ele[i].area = 0.5 * ((b_x - a_x) * (c_y - a_y) - (b_y - a_y) * (c_x - a_x));

        // if ((DS->Ele[i].BC[0] == 3)||(DS->Ele[i].BC[1] == 3)||(DS->Ele[i].BC[2] == 3))
        // {
        // DS->Ele[i].area = 10 * DS->Ele[i].area;
        // }

        if (CS->read_centroid == 0)
        {
            DS->Ele[i].zmax = (a_zmax + b_zmax + c_zmax) / 3.0;
            DS->Ele[i].zmin = (a_zmin + b_zmin + c_zmin) / 3.0;
        }
        DS->Ele[i].edge[0] = pow((b_x - c_x), 2) + pow((b_y - c_y), 2);
        DS->Ele[i].edge[1] = pow((c_x - a_x), 2) + pow((c_y - a_y), 2);
        DS->Ele[i].edge[2] = pow((a_x - b_x), 2) + pow((a_y - b_y), 2);
        //DS->zmax_init[i] = 120;
        /* calculate centroid of triangle */
        DS->Ele[i].x = (a_x + b_x + c_x) / 3.0;
        DS->Ele[i].y = (a_y + b_y + c_y) / 3.0;

        /* calculate circumcenter of triangle */
        /*
         * DS->Ele[i].x = a_x - ((b_y - a_y)*DS->Ele[i].edge[2] -
         * (c_y - a_y)*DS->Ele[i].edge[0])/(4*DS->Ele[i].area);
         * DS->Ele[i].y = a_y + ((b_x - a_x)*DS->Ele[i].edge[2] -
         * (c_x - a_x)*DS->Ele[i].edge[0])/(4*DS->Ele[i].area);
         */
        DS->Ele[i].edge[0] = sqrt(DS->Ele[i].edge[0]);
        DS->Ele[i].edge[1] = sqrt(DS->Ele[i].edge[1]);
        DS->Ele[i].edge[2] = sqrt(DS->Ele[i].edge[2]);
        DS->Ele[i].KsatH = CS->Cal.KsatH * DS->Geol[(DS->Ele[i].geol - 1)].KsatH;
        DS->Ele[i].KsatV = CS->Cal.KsatV * DS->Geol[(DS->Ele[i].geol - 1)].KsatV;
        DS->Ele[i].infKsatV = CS->Cal.infKsatV * DS->Soil[(DS->Ele[i].soil - 1)].KsatV;

        //? ? THIS IS ORIG DS->Ele[i].Porosity = CS->Cal.Porosity * (DS->Soil[(DS->Ele[i].soil - 1)].ThetaS - DS->Soil[(DS->Ele[i].soil - 1)].ThetaR);
        DS->Ele[i].Porosity = (DS->Soil[(DS->Ele[i].soil - 1)].ThetaS - DS->Soil[(DS->Ele[i].soil - 1)].ThetaR);
        //? ? XUAN FOR SHALEHILLS ONLY

        /*
         * Note above porosity statement should be replaced by
         * geologic porosity (in comments below) if the data is
         * available
         */
        // DS->Ele[i].Porosity = CS->Cal.Porosity * (DS->Geol[(DS->Ele[i].geol - 1)].ThetaS - DS->Geol[(DS->Ele[i].geol - 1)].ThetaR);
        if ((DS->Ele[i].Porosity > 1) && (DS->Ele[i].Porosity == 0))
        {
            printf("Warning: Porosity value out of bounds");
            getchar();
        }
        DS->Ele[i].Alpha = CS->Cal.Alpha * DS->Soil[(DS->Ele[i].soil - 1)].Alpha;
        DS->Ele[i].Beta = CS->Cal.Beta * DS->Soil[(DS->Ele[i].soil - 1)].Beta;
        /*
         * Note above van genuchten statement should be replaced by
         * geologic parameters (in comments below) if the data is
         * available
         */
        //DS->Ele[i].Alpha = CS->Cal.Alpha * DS->Geol[(DS->Ele[i].geol - 1)].Alpha;
        //DS->Ele[i].Beta = CS->Cal.Beta * DS->Geol[(DS->Ele[i].geol - 1)].Beta;
        DS->Ele[i].hAreaF = CS->Cal.hAreaF * DS->Soil[(DS->Ele[i].soil - 1)].hAreaF;
        DS->Ele[i].vAreaF = CS->Cal.vAreaF * DS->Geol[(DS->Ele[i].geol - 1)].vAreaF;
        DS->Ele[i].macKsatV = CS->Cal.macKsatV * DS->Soil[(DS->Ele[i].soil - 1)].macKsatV;
        DS->Ele[i].macKsatH = CS->Cal.macKsatH * DS->Geol[(DS->Ele[i].geol - 1)].macKsatH;
        DS->Ele[i].macD = CS->Cal.macD * DS->Geol[DS->Ele[i].geol - 1].macD;
        DS->Ele[i].infD = CS->Cal.infD * DS->Soil[DS->Ele[i].soil - 1].infD;

        DS->Ele[i].RzD = CS->Cal.RzD * DS->LandC[DS->Ele[i].LC - 1].RzD;
        DS->Ele[i].LAImax = DS->LandC[DS->Ele[i].LC - 1].LAImax;
        DS->Ele[i].Rmin = DS->LandC[DS->Ele[i].LC - 1].Rmin;
        DS->Ele[i].Rs_ref = DS->LandC[DS->Ele[i].LC - 1].Rs_ref;
        DS->Ele[i].Albedo = CS->Cal.Albedo * DS->LandC[DS->Ele[i].LC - 1].Albedo;
        if (DS->Ele[i].Albedo > 1)
        {
            printf("Warning: Albedo out of bounds");
            getchar();
        }
        DS->Ele[i].VegFrac = CS->Cal.VegFrac * DS->LandC[DS->Ele[i].LC - 1].VegFrac;
        DS->Ele[i].Rough = CS->Cal.Rough * DS->LandC[DS->Ele[i].LC - 1].Rough;

        DS->Ele[i].windH = DS->windH[DS->Ele[i].WindVel - 1];

        fprintf(elev_output, "%lf\t%lf\n", DS->Ele[i].zmax, DS->Ele[i].zmin);
    }
    fflush(elev_output);
    fclose(elev_output);



    for (i = 0; i < DS->NumPrep; i++)
    {
        for (j = 0; j < DS->TSD_Prep[i].length; j++)
        {
            DS->TSD_Prep[i].TS[j][1] = CS->Cal.Prep * DS->TSD_Prep[i].TS[j][1];
        }
    }
    for (i = 0; i < DS->NumTemp; i++)
    {
        for (j = 0; j < DS->TSD_Temp[i].length; j++)
        {
            DS->TSD_Temp[i].TS[j][1] = CS->Cal.Temp * DS->TSD_Temp[i].TS[j][1];
        }
    }
    /* Memory allocation of print variables */
    for (i = 0; i < 31; i++)
    {

        DS->PrintVar[i] = (double *) calloc(DS->NumEle, sizeof (double));
    }



    /*
     * Debugging artifacts in data created due to coarser resolution of
     * model elements
     */
    if (CS->Debug == 1)
    {
        for (i = 0; i < DS->NumEle; i++)
        {
            /*
             * Correction of Surf Elev (artifacts due to coarse
             * scale discretization). Not needed if there is lake
             * feature.
             */
            tmpBool = 1;
            for (j = 0; j < 3; j++)
            {
                if (DS->Ele[i].nabr[j] > 0)
                {
                    //tempvalue1=DS->Ele[i].BC[j]>-4?DS->Ele[DS->Ele[i].nabr[j]-1].zmax:DS->Riv[-(DS->Ele[i].BC[j]/4)-1].zmax;
                    tempvalue1 = DS->Ele[DS->Ele[i].nabr[j] - 1].zmax;
                    if (DS->Ele[i].zmax - tempvalue1 >= 0)
                    {
                        tmpBool = 0;
                        break;
                    }
                }
            }
            if (tmpBool == 1)
            {
                printf("\n Ele %d is sink ", i + 1);
                /*
                 * Note: Following correction is being
                 * applied for debug==1 case only
                 */
                printf("\tBfore: %lf Corrected using:", DS->Ele[i].zmax);
                tempvalue1 = 10000000;
                for (j = 0; j < 3; j++)
                {
                    if (DS->Ele[i].nabr[j] > 0)
                    {
                        //DS->Ele[i].zmax = (DS->Ele[i].BC[j] > -4 ? DS->Ele[DS->Ele[i].nabr[j] - 1].zmax : DS->Riv[-(DS->Ele[i].BC[j] / 4) - 1].zmax);
                        DS->Ele[i].zmax = DS->Ele[DS->Ele[i].nabr[j] - 1].zmax;
                        tempvalue1 = tempvalue1 > DS->Ele[i].zmax ? DS->Ele[i].zmax : tempvalue1;
                        //printf("(%d)%lf  ", j + 1, (DS->Ele[i].BC[j] > -4 ? DS->Ele[DS->Ele[i].nabr[j] - 1].zmax : DS->Riv[-(DS->Ele[i].BC[j] / 4) - 1].zmax));
                        printf("(%d)%lf  ", j + 1, DS->Ele[DS->Ele[i].nabr[j] - 1].zmax);
                    }
                }
                DS->Ele[i].zmax = tempvalue1;
                printf("=(New)%lf  ", DS->Ele[i].zmax);
            }
        }
        /* Correction of BedRck Elev. Is this needed? */
        printf("\n Do you want to correct Bed Rock Elev too (1[y]/0[n])");
        scanf("%d", &BoolBR);
        if (BoolBR == 1)
        {
            for (i = 0; i < DS->NumEle; i++)
            {
                tmpBool = 1;
                for (j = 0; j < 3; j++)
                {
                    if (DS->Ele[i].nabr[j] > 0)
                    {
                        //tempvalue1 = DS->Ele[i].BC[j] > -4 ? DS->Ele[DS->Ele[i].nabr[j] - 1].zmin : DS->Ele[-(DS->Ele[i].BC[j] / 4) - 1 + DS->NumEle].zmin;
                        tempvalue1 = DS->Ele[DS->Ele[i].nabr[j] - 1].zmin;
                        if (DS->Ele[i].zmin - tempvalue1 >= 0)
                        {
                            tmpBool = 0;
                            break;
                        }
                    }
                }
                if (tmpBool == 1)
                {
                    printf("\n Ele %d is sink ", i + 1);
                    /*
                     * Note: Following correction is
                     * being applied for debug==1 case
                     * only
                     */
                    printf("\tBfore: %lf Corrected using:", DS->Ele[i].zmin);
                    tempvalue1 = 10000000;
                    for (j = 0; j < 3; j++)
                    {
                        if (DS->Ele[i].nabr[j] > 0)
                        {
                            //DS->Ele[i].zmin = (DS->Ele[i].BC[j] > -4 ? DS->Ele[DS->Ele[i].nabr[j] - 1].zmin : DS->Ele[-(DS->Ele[i].BC[j] / 4) - 1 + DS->NumEle].zmin);
                            DS->Ele[i].zmin = (DS->Ele[DS->Ele[i].nabr[j] - 1].zmin);
                            tempvalue1 = tempvalue1 > DS->Ele[i].zmin ? DS->Ele[i].zmin : tempvalue1;
                            printf("(%d)%lf  ", j + 1, DS->Ele[DS->Ele[i].nabr[j] - 1].zmin);
                        }
                    }
                    DS->Ele[i].zmin = tempvalue1;
                    printf("=(New)%lf  ", DS->Ele[i].zmin);
                }
            }
        }
        getchar();
        printf("\nHit any key to see more details");
        // for (i = 0; i < DS->NumRiv; i++) {
        // if (DS->Riv[i].down > 0) {
        // if (DS->Riv[i].zmin < DS->Riv[DS->Riv[i].down - 1].zmin) {
        // BoolR = 1;
        // printf("\n Riv %d is lower than downstream Riv %d by %lf", i + 1, DS->Riv[i].down, DS->Riv[i].zmin - DS->Riv[DS->Riv[i].down - 1].zmin);
        // }
        // }
        // }
        if (BoolR == 1)
        {
            printf("\n\tRiver elevation correction needed");
            getchar();
        }
    }

    for (i = 0; i < DS->NumEle; i++)
    {
        a_x = DS->Node[DS->Ele[i].node[0] - 1].x;
        b_x = DS->Node[DS->Ele[i].node[1] - 1].x;
        c_x = DS->Node[DS->Ele[i].node[2] - 1].x;
        a_y = DS->Node[DS->Ele[i].node[0] - 1].y;
        b_y = DS->Node[DS->Ele[i].node[1] - 1].y;
        c_y = DS->Node[DS->Ele[i].node[2] - 1].y;
        for (j = 0; j < 3; j++)
        {
            /*
             * Note: Assumption here is that the forumulation is
             * circumcenter based
             */


            switch (j)
            {
            case 0:
                distX = (DS->Ele[i].x - 0.5 * (b_x + c_x));
                distY = (DS->Ele[i].y - 0.5 * (b_y + c_y));
                break;
            case 1:
                distX = (DS->Ele[i].x - 0.5 * (c_x + a_x));
                distY = (DS->Ele[i].y - 0.5 * (c_y + a_y));
                break;
            case 2:
                distX = (DS->Ele[i].x - 0.5 * (a_x + b_x));
                distY = (DS->Ele[i].y - 0.5 * (a_y + b_y));
                break;
            }
            // DS->Ele[i].surfH[j] = (DS->Ele[i].nabr[j] > 0) ? (DS->Ele[i].BC[j] > -4 ? (DS->Ele[DS->Ele[i].nabr[j] - 1].zmax) : DS->Riv[-(DS->Ele[i].BC[j] / 4) - 1].zmax) : DS->Ele[i].BC[j] <= -4 ? DS->Riv[-(DS->Ele[i].BC[j] / 4) - 1].zmax : (DS->Ele[i].zmax);
            // DS->Ele[i].surfX[j] = (DS->Ele[i].nabr[j] > 0) ? (DS->Ele[i].BC[j] > -4 ? DS->Ele[DS->Ele[i].nabr[j] - 1].x : DS->Riv[-(DS->Ele[i].BC[j] / 4) - 1].x) : (DS->Ele[i].x - 2 * distX);
            // DS->Ele[i].surfY[j] = DS->Ele[i].nabr[j] > 0 ? (DS->Ele[i].BC[j] > -4 ? DS->Ele[DS->Ele[i].nabr[j] - 1].y : DS->Riv[-(DS->Ele[i].BC[j] / 4) - 1].y) : (DS->Ele[i].y - 2 * distY);
            DS->Ele[i].surfH[j] = (DS->Ele[i].nabr[j] > 0) ? ((DS->Ele[DS->Ele[i].nabr[j] - 1].zmax)) : (DS->Ele[i].zmax);
            DS->Ele[i].surfX[j] = (DS->Ele[i].nabr[j] > 0) ? (DS->Ele[DS->Ele[i].nabr[j] - 1].x) : (DS->Ele[i].x - 2 * distX);
            DS->Ele[i].surfY[j] = DS->Ele[i].nabr[j] > 0 ? (DS->Ele[DS->Ele[i].nabr[j] - 1].y) : (DS->Ele[i].y - 2 * distY);
        }
        DS->Ele[i].dhBYdx = -(DS->Ele[i].surfY[2] * (DS->Ele[i].surfH[1] - DS->Ele[i].surfH[0]) + DS->Ele[i].surfY[1] * (DS->Ele[i].surfH[0] - DS->Ele[i].surfH[2]) + DS->Ele[i].surfY[0] * (DS->Ele[i].surfH[2] - DS->Ele[i].surfH[1])) / (DS->Ele[i].surfX[2] * (DS->Ele[i].surfY[1] - DS->Ele[i].surfY[0]) + DS->Ele[i].surfX[1] * (DS->Ele[i].surfY[0] - DS->Ele[i].surfY[2]) + DS->Ele[i].surfX[0] * (DS->Ele[i].surfY[2] - DS->Ele[i].surfY[1]));
        DS->Ele[i].dhBYdy = -(DS->Ele[i].surfX[2] * (DS->Ele[i].surfH[1] - DS->Ele[i].surfH[0]) + DS->Ele[i].surfX[1] * (DS->Ele[i].surfH[0] - DS->Ele[i].surfH[2]) + DS->Ele[i].surfX[0] * (DS->Ele[i].surfH[2] - DS->Ele[i].surfH[1])) / (DS->Ele[i].surfY[2] * (DS->Ele[i].surfX[1] - DS->Ele[i].surfX[0]) + DS->Ele[i].surfY[1] * (DS->Ele[i].surfX[0] - DS->Ele[i].surfX[2]) + DS->Ele[i].surfY[0] * (DS->Ele[i].surfX[2] - DS->Ele[i].surfX[1]));
    }
    /* initialize state variable */
    /* relax case */
    if (CS->init_type == 0)
    {
        for (i = 0; i < DS->NumEle; i++)
        {
            DS->EleIS[i] = 0;
            DS->EleSnow[i] = 0;
            /* Note Two components can be separately read too */
            DS->EleSnowGrnd[i] = (1 - DS->Ele[i].VegFrac) * DS->EleSnow[i];
            DS->EleSnowCanopy[i] = DS->Ele[i].VegFrac * DS->EleSnow[i];

            DS->DummyY[i] = 0;
            DS->DummyY[i + DS->NumEle] = 0;
            DS->DummyY[i + 2 * DS->NumEle] = (DS->Ele[i].zmax - DS->Ele[i].zmin) - 0.1;
            DS->DummyY[i + 3 * DS->NumEle] = 0;
        }
    }
        /* data initialization mode */
    else if (CS->init_type == 1)
    {
        if (DS->UnsatMode == 1)
        {
        }
        if (DS->UnsatMode == 2)
        {
            for (i = 0; i < DS->NumEle; i++)
            {
                DS->EleIS[i] = DS->Ele_IC[i].interception;
                DS->EleSnow[i] = DS->Ele_IC[i].snow;
                /*
                 * Note Two components can be separately read
                 * too
                 */
                DS->EleSnowGrnd[i] = (1 - DS->Ele[i].VegFrac) * DS->EleSnow[i];
                DS->EleSnowCanopy[i] = DS->Ele[i].VegFrac * DS->EleSnow[i];

                DS->DummyY[i] = DS->Ele_IC[i].surf;
                /* Note: delete 0.1 here */

                DS->DummyY[i + DS->NumEle] = DS->Ele_IC[i].unsat;
                DS->DummyY[i + 2 * DS->NumEle] = DS->Ele_IC[i].sat;
                DS->DummyY[i + 2 * DS->NumEle] = DS->Ele_IC[i].sat;
                /* Note: delete line below for general */
                //NV_Ith_S(CV_Y, i + 2 * DS->NumEle) = 0 * DS->Ele_IC[i].sat + (DS->Ele[i].zmax - DS->Ele[i].zmin) * 0.1;
                // if ((NV_Ith_S(CV_Y, i + DS->NumEle) + NV_Ith_S(CV_Y, i + 2 * DS->NumEle)) >= (DS->Ele[i].zmax - DS->Ele[i].zmin)) {
                // NV_Ith_S(CV_Y, i + DS->NumEle) = ((DS->Ele[i].zmax - DS->Ele[i].zmin) - NV_Ith_S(CV_Y, i + 2 * DS->NumEle)) * 0.98;
                // if (NV_Ith_S(CV_Y, i + DS->NumEle) < 0) {
                // NV_Ith_S(CV_Y, i + DS->NumEle) = 0;
                // }
                // }				
            }
        }
    }
        /* hot start mode */
    else
    {
        fn = (char *) malloc((strlen(filename) + 6) * sizeof (char));
        strcpy(fn, filename);
        init_file = fopen(strcat(fn, ".init"), "r");
        free(fn);
        if (init_file == NULL)
        {
            printf("\n  Fatal Error: %s.init is in use or does not exist!\n", filename);
            exit(1);
        }
        else
        {
            for (i = 0; i < DS->NumEle; i++)
            {
                fscanf(init_file, "%lf %lf %lf %lf %lf %lf", &DS->EleIS[i], &DS->EleSnow[i], &tempvalue1, &tempvalue2, &tempvalue3, &tempvalue4);
                DS->EleSnowGrnd[i] = (1 - DS->Ele[i].VegFrac) * DS->EleSnow[i];
                DS->EleSnowCanopy[i] = DS->Ele[i].VegFrac * DS->EleSnow[i];

                DS->DummyY[i] = tempvalue1;

                DS->DummyY[i + DS->NumEle] = tempvalue2;
                DS->DummyY[i + 2 * DS->NumEle] = tempvalue3;
                DS->DummyY[i + 3 * DS->NumEle] = tempvalue4;
            }
        }
        fclose(init_file);
    }

    printf("done.\n");


    for (i = 0; i < DS->NumOutlet; i++)
    {
        printf("	%d	", DS->NumOutlet);
        printf("	%d	", DS->Outlet_location[i][0]);
        printf("	%d	", DS->Outlet_location[i][1]);
        printf("	%d	", DS->Outlet_location[i][2]);
        printf("	%d	", DS->Outlet_location[i][3]);
        printf("	\n	");
    }
}


