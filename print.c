/*******************************************************************************
 * File        : print.c	                                               *
 * Version     : Nov, 2007 (2.0)                                               *
 * Function    : print out model results output files                          *
 * Developer of PIHM2.0:        Mukesh Kumar (muk139@psu.edu)                  *
 * Developer of PIHM1.0:        Yizhong Qu   (quyizhong@gmail.com)             *
 *-----------------------------------------------------------------------------*
 *                                                                             *
 *                                                                             *
 *..............MODIFICATIONS/ADDITIONS in PIHM 2.0............................*
 * a) This file is downgraded from Version 1.0, as no ancillary results are    *
 *    being output			                                       *
 * b) Only state variables and flux to/in/accross river and its bed are being  *
 *    output							               *
 * c) Addition of Average Function to output average variables at regular time *
 *    intervals								       *
 *******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/stat.h>

#include "nvector_serial.h"
#include "sundials_types.h"
#include "cvode.h"
#include "cvode_dense.h"

#include "pihm.h"
#include "lake.h"

/* Temporal average of State vectors */
void
avgResults_NV(FILE * fpin, double * tmpVarCal, Model_Data tmpDS, int tmpIntv, int tmpNumObj, double tmpt, int tmpInitObj)
{
    int j;
    int TmpIntv;

    TmpIntv = tmpIntv;

    for (j = 0; j < tmpNumObj; j++)
    {
#ifdef MEAN
        tmpVarCal[j] = tmpVarCal[j] + tmpDS->DummyY[j + tmpInitObj];
#endif
#ifdef REALTIME
        tmpVarCal[j] = tmpDS->DummyY[j + tmpInitObj];
#endif
    }
    if (((int) tmpt % tmpIntv) == 0)
    {
        fprintf(fpin, "%lf\t", tmpt);
        for (j = 0; j < tmpNumObj; j++)
        {

#ifdef MEAN
            fprintf(fpin, "%lf\t", tmpVarCal[j] / TmpIntv);
#endif
#ifdef REALTIME
            fprintf(fpin, "%lf\t", tmpVarCal[j]);
#endif

            tmpVarCal[j] = 0;
        }
        fprintf(fpin, "\n");
        fflush(fpin);
    }
}

void
avgResults_NV_lake(FILE * fpin, double * tmpVarCal, int tmpIntv, int tmpNumObj, double tmpt, double * StateVariable)
{
    int j;
    int TmpIntv;

    TmpIntv = tmpIntv;

    for (j = 0; j < tmpNumObj; j++)
    {
#ifdef MEAN
        tmpVarCal[j] = tmpVarCal[j] + StateVariable[j];
#endif
#ifdef REALTIME
        tmpVarCal[j] = StateVariable[j];
#endif
    }
    if (((int) tmpt % tmpIntv) == 0)
    {
        fprintf(fpin, "%lf\t", tmpt);
        for (j = 0; j < tmpNumObj; j++)
        {

#ifdef MEAN
            fprintf(fpin, "%lf\t", tmpVarCal[j] / TmpIntv);
#endif
#ifdef REALTIME
            fprintf(fpin, "%lf\t", tmpVarCal[j]);
#endif

            tmpVarCal[j] = 0;
        }
        fprintf(fpin, "\n");
        fflush(fpin);
    }
}

/* Temporal average of Derived states */
void
avgResults_MD(FILE * fpin, double * tmpVarCal, Model_Data tmpDS, int tmpIntv, int tmpNumObj, double tmpt, int tmpFC)
{
    int j;
    int TmpIntv;

    TmpIntv = tmpIntv;

    switch (tmpFC)
    {
    case 2:
        for (j = 0; j < tmpNumObj; j++)
        {
            tmpVarCal[j] = tmpVarCal[j] + tmpDS->EleET[j][0];
        }
        break;
    case 3:
        for (j = 0; j < tmpNumObj; j++)
        {
            tmpVarCal[j] = tmpVarCal[j] + tmpDS->EleET[j][1];
        }
        break;
    case 4:
        for (j = 0; j < tmpNumObj; j++)
        {
            tmpVarCal[j] = tmpVarCal[j] + tmpDS->EleET[j][2];
        }
        break;
    case 5:
        for (j = 0; j < tmpNumObj; j++)
        {
            tmpVarCal[j] = tmpVarCal[j] + tmpDS->EleIS[j];
        }
        break;
    case 6:
        for (j = 0; j < tmpNumObj; j++)
        {
            tmpVarCal[j] = tmpVarCal[j] + (tmpDS->EleSnowCanopy[j] + tmpDS->EleSnowGrnd[j]);
        }
        break;
    case 8:
        for (j = 0; j < tmpNumObj; j++)
        {
            tmpVarCal[j] = tmpVarCal[j] + tmpDS->Recharge[j];
        }
        break;
    case 9:
        for (j = 0; j < tmpNumObj; j++)
        {
            tmpVarCal[j] = tmpVarCal[j] + tmpDS->EleViR[j];
        }
        break;

    case 10:
        for (j = 0; j < tmpNumObj; j++)
        {
            tmpVarCal[j] = tmpVarCal[j]+(tmpDS->ElePrep[j]);
        }
        break;
    case 11:
        for (j = 0; j < tmpNumObj; j++)
        {
            tmpVarCal[j] = tmpVarCal[j] + tmpDS->FluxSurf[j][0];
        }
        break;
    case 12:
        for (j = 0; j < tmpNumObj; j++)
        {
            tmpVarCal[j] = tmpVarCal[j] + tmpDS->FluxSurf[j][1];
        }
        break;
    case 13:
        for (j = 0; j < tmpNumObj; j++)
        {
            tmpVarCal[j] = tmpVarCal[j] + tmpDS->FluxSurf[j][2];
        }
        break;
    case 14:
        for (j = 0; j < tmpNumObj; j++)
        {
            tmpVarCal[j] = tmpVarCal[j] + tmpDS->FluxSub[j][0];
        }
        break;
    case 15:
        for (j = 0; j < tmpNumObj; j++)
        {
            tmpVarCal[j] = tmpVarCal[j] + tmpDS->FluxSub[j][1];
        }
        break;
    case 16:
        for (j = 0; j < tmpNumObj; j++)
        {
            tmpVarCal[j] = tmpVarCal[j] + tmpDS->FluxSub[j][2];
        }
        break;
    case 17:
        for (j = 0; j < tmpNumObj; j++)
        {
            tmpVarCal[j] = tmpVarCal[j]+(tmpDS->EleNetPrep[j]);
        }
        break;
    case 18:
        for (j = 0; j < tmpNumObj; j++)
        {
            tmpVarCal[j] = tmpVarCal[j] + tmpDS->Courant_surf[j][0];
        }
        break;

    case 19:
        for (j = 0; j < tmpNumObj; j++)
        {
            tmpVarCal[j] = tmpVarCal[j] + tmpDS->Courant_surf[j][1];
        }
        break;
    case 20:
        for (j = 0; j < tmpNumObj; j++)
        {
            tmpVarCal[j] = tmpVarCal[j] + tmpDS->Courant_surf[j][2];
        }
        break;
    case 21:
        for (j = 0; j < tmpNumObj; j++)
        {
            tmpVarCal[j] = tmpVarCal[j]+(tmpDS->Courant_sub[j][0]);
        }
        break;
    case 22:
        for (j = 0; j < tmpNumObj; j++)
        {
            tmpVarCal[j] = tmpVarCal[j]+(tmpDS->Courant_sub[j][1]);
        }
        break;
    case 23:
        for (j = 0; j < tmpNumObj; j++)
        {
            tmpVarCal[j] = tmpVarCal[j]+(tmpDS->Courant_sub[j][2]);
        }
        break;
    case 24:
        for (j = 0; j < tmpNumObj; j++)
        {
            tmpVarCal[j] = tmpVarCal[j] + tmpDS->EleET[j][3];
        }
        break;
    case 25:
        for (j = 0; j < tmpNumObj; j++)
        {
            tmpVarCal[j] = tmpVarCal[j] + tmpDS->EleET[j][4];
        }
        break;
    case 26:
        for (j = 0; j < tmpNumObj; j++)
        {
            tmpVarCal[j] = tmpVarCal[j] + tmpDS->EleET[j][5];
        }
        break;
    case 29:
        for (j = 0; j < tmpNumObj; j++)
        {
            tmpVarCal[j] = tmpVarCal[j] + tmpDS->EleTemp[j];
        }
        break;
    case 30:
        for (j = 0; j < tmpNumObj; j++)
        {
            tmpVarCal[j] = tmpVarCal[j] + tmpDS->ElePET[j];
        }
        break;
    default:
        break;
    }
    if (((int) tmpt % tmpIntv) == 0)
    {
        fprintf(fpin, "%lf\t", tmpt);
        for (j = 0; j < tmpNumObj; j++)
        {
            fprintf(fpin, "%lf\t", tmpVarCal[j] / TmpIntv);
            tmpVarCal[j] = 0;
        }
        fprintf(fpin, "\n");
        fflush(fpin);
    }
}

void
avgResults_MD_Lake(FILE * fpin, double * tmpVarCal, Model_Data MD, LakeData LD, int tmpIntv, int tmpNumObj, double tmpt, int tmpFC)
{
    int j, k, t;
    int TmpIntv;
    int element1, element2;

    TmpIntv = tmpIntv;

    switch (tmpFC)
    {

    case 2:
        for (j = 0; j < tmpNumObj / 2; j++)
        {
            for (k = 0; k < LD->Lake[j].NumBankEle; k++)
            {
                if (LD->Lake[j].StreamIndex[k] == -1)
                {
                    if (LD->FluxSurf[j][k] > 0)
                        tmpVarCal[j * 2] = tmpVarCal[j * 2]+(LD->FluxSurf[j][k]);
                    else
                        tmpVarCal[j * 2 + 1] = tmpVarCal[j * 2 + 1]+(LD->FluxSurf[j][k]);
                }
            }
        }
        break;

    case 3:
        for (j = 0; j < tmpNumObj / 2; j++)
        {
            for (k = 0; k < LD->Lake[j].NumBankEle; k++)
            {
                if (LD->FluxSub[j][k] > 0)
                    tmpVarCal[j * 2] = tmpVarCal[j * 2]+(LD->FluxSub[j][k]);
                else
                    tmpVarCal[j * 2 + 1] = tmpVarCal[j * 2 + 1]+(LD->FluxSub[j][k]);
            }
        }
        break;

    case 4:
        for (j = 0; j < tmpNumObj / 2; j++)
        {
            for (k = 0; k < LD->Lake[j].NumBankEle; k++)
            {
                if (LD->Lake[j].StreamIndex[k] == 1)
                {
                    if (LD->FluxSurf[j][k] > 0)
                        tmpVarCal[j * 2] = tmpVarCal[j * 2]+(LD->FluxSurf[j][k]);
                    else
                        tmpVarCal[j * 2 + 1] = tmpVarCal[j * 2 + 1]+(LD->FluxSurf[j][k]);
                }
            }
        }
        break;

    case 5:
        for (j = 0; j < tmpNumObj; j++)
        {
            tmpVarCal[j] = tmpVarCal[j] + (LD->ET[j]);
        }
        break;
    case 6:
        for (j = 0; j < tmpNumObj; j++)
        {
            tmpVarCal[j] = tmpVarCal[j]+(LD->Infil[j]);
        }
        break;
    case 7:
        for (j = 0; j < tmpNumObj; j++)
        {
            tmpVarCal[j] = tmpVarCal[j]+(LD->Lake[j].Temp);
        }
        break;
    case 8:
        for (j = 0; j < tmpNumObj; j++)
        {
            tmpVarCal[j] = tmpVarCal[j]+(LD->Lake[j].Precip);
        }
        break;
    case 9:
        for (j = 0; j < tmpNumObj; j++)
        {
            for (t = 0; t < LD->NumRiverInterest; t++)//river flow with users interests
            {
                element1 = LD->RiverInterestInformation[j][t][0];
                element2 = LD->RiverInterestInformation[j][t][1];
                for (k = 0; k < LD->Lake[j].NumBankEle; k++)
                {
                    if (LD->Lake[j].BankEle[k] == element1 || LD->Lake[j].BankEle[k] == element2)
                    {
                        if (LD->FluxSurf[j][k] >= 0)
                        {
                            tmpVarCal[j * (LD->NumRiverInterest + 4) + t] = tmpVarCal[j * (LD->NumRiverInterest + 4) + t]+(LD->FluxSurf[j][k]);
                        }
                    }
                }
            }
            for (k = 0; k < LD->Lake[j].NumBankEle; k++)//all stream flow
            {
                if (LD->Lake[j].StreamIndex[k] == 1)
                {
                    if (LD->FluxSurf[j][k] > 0)
                        tmpVarCal[j * (LD->NumRiverInterest + 4) + LD->NumRiverInterest] = tmpVarCal[j * (LD->NumRiverInterest + 4) + LD->NumRiverInterest]+(LD->FluxSurf[j][k]);
                    else
                        tmpVarCal[j * (LD->NumRiverInterest + 4) + LD->NumRiverInterest + 3] = tmpVarCal[j * (LD->NumRiverInterest + 4) + LD->NumRiverInterest + 3]+(LD->FluxSurf[j][k]);
                }
                else //all overland flow
                {
                    if (LD->FluxSurf[j][k] > 0)
                        tmpVarCal[j * (LD->NumRiverInterest + 4) + LD->NumRiverInterest + 1] = tmpVarCal[j * (LD->NumRiverInterest + 4) + LD->NumRiverInterest + 1]+(LD->FluxSurf[j][k]);
                    else
                        tmpVarCal[j * (LD->NumRiverInterest + 4) + LD->NumRiverInterest + 3] = tmpVarCal[j * (LD->NumRiverInterest + 4) + LD->NumRiverInterest + 3]+(LD->FluxSurf[j][k]);
                }

                if (LD->FluxSub[j][k] > 0) //subsurface flow
                    tmpVarCal[j * (LD->NumRiverInterest + 4) + LD->NumRiverInterest + 2] = tmpVarCal[j * (LD->NumRiverInterest + 4) + LD->NumRiverInterest + 2]+(LD->FluxSub[j][k]);
                else
                    tmpVarCal[j * (LD->NumRiverInterest + 4) + LD->NumRiverInterest + 3] = tmpVarCal[j * (LD->NumRiverInterest + 4) + LD->NumRiverInterest + 3]+(LD->FluxSub[j][k]);
            }
        }
        break;
    default:
        break;
    }
    if (tmpFC == 9)
    {
        if (((int) tmpt % tmpIntv) == 0)
        {
            fprintf(fpin, "%lf\t", tmpt);
            for (j = 0; j < (LD->NumRiverInterest + 4) * tmpNumObj; j++)
            {
                fprintf(fpin, "%lf\t", tmpVarCal[j] / TmpIntv);
                tmpVarCal[j] = 0;
            }
            fprintf(fpin, "\n");
            fflush(fpin);
        }
    }
    else
    {
        if (((int) tmpt % tmpIntv) == 0)
        {
            fprintf(fpin, "%lf\t", tmpt);
            for (j = 0; j < tmpNumObj; j++)
            {
                fprintf(fpin, "%lf\t", tmpVarCal[j] / TmpIntv);
                tmpVarCal[j] = 0;
            }
            fprintf(fpin, "\n");
            fflush(fpin);
        }
    }

}

/* print individual states */
void
PrintData(FILE ** outp, FILE ** outp_Lake, Control_Data * cD, Model_Data DS, LakeData LD, double t)
{
    if (cD->gwDInt > 0)
    {
        avgResults_NV(outp[0], DS->PrintVar[0], DS, cD->gwDInt, DS->NumEle, t, 2 * DS->NumEle); //GW
    }

    if (cD->surfDInt > 0)
    {
        avgResults_NV(outp[1], DS->PrintVar[1], DS, cD->surfDInt, DS->NumEle, t, 0 * DS->NumEle); //Surf
        avgResults_NV(outp[28], DS->PrintVar[28], DS, cD->surfDInt, DS->NumEle, t, 3 * DS->NumEle); //SaltW
    }
    if (cD->etInt > 0)
    {
        avgResults_MD(outp[2], DS->PrintVar[2], DS, cD->etInt, DS->NumEle, t, 2); //et0
        avgResults_MD(outp[3], DS->PrintVar[3], DS, cD->etInt, DS->NumEle, t, 3); //et1
        avgResults_MD(outp[4], DS->PrintVar[4], DS, cD->etInt, DS->NumEle, t, 4); //et2
        avgResults_MD(outp[24], DS->PrintVar[24], DS, cD->etInt, DS->NumEle, t, 24); //et3
        avgResults_MD(outp[25], DS->PrintVar[25], DS, cD->etInt, DS->NumEle, t, 25); //et4
        avgResults_MD(outp[26], DS->PrintVar[26], DS, cD->etInt, DS->NumEle, t, 26); //et5
        avgResults_MD(outp[29], DS->PrintVar[29], DS, cD->etInt, DS->NumEle, t, 29); //Temperature, EleTemp
        avgResults_MD(outp[30], DS->PrintVar[30], DS, cD->etInt, DS->NumEle, t, 30); //Potential ET, ElePET
    }
    if (cD->IsDInt > 0)
    {
        avgResults_MD(outp[5], DS->PrintVar[5], DS, cD->IsDInt, DS->NumEle, t, 5); //Interception
    }
    if (cD->snowDInt > 0)
    {
        avgResults_MD(outp[6], DS->PrintVar[6], DS, cD->snowDInt, DS->NumEle, t, 6); //Snow
    }
    if (cD->usDInt > 0)
    {
        avgResults_NV(outp[7], DS->PrintVar[7], DS, cD->usDInt, DS->NumEle, t, 1 * DS->NumEle); //Unsat
    }
    if (cD->RechInt > 0)
    {
        avgResults_MD(outp[8], DS->PrintVar[8], DS, cD->RechInt, DS->NumEle, t, 8); //Recharge
        avgResults_MD(outp[9], DS->PrintVar[9], DS, cD->RechInt, DS->NumEle, t, 9); //Infiltration
    }
    if (cD->prepDInt > 0)
    {
        avgResults_MD(outp[10], DS->PrintVar[10], DS, cD->prepDInt, DS->NumEle, t, 10); //Precipitation
        avgResults_MD(outp[17], DS->PrintVar[17], DS, cD->prepDInt, DS->NumEle, t, 17); //Trough Fall precipitation
    }
    if (cD->surfDInt > 0)
    {
        avgResults_MD(outp[11], DS->PrintVar[11], DS, cD->surfDInt, DS->NumEle, t, 11); //Flux1
        avgResults_MD(outp[12], DS->PrintVar[12], DS, cD->surfDInt, DS->NumEle, t, 12); //Flux2
        avgResults_MD(outp[13], DS->PrintVar[13], DS, cD->surfDInt, DS->NumEle, t, 13); //Flux3
        avgResults_MD(outp[14], DS->PrintVar[14], DS, cD->surfDInt, DS->NumEle, t, 14); //SubFlux1
        avgResults_MD(outp[15], DS->PrintVar[15], DS, cD->surfDInt, DS->NumEle, t, 15); //SubFlux2
        avgResults_MD(outp[16], DS->PrintVar[16], DS, cD->surfDInt, DS->NumEle, t, 16); //SubFlux3
        avgResults_MD(outp[18], DS->PrintVar[18], DS, cD->surfDInt, DS->NumEle, t, 18); //Courant_surf1
        avgResults_MD(outp[19], DS->PrintVar[19], DS, cD->surfDInt, DS->NumEle, t, 19); //Courant_surf2
        avgResults_MD(outp[20], DS->PrintVar[20], DS, cD->surfDInt, DS->NumEle, t, 20); //Courant_surf3
        avgResults_MD(outp[21], DS->PrintVar[21], DS, cD->surfDInt, DS->NumEle, t, 21); //Courant_sub1
        avgResults_MD(outp[22], DS->PrintVar[22], DS, cD->surfDInt, DS->NumEle, t, 22); //Courant_sub2
        avgResults_MD(outp[23], DS->PrintVar[23], DS, cD->surfDInt, DS->NumEle, t, 23); //Courant_sub3
    }

    if (cD->surfDInt > 0 && cD->Lake_module != 0)
    {
        avgResults_NV_lake(outp_Lake[0], LD->PrintVar[0], cD->surfDInt, LD->NumLake, t, LD->LakeSurfDepth); //LakeSurf
        avgResults_NV_lake(outp_Lake[1], LD->PrintVar[1], cD->gwDInt, LD->NumLake, t, LD->LakeGW); //LakeGW
        avgResults_MD_Lake(outp_Lake[2], LD->PrintVar[2], DS, LD, cD->surfDInt, 2 * LD->NumLake, t, 2); //lakeFluxSurf
        avgResults_MD_Lake(outp_Lake[3], LD->PrintVar[3], DS, LD, cD->surfDInt, 2 * LD->NumLake, t, 3); //lakeFluxSub
        avgResults_MD_Lake(outp_Lake[4], LD->PrintVar[4], DS, LD, cD->surfDInt, 2 * LD->NumLake, t, 4); //lakeFluxStream
        avgResults_MD_Lake(outp_Lake[5], LD->PrintVar[5], DS, LD, cD->surfDInt, LD->NumLake, t, 5); //lakeET
        avgResults_MD_Lake(outp_Lake[6], LD->PrintVar[6], DS, LD, cD->surfDInt, LD->NumLake, t, 6); //lakeInfil
        avgResults_MD_Lake(outp_Lake[7], LD->PrintVar[7], DS, LD, cD->surfDInt, LD->NumLake, t, 7); //lakeTemp
        avgResults_MD_Lake(outp_Lake[8], LD->PrintVar[8], DS, LD, cD->surfDInt, LD->NumLake, t, 8); //lakePrecip
        avgResults_MD_Lake(outp_Lake[9], LD->PrintVar[9], DS, LD, cD->surfDInt, LD->NumLake, t, 9); //lakeFlux
    }
}

void
copy_inputs(char *directory, char *filename, Control_Data * cD)
{
    FILE *orig_file, *desti_file;
    char *temp_char, *temp_string;
    char ch;
    int flag, pos;
    temp_string = (char *) malloc(1024 * sizeof (char));
    sprintf(temp_string, "%sinput/", directory);
    flag = mkdir(temp_string, 0755);
    /********* pihm.h**********/
    orig_file = fopen("pihm.h", "r");
    desti_file = fopen(strcat(temp_string, "pihm.h"), "w");
    fseek(orig_file, 0L, SEEK_END); // file pointer at end of file
    pos = ftell(orig_file);
    fseek(orig_file, 0L, SEEK_SET); // file pointer set at start
    while (pos--)
    {
        ch = fgetc(orig_file); // copying file character by character
        fputc(ch, desti_file);
    }
    fclose(orig_file);
    fclose(desti_file);
    free(temp_string);
    /********* pihm.c**********/
    temp_string = (char *) malloc(1024 * sizeof (char));
    sprintf(temp_string, "%sinput/", directory);
    orig_file = fopen("pihm.c", "r");
    desti_file = fopen(strcat(temp_string, "pihm.c"), "w");
    fseek(orig_file, 0L, SEEK_END); // file pointer at end of file
    pos = ftell(orig_file);
    fseek(orig_file, 0L, SEEK_SET); // file pointer set at start
    while (pos--)
    {
        ch = fgetc(orig_file); // copying file character by character
        fputc(ch, desti_file);
    }
    fclose(orig_file);
    fclose(desti_file);
    free(temp_string);
    /********* f.c**********/
    temp_string = (char *) malloc(1024 * sizeof (char));
    sprintf(temp_string, "%sinput/", directory);
    orig_file = fopen("f.c", "r");
    desti_file = fopen(strcat(temp_string, "f.c"), "w");
    fseek(orig_file, 0L, SEEK_END); // file pointer at end of file
    pos = ftell(orig_file);
    fseek(orig_file, 0L, SEEK_SET); // file pointer set at start
    while (pos--)
    {
        ch = fgetc(orig_file); // copying file character by character
        fputc(ch, desti_file);
    }
    fclose(orig_file);
    fclose(desti_file);
    free(temp_string);

    /********* initialize.c**********/
    temp_string = (char *) malloc(1024 * sizeof (char));
    sprintf(temp_string, "%sinput/", directory);
    orig_file = fopen("initialize.c", "r");
    desti_file = fopen(strcat(temp_string, "initialize.c"), "w");
    fseek(orig_file, 0L, SEEK_END); // file pointer at end of file
    pos = ftell(orig_file);
    fseek(orig_file, 0L, SEEK_SET); // file pointer set at start
    while (pos--)
    {
        ch = fgetc(orig_file); // copying file character by character
        fputc(ch, desti_file);
    }
    fclose(orig_file);
    fclose(desti_file);
    free(temp_string);
    /********* is_sm_et.c**********/
    temp_string = (char *) malloc(1024 * sizeof (char));
    sprintf(temp_string, "%sinput/", directory);
    orig_file = fopen("is_sm_et.c", "r");
    desti_file = fopen(strcat(temp_string, "is_sm_et.c"), "w");
    fseek(orig_file, 0L, SEEK_END); // file pointer at end of file
    pos = ftell(orig_file);
    fseek(orig_file, 0L, SEEK_SET); // file pointer set at start
    while (pos--)
    {
        ch = fgetc(orig_file); // copying file character by character
        fputc(ch, desti_file);
    }
    fclose(orig_file);
    fclose(desti_file);
    free(temp_string);

    /********* print.c**********/
    temp_string = (char *) malloc(1024 * sizeof (char));
    sprintf(temp_string, "%sinput/", directory);
    orig_file = fopen("print.c", "r");
    desti_file = fopen(strcat(temp_string, "print.c"), "w");
    fseek(orig_file, 0L, SEEK_END); // file pointer at end of file
    pos = ftell(orig_file);
    fseek(orig_file, 0L, SEEK_SET); // file pointer set at start
    while (pos--)
    {
        ch = fgetc(orig_file); // copying file character by character
        fputc(ch, desti_file);
    }
    fclose(orig_file);
    fclose(desti_file);
    free(temp_string);
    /********* read_alloc.c**********/
    temp_string = (char *) malloc(1024 * sizeof (char));
    sprintf(temp_string, "%sinput/", directory);
    orig_file = fopen("read_alloc.c", "r");
    desti_file = fopen(strcat(temp_string, "read_alloc.c"), "w");
    fseek(orig_file, 0L, SEEK_END); // file pointer at end of file
    pos = ftell(orig_file);
    fseek(orig_file, 0L, SEEK_SET); // file pointer set at start
    while (pos--)
    {
        ch = fgetc(orig_file); // copying file character by character
        fputc(ch, desti_file);
    }
    fclose(orig_file);
    fclose(desti_file);
    free(temp_string);
    /********* update.c**********/
    temp_string = (char *) malloc(1024 * sizeof (char));
    sprintf(temp_string, "%sinput/", directory);
    orig_file = fopen("update.c", "r");
    desti_file = fopen(strcat(temp_string, "update.c"), "w");
    fseek(orig_file, 0L, SEEK_END); // file pointer at end of file
    pos = ftell(orig_file);
    fseek(orig_file, 0L, SEEK_SET); // file pointer set at start
    while (pos--)
    {
        ch = fgetc(orig_file); // copying file character by character
        fputc(ch, desti_file);
    }
    fclose(orig_file);
    fclose(desti_file);
    free(temp_string);
    /********* Makefile**********/
    temp_string = (char *) malloc(1024 * sizeof (char));
    sprintf(temp_string, "%sinput/", directory);
    orig_file = fopen("Makefile", "r");
    desti_file = fopen(strcat(temp_string, "Makefile"), "w");
    fseek(orig_file, 0L, SEEK_END); // file pointer at end of file
    pos = ftell(orig_file);
    fseek(orig_file, 0L, SEEK_SET); // file pointer set at start
    while (pos--)
    {
        ch = fgetc(orig_file); // copying file character by character
        fputc(ch, desti_file);
    }
    fclose(orig_file);
    fclose(desti_file);
    free(temp_string);
    /********* LakeDynamic.c**********/
    if (cD->Lake_module != 0)
    {
        temp_string = (char *) malloc(1024 * sizeof (char));
        sprintf(temp_string, "%sinput/", directory);
        orig_file = fopen("LakeDynamic.c", "r");
        desti_file = fopen(strcat(temp_string, "LakeDynamic.c"), "w");
        fseek(orig_file, 0L, SEEK_END); // file pointer at end of file
        pos = ftell(orig_file);
        fseek(orig_file, 0L, SEEK_SET); // file pointer set at start
        while (pos--)
        {
            ch = fgetc(orig_file); // copying file character by character
            fputc(ch, desti_file);
        }
        fclose(orig_file);
        fclose(desti_file);
        free(temp_string);
    }

    /********* initialize_lake.c**********/
    if (cD->Lake_module != 0)
    {
        temp_string = (char *) malloc(1024 * sizeof (char));
        sprintf(temp_string, "%sinput/", directory);
        orig_file = fopen("initialize_lake.c", "r");
        desti_file = fopen(strcat(temp_string, "initialize_lake.c"), "w");
        fseek(orig_file, 0L, SEEK_END); // file pointer at end of file
        pos = ftell(orig_file);
        fseek(orig_file, 0L, SEEK_SET); // file pointer set at start
        while (pos--)
        {
            ch = fgetc(orig_file); // copying file character by character
            fputc(ch, desti_file);
        }
        fclose(orig_file);
        fclose(desti_file);
        free(temp_string);
    }
    /********* read_alloc_lake.c**********/
    if (cD->Lake_module != 0)
    {
        temp_string = (char *) malloc(1024 * sizeof (char));
        sprintf(temp_string, "%sinput/", directory);
        orig_file = fopen("read_alloc_lake.c", "r");
        desti_file = fopen(strcat(temp_string, "read_alloc_lake.c"), "w");
        fseek(orig_file, 0L, SEEK_END); // file pointer at end of file
        pos = ftell(orig_file);
        fseek(orig_file, 0L, SEEK_SET); // file pointer set at start
        while (pos--)
        {
            ch = fgetc(orig_file); // copying file character by character
            fputc(ch, desti_file);
        }
        fclose(orig_file);
        fclose(desti_file);
        free(temp_string);
    }
    /********* Lake_land_exchange.c**********/
    if (cD->Lake_module != 0)
    {
        temp_string = (char *) malloc(1024 * sizeof (char));
        sprintf(temp_string, "%sinput/", directory);
        orig_file = fopen("Lake_land_exchange.c", "r");
        desti_file = fopen(strcat(temp_string, "Lake_land_exchange.c"), "w");
        fseek(orig_file, 0L, SEEK_END); // file pointer at end of file
        pos = ftell(orig_file);
        fseek(orig_file, 0L, SEEK_SET); // file pointer set at start
        while (pos--)
        {
            ch = fgetc(orig_file); // copying file character by character
            fputc(ch, desti_file);
        }
        fclose(orig_file);
        fclose(desti_file);
        free(temp_string);
    }
    /********* lake.h**********/
    if (cD->Lake_module != 0)
    {
        temp_string = (char *) malloc(1024 * sizeof (char));
        sprintf(temp_string, "%sinput/", directory);
        orig_file = fopen("lake.h", "r");
        desti_file = fopen(strcat(temp_string, "lake.h"), "w");
        fseek(orig_file, 0L, SEEK_END); // file pointer at end of file
        pos = ftell(orig_file);
        fseek(orig_file, 0L, SEEK_SET); // file pointer set at start
        while (pos--)
        {
            ch = fgetc(orig_file); // copying file character by character
            fputc(ch, desti_file);
        }
        fclose(orig_file);
        fclose(desti_file);
        free(temp_string);
    }
    /********* .calib**********/
    temp_string = (char *) malloc(1024 * sizeof (char));
    sprintf(temp_string, "%s.calib", filename);
    orig_file = fopen(temp_string, "r");
    free(temp_string);
    temp_string = (char *) malloc(1024 * sizeof (char));
    sprintf(temp_string, "%sinput/%s.calib", directory, filename);
    desti_file = fopen(temp_string, "w");
    fseek(orig_file, 0L, SEEK_END); // file pointer at end of file
    pos = ftell(orig_file);
    fseek(orig_file, 0L, SEEK_SET); // file pointer set at start
    while (pos--)
    {
        ch = fgetc(orig_file); // copying file character by character
        fputc(ch, desti_file);
    }
    fclose(orig_file);
    fclose(desti_file);
    free(temp_string);
    // /********* .forc**********/
     temp_string = (char *) malloc(1024 * sizeof(char));
     sprintf(temp_string,"%s.forc",filename);
     orig_file = fopen(temp_string,"r");
     free(temp_string);
     temp_string = (char *) malloc(1024 * sizeof(char));
     sprintf(temp_string,"%sinput/%s.forc",directory,filename);
     desti_file = fopen(temp_string, "w");  
     fseek(orig_file, 0L, SEEK_END); // file pointer at end of file
     pos = ftell(orig_file);
     fseek(orig_file, 0L, SEEK_SET); // file pointer set at start
     while (pos--)
     {
     ch = fgetc(orig_file);  // copying file character by character
     fputc(ch, desti_file);
     }    
     fclose(orig_file);
     fclose(desti_file);
     free(temp_string);
    /********* .geol**********/
    temp_string = (char *) malloc(1024 * sizeof (char));
    sprintf(temp_string, "%s.geol", filename);
    orig_file = fopen(temp_string, "r");
    free(temp_string);
    temp_string = (char *) malloc(1024 * sizeof (char));
    sprintf(temp_string, "%sinput/%s.geol", directory, filename);
    desti_file = fopen(temp_string, "w");
    fseek(orig_file, 0L, SEEK_END); // file pointer at end of file
    pos = ftell(orig_file);
    fseek(orig_file, 0L, SEEK_SET); // file pointer set at start
    while (pos--)
    {
        ch = fgetc(orig_file); // copying file character by character
        fputc(ch, desti_file);
    }
    fclose(orig_file);
    fclose(desti_file);
    free(temp_string);
    /********* .ibc**********/
    temp_string = (char *) malloc(1024 * sizeof (char));
    sprintf(temp_string, "%s.ibc", filename);
    orig_file = fopen(temp_string, "r");
    free(temp_string);
    temp_string = (char *) malloc(1024 * sizeof (char));
    sprintf(temp_string, "%sinput/%s.ibc", directory, filename);
    desti_file = fopen(temp_string, "w");
    fseek(orig_file, 0L, SEEK_END); // file pointer at end of file
    pos = ftell(orig_file);
    fseek(orig_file, 0L, SEEK_SET); // file pointer set at start
    while (pos--)
    {
        ch = fgetc(orig_file); // copying file character by character
        fputc(ch, desti_file);
    }
    fclose(orig_file);
    fclose(desti_file);
    free(temp_string);
    /********* .init**********/
    temp_string = (char *) malloc(1024 * sizeof (char));
    sprintf(temp_string, "%s.init", filename);
    orig_file = fopen(temp_string, "r");
    if (orig_file != NULL)
    {
        free(temp_string);
        temp_string = (char *) malloc(1024 * sizeof (char));
        sprintf(temp_string, "%sinput/%s.init", directory, filename);
        desti_file = fopen(temp_string, "w");
        fseek(orig_file, 0L, SEEK_END); // file pointer at end of file
        pos = ftell(orig_file);
        fseek(orig_file, 0L, SEEK_SET); // file pointer set at start
        while (pos--)
        {
            ch = fgetc(orig_file); // copying file character by character
            fputc(ch, desti_file);
        }
        fclose(orig_file);
        fclose(desti_file);
        free(temp_string);
    }

    /********* .lc**********/
    temp_string = (char *) malloc(1024 * sizeof (char));
    sprintf(temp_string, "%s.lc", filename);
    orig_file = fopen(temp_string, "r");
    free(temp_string);
    temp_string = (char *) malloc(1024 * sizeof (char));
    sprintf(temp_string, "%sinput/%s.lc", directory, filename);
    desti_file = fopen(temp_string, "w");
    fseek(orig_file, 0L, SEEK_END); // file pointer at end of file
    pos = ftell(orig_file);
    fseek(orig_file, 0L, SEEK_SET); // file pointer set at start
    while (pos--)
    {
        ch = fgetc(orig_file); // copying file character by character
        fputc(ch, desti_file);
    }
    fclose(orig_file);
    fclose(desti_file);
    free(temp_string);
    /********* .mesh**********/
    temp_string = (char *) malloc(1024 * sizeof (char));
    sprintf(temp_string, "%s.mesh", filename);
    orig_file = fopen(temp_string, "r");
    free(temp_string);
    temp_string = (char *) malloc(1024 * sizeof (char));
    sprintf(temp_string, "%sinput/%s.mesh", directory, filename);
    desti_file = fopen(temp_string, "w");
    fseek(orig_file, 0L, SEEK_END); // file pointer at end of file
    pos = ftell(orig_file);
    fseek(orig_file, 0L, SEEK_SET); // file pointer set at start
    while (pos--)
    {
        ch = fgetc(orig_file); // copying file character by character
        fputc(ch, desti_file);
    }
    fclose(orig_file);
    fclose(desti_file);
    free(temp_string);
    /********* .att**********/
    temp_string = (char *) malloc(1024 * sizeof (char));
    sprintf(temp_string, "%s.att", filename);
    orig_file = fopen(temp_string, "r");
    free(temp_string);
    temp_string = (char *) malloc(1024 * sizeof (char));
    sprintf(temp_string, "%sinput/%s.att", directory, filename);
    desti_file = fopen(temp_string, "w");
    fseek(orig_file, 0L, SEEK_END); // file pointer at end of file
    pos = ftell(orig_file);
    fseek(orig_file, 0L, SEEK_SET); // file pointer set at start
    while (pos--)
    {
        ch = fgetc(orig_file); // copying file character by character
        fputc(ch, desti_file);
    }
    fclose(orig_file);
    fclose(desti_file);
    free(temp_string);
    /********* .initelev**********/

    temp_string = (char *) malloc(1024 * sizeof (char));
    sprintf(temp_string, "%s.initelev", filename);
    orig_file = fopen(temp_string, "r");
    if (orig_file != NULL)
    {
        free(temp_string);
        temp_string = (char *) malloc(1024 * sizeof (char));
        sprintf(temp_string, "%sinput/%s.initelev", directory, filename);
        desti_file = fopen(temp_string, "w");
        fseek(orig_file, 0L, SEEK_END); // file pointer at end of file
        pos = ftell(orig_file);
        fseek(orig_file, 0L, SEEK_SET); // file pointer set at start
        while (pos--)
        {
            ch = fgetc(orig_file); // copying file character by character
            fputc(ch, desti_file);
        }
        fclose(orig_file);
        fclose(desti_file);
        free(temp_string);
    }

    /********* .newelev**********/
    temp_string = (char *) malloc(1024 * sizeof (char));
    sprintf(temp_string, "%s.newelev", filename);
    orig_file = fopen(temp_string, "r");
    if (orig_file != NULL)
    {
        free(temp_string);
        temp_string = (char *) malloc(1024 * sizeof (char));
        sprintf(temp_string, "%sinput/%s.newelev", directory, filename);
        desti_file = fopen(temp_string, "w");
        fseek(orig_file, 0L, SEEK_END); // file pointer at end of file
        pos = ftell(orig_file);
        fseek(orig_file, 0L, SEEK_SET); // file pointer set at start
        while (pos--)
        {
            ch = fgetc(orig_file); // copying file character by character
            fputc(ch, desti_file);
        }
        fclose(orig_file);
        fclose(desti_file);
        free(temp_string);
    }
    /********* .soil**********/
    temp_string = (char *) malloc(1024 * sizeof (char));
    sprintf(temp_string, "%s.soil", filename);
    orig_file = fopen(temp_string, "r");
    free(temp_string);
    temp_string = (char *) malloc(1024 * sizeof (char));
    sprintf(temp_string, "%sinput/%s.soil", directory, filename);
    desti_file = fopen(temp_string, "w");
    fseek(orig_file, 0L, SEEK_END); // file pointer at end of file
    pos = ftell(orig_file);
    fseek(orig_file, 0L, SEEK_SET); // file pointer set at start
    while (pos--)
    {
        ch = fgetc(orig_file); // copying file character by character
        fputc(ch, desti_file);
    }
    fclose(orig_file);
    fclose(desti_file);
    free(temp_string);

    /********* projectName.txt**********/
    temp_string = (char *) malloc(1024 * sizeof (char));
    sprintf(temp_string, "%sinput/", directory);
    orig_file = fopen("projectName.txt", "r");
    desti_file = fopen(strcat(temp_string, "projectName.txt"), "w");
    fseek(orig_file, 0L, SEEK_END); // file pointer at end of file
    pos = ftell(orig_file);
    fseek(orig_file, 0L, SEEK_SET); // file pointer set at start
    while (pos--)
    {
        ch = fgetc(orig_file); // copying file character by character
        fputc(ch, desti_file);
    }
    fclose(orig_file);
    fclose(desti_file);
    free(temp_string);

    /********* qsub.txt**********/
    temp_string = (char *) malloc(1024 * sizeof (char));
    sprintf(temp_string, "%sinput/", directory);
    orig_file = fopen("qsub.txt", "r");
    desti_file = fopen(strcat(temp_string, "qsub.txt"), "w");
    fseek(orig_file, 0L, SEEK_END); // file pointer at end of file
    pos = ftell(orig_file);
    fseek(orig_file, 0L, SEEK_SET); // file pointer set at start
    while (pos--)
    {
        ch = fgetc(orig_file); // copying file character by character
        fputc(ch, desti_file);
    }
    fclose(orig_file);
    fclose(desti_file);
    free(temp_string);

    /********* .para**********/
    temp_string = (char *) malloc(1024 * sizeof (char));
    sprintf(temp_string, "%s.para", filename);
    orig_file = fopen(temp_string, "r");
    free(temp_string);
    temp_string = (char *) malloc(1024 * sizeof (char));
    sprintf(temp_string, "%sinput/%s.para", directory, filename);
    desti_file = fopen(temp_string, "w");
    fseek(orig_file, 0L, SEEK_END); // file pointer at end of file
    pos = ftell(orig_file);
    fseek(orig_file, 0L, SEEK_SET); // file pointer set at start
    while (pos--)
    {
        ch = fgetc(orig_file); // copying file character by character
        fputc(ch, desti_file);
    }
    fclose(orig_file);
    fclose(desti_file);
    free(temp_string);
    /********* .bid**********/
    temp_string = (char *) malloc(1024 * sizeof (char));
    sprintf(temp_string, "%s.bid", filename);
    orig_file = fopen(temp_string, "r");
    free(temp_string);
    temp_string = (char *) malloc(1024 * sizeof (char));
    sprintf(temp_string, "%sinput/%s.bid", directory, filename);
    desti_file = fopen(temp_string, "w");
    fseek(orig_file, 0L, SEEK_END); // file pointer at end of file
    pos = ftell(orig_file);
    fseek(orig_file, 0L, SEEK_SET); // file pointer set at start
    while (pos--)
    {
        ch = fgetc(orig_file); // copying file character by character
        fputc(ch, desti_file);
    }
    fclose(orig_file);
    fclose(desti_file);
    free(temp_string);
    /********* .tide**********/
    temp_string = (char *) malloc(1024 * sizeof (char));
    sprintf(temp_string, "%s.tide", filename);
    orig_file = fopen(temp_string, "r");
    free(temp_string);
    temp_string = (char *) malloc(1024 * sizeof (char));
    sprintf(temp_string, "%sinput/%s.tide", directory, filename);
    desti_file = fopen(temp_string, "w");
    fseek(orig_file, 0L, SEEK_END); // file pointer at end of file
    pos = ftell(orig_file);
    fseek(orig_file, 0L, SEEK_SET); // file pointer set at start
    while (pos--)
    {
        ch = fgetc(orig_file); // copying file character by character
        fputc(ch, desti_file);
    }
    fclose(orig_file);
    fclose(desti_file);
    free(temp_string);
    
    /********* .lakeatt**********/
    if (cD->Lake_module != 0)
    {
        temp_string = (char *) malloc(1024 * sizeof (char));
        sprintf(temp_string, "%s.lakeatt", filename);
        orig_file = fopen(temp_string, "r");
        free(temp_string);
        temp_string = (char *) malloc(1024 * sizeof (char));
        sprintf(temp_string, "%sinput/%s.lakeatt", directory, filename);
        desti_file = fopen(temp_string, "w");
        fseek(orig_file, 0L, SEEK_END); // file pointer at end of file
        pos = ftell(orig_file);
        fseek(orig_file, 0L, SEEK_SET); // file pointer set at start
        while (pos--)
        {
            ch = fgetc(orig_file); // copying file character by character
            fputc(ch, desti_file);
        }
        fclose(orig_file);
        fclose(desti_file);
        free(temp_string);
    }
    /********* .lakebathy**********/
    if (cD->Lake_module != 0)
    {
        temp_string = (char *) malloc(1024 * sizeof (char));
        sprintf(temp_string, "%s.lakebathy", filename);
        orig_file = fopen(temp_string, "r");
        free(temp_string);
        temp_string = (char *) malloc(1024 * sizeof (char));
        sprintf(temp_string, "%sinput/%s.lakebathy", directory, filename);
        desti_file = fopen(temp_string, "w");
        fseek(orig_file, 0L, SEEK_END); // file pointer at end of file
        pos = ftell(orig_file);
        fseek(orig_file, 0L, SEEK_SET); // file pointer set at start
        while (pos--)
        {
            ch = fgetc(orig_file); // copying file character by character
            fputc(ch, desti_file);
        }
        fclose(orig_file);
        fclose(desti_file);
        free(temp_string);
    }
    /********* .lakebathy**********/
    if (cD->Lake_module != 0)
    {
        temp_string = (char *) malloc(1024 * sizeof (char));
        sprintf(temp_string, "%s.lakebathy", filename);
        orig_file = fopen(temp_string, "r");
        free(temp_string);
        temp_string = (char *) malloc(1024 * sizeof (char));
        sprintf(temp_string, "%sinput/%s.lakebathy", directory, filename);
        desti_file = fopen(temp_string, "w");
        fseek(orig_file, 0L, SEEK_END); // file pointer at end of file
        pos = ftell(orig_file);
        fseek(orig_file, 0L, SEEK_SET); // file pointer set at start
        while (pos--)
        {
            ch = fgetc(orig_file); // copying file character by character
            fputc(ch, desti_file);
        }
        fclose(orig_file);
        fclose(desti_file);
        free(temp_string);
    }
    /********* .lakegeom**********/
    if (cD->Lake_module != 0)
    {
        temp_string = (char *) malloc(1024 * sizeof (char));
        sprintf(temp_string, "%s.lakegeom", filename);
        orig_file = fopen(temp_string, "r");
        free(temp_string);
        temp_string = (char *) malloc(1024 * sizeof (char));
        sprintf(temp_string, "%sinput/%s.lakegeom", directory, filename);
        desti_file = fopen(temp_string, "w");
        fseek(orig_file, 0L, SEEK_END); // file pointer at end of file
        pos = ftell(orig_file);
        fseek(orig_file, 0L, SEEK_SET); // file pointer set at start
        while (pos--)
        {
            ch = fgetc(orig_file); // copying file character by character
            fputc(ch, desti_file);
        }
        fclose(orig_file);
        fclose(desti_file);
        free(temp_string);
    }
    /********* .lakesoil**********/
    if (cD->Lake_module != 0)
    {
        temp_string = (char *) malloc(1024 * sizeof (char));
        sprintf(temp_string, "%s.lakesoil", filename);
        orig_file = fopen(temp_string, "r");
        free(temp_string);
        temp_string = (char *) malloc(1024 * sizeof (char));
        sprintf(temp_string, "%sinput/%s.lakesoil", directory, filename);
        desti_file = fopen(temp_string, "w");
        fseek(orig_file, 0L, SEEK_END); // file pointer at end of file
        pos = ftell(orig_file);
        fseek(orig_file, 0L, SEEK_SET); // file pointer set at start
        while (pos--)
        {
            ch = fgetc(orig_file); // copying file character by character
            fputc(ch, desti_file);
        }
        fclose(orig_file);
        fclose(desti_file);
        free(temp_string);
    }
    /********* .riverInterest**********/
    if (cD->Lake_module != 0)
    {
        temp_string = (char *) malloc(1024 * sizeof (char));
        sprintf(temp_string, "%s.riverInterest", filename);
        orig_file = fopen(temp_string, "r");
        free(temp_string);
        temp_string = (char *) malloc(1024 * sizeof (char));
        sprintf(temp_string, "%sinput/%s.riverInterest", directory, filename);
        desti_file = fopen(temp_string, "w");
        fseek(orig_file, 0L, SEEK_END); // file pointer at end of file
        pos = ftell(orig_file);
        fseek(orig_file, 0L, SEEK_SET); // file pointer set at start
        while (pos--)
        {
            ch = fgetc(orig_file); // copying file character by character
            fputc(ch, desti_file);
        }
        fclose(orig_file);
        fclose(desti_file);
        free(temp_string);
    }
    /********* .lakeinit**********/
    if (cD->Lake_module != 0)
    {
        temp_string = (char *) malloc(1024 * sizeof (char));
        sprintf(temp_string, "%s.lakeinit", filename);
        orig_file = fopen(temp_string, "r");
        free(temp_string);
        temp_string = (char *) malloc(1024 * sizeof (char));
        sprintf(temp_string, "%sinput/%s.lakeinit", directory, filename);
        desti_file = fopen(temp_string, "w");
        fseek(orig_file, 0L, SEEK_END); // file pointer at end of file
        pos = ftell(orig_file);
        fseek(orig_file, 0L, SEEK_SET); // file pointer set at start
        while (pos--)
        {
            ch = fgetc(orig_file); // copying file character by character
            fputc(ch, desti_file);
        }
        fclose(orig_file);
        fclose(desti_file);
        free(temp_string);
    }
}




