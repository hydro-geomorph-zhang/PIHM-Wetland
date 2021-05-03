/* File        : pihm.c                                                        *
 * Version     : Nov, 2007 (2.0)                                               *
 * Developer of PIHM2.0:        Mukesh Kumar (muk139@psu.edu)                  *
 * Developer of PIHM1.0:        Yizhong Qu   (quyizhong@gmail.com)             *
 *-----------------------------------------------------------------------------*
 *                                                                             *
 *..............MODIFICATIONS/ADDITIONS in PIHM 2.0............................*
 * a) All modifications in physical process representations  in this version   *
 *    are listed as header in f.c and is_sm_et.c.     			       *
 * b) All addition/modifications in variable and structure definition/declarat-*
 *    -ion are listed as header in read_alloc.c and initialize.c	       *
 * c) 3 new input files have been added for geology, landcover and calibration *
 *    data								       *
 * d) Ported to Sundials 2.1.0                                                 *
 *-----------------------------------------------------------------------------*
 *                                                                             *
 * PIHM is an integrated finite volume based hydrologic model. It simulates    *
 * channel routing, overland flow, groundwater flow, macropore based infiltra- *
 * tion and stormflow, throughfall, evaporation from overlandflow-subsurface-  *
 * canopy, transpiration and  snowmelt by full coupling of processes.          *
 * It uses semi-discrete finite volume approach to discretize PDEs into ODEs,  *
 * and henceforth solving the global system of ODEs using CVODE. Global ODEs   *
 * are created in f.c. Any modifications in the process equations has to be    *
 * performed in f.c
 *                                                                             *
 *-----------------------------------------------------------------------------*
 * For questions or comments, please contact                                   *
 *      --> Mukesh Kumar (muk139@psu.edu)                                      *
 *      --> Prof. Chris Duffy (cxd11@psu.edu)                                  *
 * This code is free for research purpose only.                                *
 * Please provide relevant references if you use this code in your research work*
 *-----------------------------------------------------------------------------*
 *									       *
 * DEVELOPMENT RELATED REFERENCES:					       *
 * PIHM2.0:								       *
 *	a) Kumar, M., 2008, "Development and Implementation of a Multiscale,   *
 *	Multiprocess Hydrologic Model". PhD Thesis, Penn State University      *
 *	b) Kumar, M, G.Bhatt & C.Duffy, "Coupling of Data and Processes in     *
 *	Mesoscale Watershed", Advances in Water Resources (submitted)          *
 * PIHM1.0:								       *
 *	a) Qu, Y., 2005, "An Integrated hydrologic model for multiproces       *
 *	simulation using semi-discrete finite volume approach".PhD Thesis, PSU *
 *	b) Qu, Y. & C. Duffy, 2007, "A semidiscrete finite volume formulation  *
 *	for multiprocess watershed simulation". Water Resources Research       *
 *-----------------------------------------------------------------------------*
 * LICENSE:
 *******************************************************************************/



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <sys/stat.h>

/* SUNDIAL Header Files */
#include "sundials_types.h" /* double, integertype, booleantype
				 * defination */
#include "cvode.h"  /* CVODE header file                             */
#include "cvode_spgmr.h" /* CVSPGMR linear header file                    */
#include "sundials_smalldense.h"/* use generic DENSE linear solver for
				 * "small"   */
#include "nvector_serial.h" /* contains the definition of type N_Vector      */
#include "sundials_math.h" /* contains UnitRoundoff, RSqrt, SQR
				 * functions   */
#include "cvode_dense.h" /* CVDENSE header file                           */
#include "sundials_dense.h" /* generic dense solver header file              */
#include "pihm.h"  /* Data Model and Variable Declarations     */
#include "lake.h"
#define UNIT_C 1440  /* Unit Conversions */


/* Function Declarations */

void initialize(char *, Model_Data, Control_Data *, FILE*);
void is_sm_et(double, double, Model_Data);

/* Function to calculate right hand side of ODE systems */
int f_explicit(double, Model_Data);
void read_alloc(char *, Model_Data, Control_Data *); /* Variable definition */
void update(double, void *);
void PrintData(FILE **, FILE **, Control_Data *, Model_Data, LakeData, double);
void FreeData(Model_Data, Control_Data *, LakeData);
void copy_inputs(char *, char *, Control_Data *);
double Interpolation(TSD * Data, double t);
double avgY(double diff, double yi, double yinabr);
double effKV(double ksatFunc, double gradY, double macKV, double KV, double areaF);
double effKH(int mp, double tmpY, double aqDepth, double MacD, double MacKsatH, double areaF, double ksatH);
void Courant(Model_Data);
void HydroSurfLateral(Model_Data, int, int, int, double);
void HydroSubLateral(Model_Data, int, int, int, double);
void ET(Model_Data, int, double);
void Infiltration(Model_Data, int);
double* initiation(Model_Data MD, double *, double);
void Adjust_water(Model_Data);
double* Calculate_DY(Model_Data, double *);



void read_alloc_lake(char *filename, LakeData LD, Model_Data DS);
void initialize_lake(char *filename, LakeData LD, Model_Data DS, Control_Data * CS, FILE *, FILE *);
void Lake_land_exchange(Model_Data DS, LakeData LD);
void LakeDynamic(double t, LakeData LD, Model_Data MD);

/* Main Function */
int
main(int argc, char *argv[])
{
    Model_Data mData; /* Model Data                */
    Control_Data cData; /* Solver Control Data       */
    LakeData LD;
    void *cvode_mem; /* Model Data Pointer        */
    int flag; /* flag to test return value */
    FILE * Ofile[31]; /* Output file     */
    char *ofn[31];
    FILE * Ofile_Lake[10]; /* Output file     */
    char *ofn_Lake[10];
    FILE *init1, *init2, *para_file, *lakeinit1, *lakeinit2, *Lakeboundary, *Lakeboundarynode;
    char *init_char1, *init_char2, *init_char3, *init_char4, *para_file_char, *directory, *path, *temp_path;
    char *Lakeboundary_char, *Lakeboundarynode_char, *LakeReference_char, *Elev_output_char;
    char *init_file;
    int IsExist = 0;
    int OutputID = 1;
    double total_t, P;
    FILE *iproj; /* Project File */
    int N; /* Problem size              */
    int i, k; /* loop index                */
    double NextPtr, StepSize, total_NextPtr, total_LE_NextPtr; /* stress period & step size */
    clock_t start, finish, diff_t; /* system clock at points    */
    char *filename;
    double T;
    int q;
    int EleID;
    time_t current_time;
    FILE *centroidfile, *LakeReference, *Elev_output;
    char *tempchar;
    int t_count = 5;

    /* Project Input Name */
    if (argc != 2)
    {
        iproj = fopen("projectName.txt", "r");
        if (iproj == NULL)
        {
            printf("\t\nUsage ./pihm project_name");
            printf("\t\n         OR              ");
            printf("\t\nUsage ./pihm, and have a file in the current directory named projectName.txt with the project name in it");
            exit(0);
        }
        else
        {
            filename = (char *) malloc(15 * sizeof (char));
            fscanf(iproj, "%s", filename);
            fclose(iproj);
        }
    }
    else
    {
        /* get user specified file name in command line */
        filename = (char *) malloc(strlen(argv[1]) * sizeof (char));
        strcpy(filename, argv[1]);
    }
    /* Open Output Files */
    directory = (char *) malloc(1024 * sizeof (char));
    path = (char *) malloc(1024 * sizeof (char));
    temp_path = (char *) malloc(1024 * sizeof (char));

    strcpy(path, "output/");
    i = mkdir(path, 0755);
    sprintf(temp_path, "%s%d/", path, OutputID);
    i = mkdir(temp_path, 0755);
    while (i == -1)
    {
        free(temp_path);
        temp_path = (char *) malloc(1024 * sizeof (char));
        OutputID += 1;
        sprintf(temp_path, "%s%d/", path, OutputID);
        i = mkdir(temp_path, 0755);
    }

    strcpy(directory, temp_path);

    /* allocate memory for model data structure */
    mData = (Model_Data) malloc(sizeof *mData);
    read_alloc(filename, mData, &cData);

    copy_inputs(directory, filename, &cData);

    ofn[0] = (char *) malloc((strlen(filename) + 1024) * sizeof (char));
    strcpy(ofn[0], directory);
    strcat(ofn[0], filename);
    Ofile[0] = fopen(strcat(ofn[0], ".GW"), "w");
    ofn[1] = (char *) malloc((strlen(filename) + 1024) * sizeof (char));
    strcpy(ofn[1], directory);
    strcat(ofn[1], filename);
    Ofile[1] = fopen(strcat(ofn[1], ".surf"), "w");
    ofn[2] = (char *) malloc((strlen(filename) + 1024) * sizeof (char));
    strcpy(ofn[2], directory);
    strcat(ofn[2], filename);
    Ofile[2] = fopen(strcat(ofn[2], ".et0"), "w");
    ofn[3] = (char *) malloc((strlen(filename) + 1024) * sizeof (char));
    strcpy(ofn[3], directory);
    strcat(ofn[3], filename);
    Ofile[3] = fopen(strcat(ofn[3], ".et1"), "w");
    ofn[4] = (char *) malloc((strlen(filename) + 1024) * sizeof (char));
    strcpy(ofn[4], directory);
    strcat(ofn[4], filename);
    Ofile[4] = fopen(strcat(ofn[4], ".et2"), "w");
    ofn[5] = (char *) malloc((strlen(filename) + 1024) * sizeof (char));
    strcpy(ofn[5], directory);
    strcat(ofn[5], filename);
    Ofile[5] = fopen(strcat(ofn[5], ".is"), "w");
    ofn[6] = (char *) malloc((strlen(filename) + 1024) * sizeof (char));
    strcpy(ofn[6], directory);
    strcat(ofn[6], filename);
    Ofile[6] = fopen(strcat(ofn[6], ".snow"), "w");
    ofn[7] = (char *) malloc((strlen(filename) + 1024) * sizeof (char));
    strcpy(ofn[7], directory);
    strcat(ofn[7], filename);
    Ofile[7] = fopen(strcat(ofn[7], ".unsat"), "w");
    ofn[8] = (char *) malloc((strlen(filename) + 1024) * sizeof (char));
    strcpy(ofn[8], directory);
    strcat(ofn[8], filename);
    Ofile[8] = fopen(strcat(ofn[8], ".Rech"), "w");
    ofn[9] = (char *) malloc((strlen(filename) + 1024) * sizeof (char));
    strcpy(ofn[9], directory);
    strcat(ofn[9], filename);
    Ofile[9] = fopen(strcat(ofn[9], ".infil"), "w");
    ofn[10] = (char *) malloc((strlen(filename) + 1024) * sizeof (char));
    strcpy(ofn[10], directory);
    strcat(ofn[10], filename);
    Ofile[10] = fopen(strcat(ofn[10], ".prep"), "w");
    ofn[11] = (char *) malloc((strlen(filename) + 1024) * sizeof (char));
    strcpy(ofn[11], directory);
    strcat(ofn[11], filename);
    Ofile[11] = fopen(strcat(ofn[11], ".Flux1"), "w");
    ofn[12] = (char *) malloc((strlen(filename) + 1024) * sizeof (char));
    strcpy(ofn[12], directory);
    strcat(ofn[12], filename);
    Ofile[12] = fopen(strcat(ofn[12], ".Flux2"), "w");
    ofn[13] = (char *) malloc((strlen(filename) + 1024) * sizeof (char));
    strcpy(ofn[13], directory);
    strcat(ofn[13], filename);
    Ofile[13] = fopen(strcat(ofn[13], ".Flux3"), "w");
    ofn[14] = (char *) malloc((strlen(filename) + 1024) * sizeof (char));
    strcpy(ofn[14], directory);
    strcat(ofn[14], filename);
    Ofile[14] = fopen(strcat(ofn[14], ".SubFlux1"), "w");
    ofn[15] = (char *) malloc((strlen(filename) + 1024) * sizeof (char));
    strcpy(ofn[15], directory);
    strcat(ofn[15], filename);
    Ofile[15] = fopen(strcat(ofn[15], ".SubFlux2"), "w");
    ofn[16] = (char *) malloc((strlen(filename) + 1024) * sizeof (char));
    strcpy(ofn[16], directory);
    strcat(ofn[16], filename);
    Ofile[16] = fopen(strcat(ofn[16], ".SubFlux3"), "w");
    ofn[17] = (char *) malloc((strlen(filename) + 1024) * sizeof (char));
    strcpy(ofn[17], directory);
    strcat(ofn[17], filename);
    Ofile[17] = fopen(strcat(ofn[17], ".NetPrep"), "w");

    ofn[18] = (char *) malloc((strlen(filename) + 1024) * sizeof (char));
    strcpy(ofn[18], directory);
    strcat(ofn[18], filename);
    Ofile[18] = fopen(strcat(ofn[18], ".Courant_surf1"), "w");
    ofn[19] = (char *) malloc((strlen(filename) + 1024) * sizeof (char));
    strcpy(ofn[19], directory);
    strcat(ofn[19], filename);
    Ofile[19] = fopen(strcat(ofn[19], ".Courant_surf2"), "w");
    ofn[20] = (char *) malloc((strlen(filename) + 1024) * sizeof (char));
    strcpy(ofn[20], directory);
    strcat(ofn[20], filename);
    Ofile[20] = fopen(strcat(ofn[20], ".Courant_surf3"), "w");
    ofn[21] = (char *) malloc((strlen(filename) + 1024) * sizeof (char));
    strcpy(ofn[21], directory);
    strcat(ofn[21], filename);
    Ofile[21] = fopen(strcat(ofn[21], ".Courant_sub1"), "w");
    ofn[22] = (char *) malloc((strlen(filename) + 1024) * sizeof (char));
    strcpy(ofn[22], directory);
    strcat(ofn[22], filename);
    Ofile[22] = fopen(strcat(ofn[22], ".Courant_sub2"), "w");
    ofn[23] = (char *) malloc((strlen(filename) + 1024) * sizeof (char));
    strcpy(ofn[23], directory);
    strcat(ofn[23], filename);
    Ofile[23] = fopen(strcat(ofn[23], ".Courant_sub3"), "w");
    ofn[24] = (char *) malloc((strlen(filename) + 1024) * sizeof (char));
    strcpy(ofn[24], directory);
    strcat(ofn[24], filename);
    Ofile[24] = fopen(strcat(ofn[24], ".et3"), "w");
    ofn[25] = (char *) malloc((strlen(filename) + 1024) * sizeof (char));
    strcpy(ofn[25], directory);
    strcat(ofn[25], filename);
    Ofile[25] = fopen(strcat(ofn[25], ".et4"), "w");
    ofn[26] = (char *) malloc((strlen(filename) + 1024) * sizeof (char));
    strcpy(ofn[26], directory);
    strcat(ofn[26], filename);
    Ofile[26] = fopen(strcat(ofn[26], ".et5"), "w");
    ofn[27] = (char *) malloc((strlen(filename) + 1024) * sizeof (char));
    strcpy(ofn[27], directory);
    strcat(ofn[27], filename);
    Ofile[27] = fopen(strcat(ofn[27], ".Rn"), "w");
    ofn[28] = (char *) malloc((strlen(filename) + 1024) * sizeof (char));
    strcpy(ofn[28], directory);
    strcat(ofn[28], filename);
    Ofile[28] = fopen(strcat(ofn[28], ".SaltW"), "w");
    ofn[29] = (char *) malloc((strlen(filename) + 1024) * sizeof (char));
    strcpy(ofn[29], directory);
    strcat(ofn[29], filename);
    Ofile[29] = fopen(strcat(ofn[29], ".Temperature"), "w");
    ofn[30] = (char *) malloc((strlen(filename) + 1024) * sizeof (char));
    strcpy(ofn[30], directory);
    strcat(ofn[30], filename);
    Ofile[30] = fopen(strcat(ofn[30], ".PET"), "w");



    if (cData.Lake_module != 0)
    {
        ofn_Lake[0] = (char *) malloc((strlen(filename) + 1024) * sizeof (char));
        strcpy(ofn_Lake[0], directory);
        strcat(ofn_Lake[0], filename);
        Ofile_Lake[0] = fopen(strcat(ofn_Lake[0], ".lakeSurf"), "w");
        ofn_Lake[1] = (char *) malloc((strlen(filename) + 1024) * sizeof (char));
        strcpy(ofn_Lake[1], directory);
        strcat(ofn_Lake[1], filename);
        Ofile_Lake[1] = fopen(strcat(ofn_Lake[1], ".lakeGW"), "w");
        ofn_Lake[2] = (char *) malloc((strlen(filename) + 1024) * sizeof (char));
        strcpy(ofn_Lake[2], directory);
        strcat(ofn_Lake[2], filename);
        Ofile_Lake[2] = fopen(strcat(ofn_Lake[2], ".lakeFluxSurf"), "w");
        ofn_Lake[3] = (char *) malloc((strlen(filename) + 1024) * sizeof (char));
        strcpy(ofn_Lake[3], directory);
        strcat(ofn_Lake[3], filename);
        Ofile_Lake[3] = fopen(strcat(ofn_Lake[3], ".lakeFluxSub"), "w");
        ofn_Lake[4] = (char *) malloc((strlen(filename) + 1024) * sizeof (char));
        strcpy(ofn_Lake[4], directory);
        strcat(ofn_Lake[4], filename);
        Ofile_Lake[4] = fopen(strcat(ofn_Lake[4], ".lakeFluxStream"), "w");
        ofn_Lake[5] = (char *) malloc((strlen(filename) + 1024) * sizeof (char));
        strcpy(ofn_Lake[5], directory);
        strcat(ofn_Lake[5], filename);
        Ofile_Lake[5] = fopen(strcat(ofn_Lake[5], ".lakeET"), "w");
        ofn_Lake[6] = (char *) malloc((strlen(filename) + 1024) * sizeof (char));
        strcpy(ofn_Lake[6], directory);
        strcat(ofn_Lake[6], filename);
        Ofile_Lake[6] = fopen(strcat(ofn_Lake[6], ".lakeInfil"), "w");
        ofn_Lake[7] = (char *) malloc((strlen(filename) + 1024) * sizeof (char));
        strcpy(ofn_Lake[7], directory);
        strcat(ofn_Lake[7], filename);
        Ofile_Lake[7] = fopen(strcat(ofn_Lake[7], ".lakeTemp"), "w");
        ofn_Lake[8] = (char *) malloc((strlen(filename) + 1024) * sizeof (char));
        strcpy(ofn_Lake[8], directory);
        strcat(ofn_Lake[8], filename);
        Ofile_Lake[8] = fopen(strcat(ofn_Lake[8], ".lakePrecip"), "w");
        ofn_Lake[9] = (char *) malloc((strlen(filename) + 1024) * sizeof (char));
        strcpy(ofn_Lake[9], directory);
        strcat(ofn_Lake[9], filename);
        Ofile_Lake[9] = fopen(strcat(ofn_Lake[9], ".lakeFlux"), "w");
    }


    init_char1 = (char *) malloc((strlen(filename) + 1024) * sizeof (char));
    init_char2 = (char *) malloc((strlen(filename) + 1024) * sizeof (char));
    strcpy(init_char1, directory);
    strcat(init_char1, filename);
    strcat(init_char1, ".init");
    // strcat(init_char2, "/input/");
    strcpy(init_char2, filename);
    strcat(init_char2, ".init");

    if (cData.Lake_module != 0)
    {
        init_char3 = (char *) malloc((strlen(filename) + 1024) * sizeof (char));
        init_char4 = (char *) malloc((strlen(filename) + 1024) * sizeof (char));
        strcpy(init_char3, directory);
        strcat(init_char3, filename);
        strcat(init_char3, ".lakeinit");
        // strcat(init_char4, "/input/");
        strcpy(init_char4, filename);
        strcat(init_char4, ".lakeinit");

        Lakeboundarynode_char = (char *) malloc((strlen(filename) + 1024) * sizeof (char));
        strcpy(Lakeboundarynode_char, directory);
        strcat(Lakeboundarynode_char, filename);
        strcat(Lakeboundarynode_char, ".lakeboundnode");
        Lakeboundarynode = fopen(Lakeboundarynode_char, "w");

        LakeReference_char = (char *) malloc((strlen(filename) + 1024) * sizeof (char));
        strcpy(LakeReference_char, directory);
        strcat(LakeReference_char, filename);
        strcat(LakeReference_char, ".lakereference");
        LakeReference = fopen(LakeReference_char, "w");
    }

    Elev_output_char = (char *) malloc((strlen(filename) + 1024) * sizeof (char));
    strcpy(Elev_output_char, directory);
    strcat(Elev_output_char, filename);
    strcat(Elev_output_char, ".elevation");
    Elev_output = fopen(Elev_output_char, "w");


    para_file_char = (char *) malloc((strlen(filename) + 1024) * sizeof (char));
    // strcat(para_file_char, "/input/");
    strcpy(para_file_char, filename);
    strcat(para_file_char, ".para");


    if (cData.Lake_module != 0)
    {
        LD = (LakeData) malloc(sizeof *LD);
        read_alloc_lake(filename, LD, mData);
        mData->NumLake = LD->NumLake;
    }

    printf("\n ...  PIHM 2.2 is starting ... \n");

    /* read in 9 input files with "filename" as prefix */



    /*
     * if(mData->UnsatMode ==1) {    }
     */
    if (mData->UnsatMode == 2)
    {
        /* problem size */
        N = 4 * mData->NumEle;
        mData->DummyY = (double *) malloc(N * sizeof (double));
        mData->DY = (double *) malloc(N * sizeof (double));
        /* initial state variable depending on machine */


        /* initialize mode data structure */
        initialize(filename, mData, &cData, Elev_output);
    }
    if (cData.Lake_module != 0)
    {
        /*Initialize Lake module*/
        initialize_lake(filename, LD, mData, &cData, Lakeboundarynode, LakeReference);

        /****************Lakeboundary geometry*******************************************/

        /***************************************************************************/
        tempchar = (char *) malloc((strlen(filename) + 1024) * sizeof (char));
        strcpy(tempchar, directory);
        strcat(tempchar, filename);
        strcat(tempchar, ".centroid");
        centroidfile = fopen(tempchar, "w");

        for (k = 0; k < mData->NumEle; k++)
        {
            fprintf(centroidfile, "%lf\t%lf\t%lf\t%lf\n", mData->Ele[k].x, mData->Ele[k].y, mData->Ele[k].zmax, mData->Ele[k].zmin);
        }
        fflush(centroidfile);
        fclose(centroidfile);
    }

    /*Finish initializing Lake module*/

    printf("\nSolving ODE system ... \n");


    start = clock();
    finish = clock();
    diff_t = (double) (finish - start) / CLOCKS_PER_SEC;
    /* start solver in loops */
    total_t = cData.StartTime;
    while (total_t < cData.FinishTime)
    {
        /* set start time */
        mData->t = cData.StartTime;
        for (i = 0; i < cData.NumSteps; i++)
        {
            /*
             * inner loops to next output points with ET step size
             * control
             */
            while (mData->t < cData.Tout[i + 1])
            {
                if (mData->t + cData.ETStep >= cData.Tout[i + 1])
                {
                    NextPtr = cData.Tout[i + 1];
                }
                else
                {
                    NextPtr = mData->t + cData.ETStep;
                }
                StepSize = NextPtr - mData->t;
                total_NextPtr = total_t + StepSize;
                /* calculate Interception Storage */

                printf("\n Tsteps = %f	total t = %f", mData->t, total_t);

                if (cData.Lake_module != 0)
                    Lake_land_exchange(mData, LD);

                is_sm_et(mData->t, StepSize, mData);

                f_explicit(mData->t, mData);

                if (cData.Lake_module != 0)
                    LakeDynamic(mData->t, LD, mData);
                mData->t = NextPtr;
                total_t = total_NextPtr;
                update(mData->t, mData);
            }
            Courant(mData);
            PrintData(Ofile, Ofile_Lake, &cData, mData, LD, total_t);

            if ((int) total_t % cData.surfDInt == 0)
            {
                printf("\n Tsteps = %f	total t = %f", mData->t, total_t);
                init1 = fopen(init_char1, "w");

                init2 = fopen(init_char2, "w");

                if (cData.Lake_module != 0)
                {
                    lakeinit1 = fopen(init_char3, "w");

                    lakeinit2 = fopen(init_char4, "w");
                }

                para_file = fopen(para_file_char, "w");

                for (k = 0; k < mData->NumEle; k++)
                {
                    fprintf(init1, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf", mData->EleIS[k], mData->EleSnowCanopy[k] + mData->EleSnowGrnd[k], mData->DummyY[k], mData->DummyY[k + mData->NumEle], mData->DummyY[k + 2 * mData->NumEle], mData->DummyY[k + 3 * mData->NumEle]);
                    fprintf(init1, "\n");

                    fprintf(init2, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf", mData->EleIS[k], mData->EleSnowCanopy[k] + mData->EleSnowGrnd[k], mData->DummyY[k], mData->DummyY[k + mData->NumEle], mData->DummyY[k + 2 * mData->NumEle], mData->DummyY[k + 3 * mData->NumEle]);
                    fprintf(init2, "\n");

                }

                if (cData.Lake_module != 0)
                {
                    for (k = 0; k < LD->NumLake; k++)
                    {
                        fprintf(lakeinit1, "%lf\t%lf", LD->LakeSurfDepth[k], LD->LakeGW[k]);
                        fprintf(lakeinit1, "\n");

                        fprintf(lakeinit2, "%lf\t%lf", LD->LakeSurfDepth[k], LD->LakeGW[k]);
                        fprintf(lakeinit2, "\n");
                    }
                }

                fprintf(para_file, "Verbose\t%d\n", cData.Verbose);
                fprintf(para_file, "Debug\t%d\n", cData.Debug);
                fprintf(para_file, "Hydro_init_type\t2\n");
                fprintf(para_file, "Lake_module\t%d\n", cData.Lake_module);
                fprintf(para_file, "Tide_module\t%d\n", mData->Tide_mode);
                fprintf(para_file, "Sea_level_rise\t%d\n", mData->Sea_level_rise_mode);
                fprintf(para_file, "gwDInt\t%d\n", cData.gwDInt);
                fprintf(para_file, "surfDInt\t%d\n", cData.surfDInt);
                fprintf(para_file, "snowDInt\t%d\n", cData.snowDInt);
                fprintf(para_file, "prepDInt\t%d\n", cData.prepDInt);
                fprintf(para_file, "RechInt\t%d\n", cData.RechInt);
                fprintf(para_file, "IsDInt\t%d\n", cData.IsDInt);
                fprintf(para_file, "usDInt\t%d\n", cData.usDInt);
                fprintf(para_file, "etInt\t%d\n", cData.etInt);
                fprintf(para_file, "UnsatMode\t%d\n", mData->UnsatMode);
                fprintf(para_file, "SurfMode\t%d\n", mData->SurfMode);
                fprintf(para_file, "Solver\t%d\n", cData.Solver);
                fprintf(para_file, "GSType\t%d\n", cData.GSType);
                fprintf(para_file, "MaxK\t%d\n", cData.MaxK);
                fprintf(para_file, "delt\t%f\n", cData.delt);
                fprintf(para_file, "abstol\t%f\n", cData.abstol);
                fprintf(para_file, "reltol\t%f\n", cData.reltol);
                fprintf(para_file, "InitStep\t%f\n", cData.InitStep);
                fprintf(para_file, "MaxStep\t%f\n", cData.MaxStep);
                fprintf(para_file, "ETStep\t%f\n", cData.ETStep);
                if (mData->t < cData.EndTime)
                {
                    fprintf(para_file, "StartTime\t%f\n", mData->t);
                }
                else
                {
                    fprintf(para_file, "StartTime\t%f\n", cData.StartTime);
                }
                fprintf(para_file, "EndTime\t%f\n", cData.EndTime);
                fprintf(para_file, "FinishTime\t%f\n", cData.FinishTime);
                fprintf(para_file, "outtype\t%d\n", cData.outtype);
                fprintf(para_file, "a\t%f\n", cData.a);
                fprintf(para_file, "b\t%f\n", cData.b);



                fflush(init1);
                fclose(init1);
                fflush(init2);
                fclose(init2);
                fflush(para_file);
                fclose(para_file);
                if (cData.Lake_module != 0)
                {
                    fflush(lakeinit1);
                    fclose(lakeinit1);
                    fflush(lakeinit2);
                    fclose(lakeinit2);
                }

            }
        }
        for (k = 0; k < mData->NumPrep; k++)
        {
            mData->TSD_Prep[k].iCounter = 0;
        }

        for (k = 0; k < mData->NumTemp; k++)
        {
            mData->TSD_Temp[k].iCounter = 0;
        }

        for (k = 0; k < mData->NumHumidity; k++)
        {
            mData->TSD_Humidity[k].iCounter = 0;
        }

        for (k = 0; k < mData->NumWindVel; k++)
        {
            mData->TSD_WindVel[k].iCounter = 0;
        }

        for (k = 0; k < mData->NumRn; k++)
        {
            mData->TSD_Rn[k].iCounter = 0;
        }

        for (k = 0; k < mData->NumG; k++)
        {
            mData->TSD_G[k].iCounter = 0;
        }

        for (k = 0; k < mData->NumP; k++)
        {
            mData->TSD_Pressure[k].iCounter = 0;
        }

        for (k = 0; k < mData->NumLC; k++)
        {
            mData->TSD_LAI[k].iCounter = 0;
        }

        for (k = 0; k < mData->NumLC; k++)
        {
            mData->TSD_RL[k].iCounter = 0;
        }

        for (k = 0; k < mData->NumMeltF; k++)
        {
            mData->TSD_MeltF[k].iCounter = 0;
        }

        for (k = 0; k < mData->NumTide; k++)
        {
            mData->TSD_Tide[k].iCounter = 0;
        }
    }


    /* Free memory */

    //N_VDestroy_Serial(CV_Ydot);
    /* Free integrator memory */

    FreeData(mData, &cData, LD);

    for (i = 0; i < 28; i++)fclose(Ofile[i]);

    for (i = 0; i < 28; i++)free(ofn[i]);

    if (cData.Lake_module != 0)
    {
        for (i = 0; i < 10; i++)fclose(Ofile_Lake[i]);

        for (i = 0; i < 10; i++)free(ofn_Lake[i]);
    }
    free(filename);
    free(mData);
    free(LD);
    free(init_char1);
    free(init_char2);
    free(init_char3);
    free(init_char4);
    free(directory);
    return 0;
}

void Courant(Model_Data MD)
{
    int i, j;
    double u;

    for (i = 0; i < MD->NumEle; i++)
    {
        for (j = 0; j < 3; j++)
        {

            if (MD->CrossA_surf[i][j] > 0 && MD->Distance[i][j] > 0)
            {
                u = MD->FluxSurf[i][j] / MD->CrossA_surf[i][j];
                MD->Courant_surf[i][j] = fabs(u / (MD->Distance[i][j]*1440));
            }
            else
            {
                MD->Courant_surf[i][j] = 0;
            }

            if (MD->CrossA_sub[i][j] > 0 && MD->Distance[i][j] > 0)
            {
                u = MD->FluxSub[i][j] / MD->CrossA_sub[i][j];
                MD->Courant_sub[i][j] = fabs(u / (MD->Distance[i][j]*1440));
            }
            else
            {
                MD->Courant_sub[i][j] = 0;
            }

        }
    }
}
