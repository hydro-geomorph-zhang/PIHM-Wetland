/**********************************************************************************
 * File        : pihm.h                                                           *
 * Function    : Declaration and Definition of global variables and data structure*
 * Developer of PIHM 2.0: Mukesh Kumar (muk139@psu.edu)                           *
 * Developer of PIHM 1.0: Yizhong Qu   (quyizhong@gmail.com)                      *
 * Version     : Nov, 2007 (2.0)                                                  *
 *--------------------------------------------------------------------------------*
 *                                                                                *
 *..............MODIFICATIONS/ADDITIONS in PIHM 2.0...............................*
 * a) Definition of new variables for ELEMENT data model related to 		  *
 *	i. Subsurface: KsatH, KsatV, infKsatV,infD, RzD, macD, macKsatH, macKsatV *
 *	vAreaF, hAreaF								  *
 *	ii. ET: LAImax, vegFrac, Albedo, Rs_ref, Rmin				  *
 *	iii. Snow: meltF							  *
 *	iv. Surface: dhBydx, dhBYdy, surfH					  *
 * b) Definition of new variables for RIVER data model				  *
 *	i. KsatH, KsatV, bedThick 						  *
 * c) Definition of new structures:						  *
 *	i. geology								  *
 *	ii. Land Cover								  *
 *	iii. Calibration							  *
 * d) Definition of New Control Parameters for each file			  *
 *--------------------------------------------------------------------------------*
 * For questions or comments, please contact                                      *
 *      --> Mukesh Kumar (muk139@psu.edu)                                	  *
 *      --> Prof. Chris Duffy (cxd11@psu.edu)                                     *
 * This code is free for research purpose only.                                   *
 * Please provide relevant references if you use this code in your research work  *
 *--------------------------------------------------------------------------------*
 **********************************************************************************/

#include "sundials_types.h"
#include "nvector_serial.h"

/* Definition of Global Variable Types */
#define FULLY_COUPLE
// #define LE_PIHM
// #define BEDROCK
//#define LE_PIHM_SED
#define LE_PIHM_HYDRO
#define LE_PIHM_SUBHYDRO
#define REALTIME
// #define MEAN

typedef struct element_type { /* Data model for a triangular element */
    int index; /* Element No. */
    int node[3]; /* anti-clock-wise */
    int nabr[3]; /* neighbor i shares edge i (0: on boundary) */

    double edge[3]; /* edge i is from node i to node i+1 */
    double area; /* area of element */

    double x; /* x of centroid */
    double y; /* y of centroid */
    double zmin; /* z_min of centroid */
    double zmax; /* z_max of centroid */


    double KsatH; /* horizontal geologic saturated hydraulic
				 * conductivity */
    double KsatV; /* vertical geologic saturated hydraulic
				 * conductivity */
    double infKsatV; /* vertical surface saturated
					 * hydraulic conductivity */
    double Porosity;
    double infD; /* depth from ground surface accross which
				 * head is calculated during infiltration */
    double Alpha; /* Alpha from van-genuchten eqn which is
				 * given by satn =
				 * 1/pow(1+pow(abs(Alpha*psi),Beta),1-1/Beta) */
    double Beta;
    double RzD; /* Root zone depth */
    double macD; /* macropore Depth */
    double macKsatH; /* macropore horizontal saturated
					 * hydraulic conductivity */
    double macKsatV; /* macropore vertical saturated
					 * hydraulic conductivity */
    double vAreaF; /* macropore area fraction on a vertical
				 * cross-section */
    double hAreaF; /* macropore area fraction on a horizontal
				 * cross-section */
    int Macropore; /* 1: macropore; 0: regular soil */

    double LAImax; /* maxm. LAI accross all seasons for a
				 * vegetation type */
    double VegFrac; /* areal vegetation fraction in a triangular
				 * element */
    double Albedo; /* albedo of a triangular element */
    double Rs_ref; /* reference incoming solar flux for
				 * photosynthetically active canopy */
    double Rmin; /* minimum canopy resistance */
    double Rough; /* surface roughness of an element */

    double windH; /* wind measurement height */


    int soil; /* soil type */
    int geol; /* geology type */
    int LC; /* Land Cover type  */
    int IC; /* initial condition type */
    int BC[3]; /* boundary type. 0:natural bc (no flow);
				 * 1:Dirichlet BC; 2:Neumann BC */
    int prep; /* precipitation (forcing) type */
    int temp; /* temperature (forcing) type   */
    int humidity; /* humidity type */
    int WindVel; /* wind velocity type  */
    int Rn; /* net radiation input */
    int G; /* radiation into ground */
    int pressure; /* pressure type */
    int source; /* source (well) type */
    int meltF; /* meltFactor */
    /* for calculation of dh/ds */
    double surfH[3]; /* Total head in neighboring cells */
    double surfX[3]; /* Center X location of neighboring
					 * cells */
    double surfY[3]; /* Center Y location of neighboring
					 * cells */
    double dhBYdx; /* Head gradient in x dirn. */
    double dhBYdy; /* Head gradient in y dirn. */
    int *Nabr_index; /* the index of the Nabour element, 0,1 or 2*/

    //#ifdef LE_PIHM
    double DReg; /*diffusivity of regolith*/
    double P0; /*Regolith production rate in the absence of soil above bedrock*/
    double Uplift; /*Bedrock uplift rate*/
    double CoefP0; /*Fitting constant in weather equation*/
    double RhoReg;
    double RhoBed;
    double Diameter; //The diameter of sediment particle
    //#endif


} element;

typedef struct nodes_type { /* Data model for a node */
    int index; /* Node no. */

    double x; /* x coordinate */
    double y; /* y coordinate */
    double zmin; /* z bed rock elevation */
    double zmax; /* z surface elevation  */

} nodes;

typedef struct element_IC_type {/* Initial state variable conditions on each
				 * element */
    int index;

    double interception; /* Interception storage (Note all
					 * these variables have dimension of
					 * L */
    double snow; /* Snow depth */
    double surf; /* Overland flow depth */
    double unsat; /* unsaturated zone depth */
    double sat; /* saturated zone depth */

    //#ifdef LE_PIHM
    double groundelev; /* elevation of ground surface*/
    double bedrockelev; /* elevation of bedrock */
    //#endif

#ifdef LE_PIHM_SED
    double groundelev; /* elevation of ground surface*/
#endif

} element_IC;

typedef struct soils_type {
    int index; /* index */

    double KsatV; /* vertical saturated soil conductivity */
    double ThetaS; /* soil porosity */
    double ThetaR; /* soil moisture residual */
    double Alpha; /* soil curve parameter 1 */
    double Beta; /* soil curve parameter 2 */

    double hAreaF; /* macroporous area fraction on horizontal
				 * section */
    double macKsatV; /* macroporous saturated vertical
					 * conductivity */

    double infD; /* depth from ground surface accross which
				 * head is calculated during infiltration */
} soils;

typedef struct LE_soils {
    int index; /* index */

    //#ifdef LE_PIHM
    double RhoReg; /*Density of regolith*/
    double DReg; /*Diffusivity*/
    double Diameter; //The diameter of sediment particle
    //#endif

#ifdef LE_PIHM_SED
    double RhoReg; /*Density of regolith*/
    double Diameter; //The diameter of sediment particle
#endif
} LE_soils;

typedef struct geol_type {
    int index; /* index */

    double KsatH; /* horizontal saturated geology conductivity */
    double KsatV; /* vertical saturated geology conductivity */
    double ThetaS; /* geology porosity */
    double ThetaR; /* residual porosity */
    double Alpha; /* van genuchten parameter */
    double Beta; /* van genuchten parameter */

    double vAreaF; /* macroporous area fraction on vertical
				 * section */
    double macKsatH; /* macroporous saturated horizontal
					 * conductivity */
    double macD;

#ifdef LE_PIHM_SED

#endif



} geol;

typedef struct LE_bedrock
 {
    int index; /* index */
    //#ifdef LE_PIHM
    double RhoBed; /*Density of bedrock*/
    double P0; /*Regolith production rate in the absence of soil above bedrock*/
    double Uplift; /*Bedrock uplift rate*/
    double CoefP0; /*Fitting constant in weather equation*/
    //#endif	

} LE_bedrock;

typedef struct lc_type {
    int index; /* index */

    double LAImax; /* max LAI */
    double VegFrac; /* Canopy Fracn */
    double Albedo; /* Albedo */
    double Rs_ref;
    double Rmin; /* Minimum stomatal resistance */
    double Rough; /* surface roughness factor  */
    double RzD; /* rootZone Depth */
} LC;

typedef struct river_segment_type {
    int index;
    int FromNode; /* Upstream Node no. */
    int ToNode; /* Dnstream Node no. */
    int LeftEle; /* Left neighboring element */
    int RightEle; /* Right neighboring element */
    int down; /* down stream segment */
} river_segment;

typedef struct TSD_type {
    char name[5];
    int index;
    int length; /* length of time series */
    int iCounter; /* interpolation counter */
    double **TS; /* 2D time series data */

} TSD;

typedef struct global_calib {
    double KsatH; /* For explanation of each calibration
				 * variable, look for corresponding variables
				 * above */
    double KsatV;
    double infKsatV;
    double macKsatH;
    double macKsatV;
    double infD;
    double RzD;
    double macD;
    double Porosity;
    double Alpha;
    double Beta;
    double vAreaF;
    double hAreaF;
    double Temp;
    double Prep;
    double VegFrac;
    double Albedo;
    double Rough;
    double LakeKh;
    double LakeKv;
    double LakeThetaS;
    double LakeThetaR;
    double LakemacKh;

    //#ifdef LE_PIHM
    double DReg; /*diffusivity of regolith*/
    double P0; /*Regolith production rate in the absence of soil above bedrock*/
    double CoefP0; /*Fitting constant in weather equation*/
    //#endif

#ifdef LE_PIHM_SED

#endif

} globalCal;

typedef struct process_control {
    double Et0;
    double Et1;
    double Et2;
} processCal;

typedef struct model_data_structure { /* Model_data definition */
    int UnsatMode; /* Unsat Mode */
    int SurfMode; /* Surface Overland Flow Mode */

    int NumEle; /* Number of Elements */
    int NumNode; /* Number of Nodes    */
    int NumRiv; /* Number of Rivere Segments */

    int NumPrep; /* Number of Precipatation time series types  */
    int NumTemp; /* Number of Temperature time series types      */
    int NumHumidity; /* Number of Humidity time series
					 * types         */
    int NumWindVel; /* Number of Wind Velocity time
					 * series types    */
    int NumRn; /* Number of Net Radiation time series types    */
    int NumG; /* Number of Ground Heat time series types      */
    int NumP; /* Number of Pressure time series types         */
    int NumSource; /* Number of Source time series types           */

    int NumSoil; /* Number of Soils           */
    int NumLESoil; /* Number of Soils           */
    int NumLEBedrock; /* Number of Soils           */
    int NumGeol; /* Number of Geologies           */
    int NumRes; /* Number of Reservoir       */
    int NumLC; /* Number of Land Cover Index Data */

    int NumMeltF; /* Number of Melt Factor Time series */

    int Num1BC; /* Number of Dirichlet BC    */
    int Num2BC; /* Number of Numann BC       */
    int NumEleIC; /* Number of Element Initial Condtion */

    int NumOutlet; /*Number of element at outlet */
    double MSF; /* Morphological Scaling Factor */
    int computtime;
    int NumLake;
    int NumTide;

    element *Ele; /* Store Element Information  */
    nodes *Node; /* Store Node Information     */
    element_IC *Ele_IC; /* Store Element Initial Condtion */
    soils *Soil; /* Store Soil Information     */
    LE_soils *LE_soil;
    geol *Geol; /* Store Geology Information     */
    LE_bedrock *LE_bedrock;
    LC *LandC; /* Store Land Cover Information */

    river_segment *Riv; /* Store River Segment Information */

    TSD *TSD_Inc; /* Infiltration Capacity Time Series Data */
    TSD *TSD_LAI; /* Leaves Area Index Time Series Data     */
    //TSD * TSD_DH;	/* Zero plane Displacement Height */
    TSD *TSD_RL; /* Roughness Length */
    double *ISFactor; /* ISFactor is used to calculate
					 * ISMax from LAI */
    double *windH; /* Height at which wind velocity is measured */
    TSD *TSD_MeltF; /* Monthly Varying Melt Factor for
					 * Temperature Index model */

    TSD *TSD_EleBC; /* Element Boundary Condition Time
					 * Series Data  */
    TSD *TSD_Prep; /* RainFall Time Series Data       */
    TSD *TSD_Temp; /* Temperature Time Series Data    */
    TSD *TSD_Humidity; /* Humidity Time Series Data       */
    TSD *TSD_WindVel; /* Wind Velocity Time Series Data  */
    TSD *TSD_Rn; /* Net Radiation Time Series Data  */
    TSD *TSD_G; /* Radiation into Ground Time Series Data */
    TSD *TSD_Pressure; /* Vapor Pressure Time Series data       */
    TSD *TSD_Source; /* Source (well) Time Series data  */
    TSD *TSD_Tide; /* tide information*/

    double **FluxSurf; /* Overland Flux   */
    double **FluxSub; /* Subsurface flow Flux   */
    double **FluxSalt; /* Subsurface Saltwater flux*/
    double **Seepage; /* GW seepage Flux   */
    double **Courant_surf; /* Courant Flux   */
    double **Courant_sub; /* Courant Flux   */
    double **CrossA_surf;
    double **CrossA_sub;
    double **Distance;

    double *ElePrep; /* Precep. on each element */
    double *EleETloss;
    double *EleNetPrep; /* Net precep. on each elment */
    double *EleViR; /* Variable infiltration rate */
    double *EleViR_Salt; /* Variable infiltration rate */
    double *Recharge; /* Recharge rate to GW */
    double *EleSnow; /* Snow depth on each element */
    double *EleSnowGrnd; /* Snow depth on ground element */
    double *EleSnowCanopy; /* Snow depth on canopy element */
    double *EleIS; /* Interception storage */
    double *EleISmax; /* Maximum interception storage
					 * (liquid precep) */
    double *EleISsnowmax; /* Maximum interception storage
					 * (snow) */
    double *EleTF; /* Through Fall */
    double **EleET; /* Evapo-transpiration (from canopy, ground,
				 * subsurface, transpiration) */
    double *EleTemp;
    double *ElePET; 
    
    double *DY; /*store the DY for surface water, surface elevation and bed elevation*/
    double *DummyY;
    double *DummyY_Hydro;
    double *DummyY_LE;
    double *PrintVar[31];
    double *Total_Surfwater_out;
    double *Total_Unsatwater_out;
    double *Total_Satwater_out;
    double *Total_Saltwater_out;
    double *infil_mode;
    int Tide_mode;
    int Sea_level_rise_mode;

    int **Outlet_location;
    int *Land_lake_index; /*stores the lake ID surrounded by an element*/
    processCal pcCal;
    double *SurfDepth; /*Surface water*/
    double *GW; /*Ground water*/
    double *UGW; /*Unsaturated water*/
    double *grdelev; /*Unsaturated water*/
    double *bedelev; /*Unsaturated water*/
    double *LakeSurfElev; /*Unsaturated water*/
    double *LakeGWElev; /*Unsaturated water*/
    double *LakeInitialElev;
    double *LakeOutletDecline;
    double *BankElevAdjust;
    double *LakeBed;
    double *Weir_length;
    double *Orifice_height;
    double t; /*current simulation time*/
} *Model_Data;

typedef struct control_data_structure {
    int Verbose;
    int Debug;

    int Solver; /* Solver type */
    int NumSteps; /* Number of external time steps
					 * (when results can be printed) for
					 * the whole simulation */

    int gwD; /* File boolean, Choose 1 if your want to
				 * print ground water */
    int prepD; /* File boolean, Choose 1 if your want to print precipitation */
    int surfD; /* File boolean, Choose 1 if your want to
				 * print overland flow */
    int snowD; /* File boolean, Choose 1 if your want to
				 * print snow Depth */
    int Rech; /* File boolean, Choose 1 if your want to
				 * print recharge to ground water */
    int IsD; /* File boolean, Choose 1 if your want to
				 * print interception depth */
    int usD; /* File boolean, Choose 1 if your want to
				 * print unsaturated depth */
    int et[3]; /* File boolean, Choose 1 if your want to
				 * print individual evapo-transpiration
				 * components */

    int rivFlx; /* File boolean, Choose 1 if your
					 * want to print river/river bed
					 * fluxes */

    int gwDInt; /* Time interval to output average val of
				 * variables */
    int prepDInt;
    int surfDInt;
    int snowDInt;
    int RechInt;
    int IsDInt;
    int usDInt;
    int etInt;
    int computtimeDInt;

    int init_type; /* initialization mode */



    double abstol; /* absolute tolerance */
    double reltol; /* relative tolerance */
    double InitStep; /* initial step size */
    double MaxStep; /* Maximum step size */
    double ETStep; /* Step for et from interception */

    int GSType, MaxK; /* Maximum Krylov order */
    double delt;

    double StartTime; /* Start time of simulation */
    double EndTime; /* End time of simulation */
    double FinishTime; /*Time of the total run, if FinishTim is greater than EndTime, time t will start from StartTime again. */


    int outtype;
    double a; /* External time stepping controls */
    double b;

    double *Tout;
    int Lake_module;


    globalCal Cal; /* Convert this to pointer for localized
				 * calibration */
    //#ifdef LE_PIHM
    
    int read_centroid;
    
    int groundelevDInt;
    int rockelevDInt;
    int bedloadDInt;
    int soilsubDInt;
    int upliftDInt;
    int weatheringDInt;
    int continue_LEM; /* Use previous results as initial condition running Landscape evolution model again. 0 don't use, 1 use.*/
    //#endif

#ifdef LE_PIHM_SED

#endif
} Control_Data;





