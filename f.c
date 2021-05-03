/*******************************************************************************
 *-----------------------------------------------------------------------------*
 * File        : f.c   (PIHM v.2.0)                                            *
 * Function    : Model Kernel: Building ODE system for each physical process   *
 *-----------------------------------------------------------------------------*
 *                                                                             *
 *-----------------------------------------------------------------------------*
 * Developer of PIHM v.2.0:  Mukesh Kumar (muk139@psu.edu)		       *
 * Developer of PIHM v.1.0:  Yizhong Qu	(quyizhong@gmail.com)		       *
 *-----------------------------------------------------------------------------*
 *                                                                             *
 *-----------------------------------------------------------------------------*
 * NOTE: f.c has gone a massive revamp (essentially rewritten) since PIHM v.1.0*
 *                                                                             *
 *...........MODIFICATIONS/ADDITIONS incorporated in f.c (PIHM v.2.0)...........*
 * a) Surface Flow: 							       *
 *	--> Correction of diffusion wave approximation (calculation of dh/ds)  *
 *              i. Calculation of dh/ds performed using planar slope connecting*
 *                 neighboring centroids				       *
 *              ii.Reflection of elements at boundaries and rivers for dh/ds   *
 *		   calculation
 *	--> Correction of kinematic wave approximation (dh/ds calculation based*
 *	    on elevation only instead of head				       *
 *	--> Correction of gradient for cases with steep change in topography   *
 * b) Subsurface Flow:							       *
 *	--> Addition of macropore phenomena				       *
 *	--> Addition of rectangular cell beneath a river element	       *
 *	--> Implementation of two layered subsurface model(sat/unsat) based on *
 *	Richard's eqn							       *
 *	--> Incorporation of Vertical and Horizontal Anisotropy                *
 *	--> Use of geologic data					       *
 * c) River Flow:							       *
 *	--> Correction of kinematic and diff. wave approximation of SV eqn     *
 *	--> Incorporation of flexible river shapes			       *
 *	--> Separate incorporation of leakage and lateral flow		       *
 *	--> Correction of bank overland flow for extreme cases		       *
 *	--> Addition of aquifer cells below river elements		       *
 * c) Surface/Subsurface Coupling:					       *
 *	--> Implementation of First order coupling through (in/ex)filtration   *
 *		based on head continuity 				       *
 * d) Evaporation:							       *
 *	--> Incorporation of ET from ground/subsurface/vegetation	       *
 *	--> Incorporation of landcover properties for calculation of each ET   *
 *	    component							       *
 * e) Computational:							       *
 *	--> Use of temporary state variables in calculation. Note: Never change*
 *		core state variables					       *
 * f) Miscellaneous (other advantages realtive to PIHM1.0): No maximum         *
 *    constraint on gw level. Accordingly, no numerical constraints on subsur- *
 *    face flux terms.Faster Implementation. Led to first large scale model    *
 *    application.
 *-----------------------------------------------------------------------------*
 *									       *
 *-----------------------------------------------------------------------------*
 * For questions or comments, please contact 				       *
 *	--> Mukesh Kumar (muk139@psu.edu)				       *
 *	--> Prof. Chris Duffy (cxd11@psu.edu)				       *
 * This code is free for research purpose only.		       		       *
 * Please provide relevant references if you use this code in your research work*
 *-----------------------------------------------------------------------------*
 *                                                                             *
 * REFERENCES:                                                                 *
 * PIHM2.0:                                                                    *
 *      a) Kumar, M., 2008, "Development and Implementation of a Multiscale,   *
 *      Multiprocess Hydrologic Model". PhD Thesis, Penn State University      *
 *      b) Kumar, M, G.Bhatt & C.Duffy, "Coupling of Data and Processes in     *
 *      a Mesoscale Watershed", Advances in Water Resources (submitted)        *
 * PIHM1.0:                                                                    *
 *      a) Qu, Y., 2005, "An Integrated hydrologic model for multiproces       *
 *      simulation using semi-discrete finite volume approach".PhD Thesis, PSU *
 *      b) Qu, Y. & C. Duffy, 2007, "A semidiscrete finite volume formulation  *
 *      for multiprocess watershed simulation". Water Resources Research       *
 *******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


#include "nvector_serial.h"
#include "sundials_types.h"

#include "pihm.h"
#include "lake.h"
#define multF	2
#define MINpsi	-70
#define EPS 0.05
#define THRESH 0.0
#define UNIT_C 1440		/* Note 60*24 for calculation of yDot in
				 * m/min units while forcing is in m/day. */
#define GRAV 9.8*60*60		/* Note the dependence on physical units */

#define C_air 1004.0
#define Lv (2.503*pow(10,6))
#define SIGMA (5.67*pow(10,-8)*60)
#define R_dry 287.04
#define R_v 461.5
#define TIME_CONVERT 525600

#define Niu 0.000001004 //Unit: (m^2/s)

double        Interpolation(TSD * Data, double t);


/* Note: above formulation doesnt' satisfies upwind/downwind scheme. */
double
avgY(double diff, double yi, double yinabr)
{
	if (diff > 0) {
		if (yi > 1 * EPS / 100) {
			//return 0.5 * (yi + yinabr);
			//return ((yinabr > yi) ? 0 : 1.0 * yi);
			/*
			 * Note the if-else TRUE case can be possible only
			 * for Kinematic case
			 */
			return 1.0 * yi;
		} else {
			return 0;
		}
	} else {
		if (yinabr > 1 * EPS / 100) {
			//return 0.5 * (yi + yinabr);
			//return ((yi > yinabr) ? 0 : 1.0 * yinabr);
			/*
			 * Note the if-else TRUE case can be possible only
			 * for Kinematic case
			 */
			return 1.0 * yinabr;
		} else {
			return 0;
		}
	}
}
double
effKV(double ksatFunc, double gradY, double macKV, double KV, double areaF)
{
	if (ksatFunc >= 0.98) {
		return (macKV * areaF + KV * (1 - areaF) * ksatFunc);
	} else {
		if (fabs(gradY) * ksatFunc * KV <= 1 * KV * ksatFunc) {
			return KV * ksatFunc;
		} else {
			if (fabs(gradY) * ksatFunc * KV < (macKV * areaF + KV * (1 - areaF) * ksatFunc)) {
				return (macKV * areaF * ksatFunc + KV * (1 - areaF) * ksatFunc);
			} else {
				return (macKV * areaF + KV * (1 - areaF) * ksatFunc);
			}
		}
	}
}
double
effKH(int mp, double tmpY, double aqDepth, double MacD, double MacKsatH, double areaF, double ksatH)
{
	if (mp == 1) {
		if (tmpY > aqDepth - MacD) {
			if (tmpY > aqDepth) {
				return (MacKsatH * MacD * areaF + ksatH * (aqDepth - MacD * areaF)) / aqDepth;
			} else {
				return (MacKsatH * (tmpY - (aqDepth - MacD)) * areaF + ksatH * (aqDepth - MacD + (tmpY - (aqDepth - MacD)) * (1 - areaF))) / tmpY;
			}
		} else {
			return ksatH;
		}
	} else {
		return ksatH;
	}
}


double* Calculate_DY(Model_Data MD,double *DY)
{
	int             i,j, k, inabr,boundary_index,q,is_outlet,out_i,out_j;
	double        Delta, Gamma;
	double        Rn, G, T, Vel, RH, VP, P, LAI, zero_dh, cnpy_h, rl,
	                r_a, r_s, alpha_r, f_r, eta_s, beta_s, gamma_s, Rmax,
	                Lambda, P_c, qv, qv_sat, ETp;
	double        ThetaRef, ThetaW;
	double        Avg_Y_Surf, Dif_Y_Surf, Grad_Y_Surf, Avg_Sf, Distance,Dif_Y_Surf1,Dif_Y_Surf2;
	double        Cwr, TotalY_Riv, TotalY_Riv_down, CrossA, CrossAdown, AvgCrossA,
	                Perem, Perem_down, Avg_Rough, Avg_Perem, Avg_Y_Riv,
	                Dif_Y_Riv, Grad_Y_Riv, Wid, Wid_down, Avg_Wid;
	double        Avg_Y_Sub, Dif_Y_Sub, Avg_Ksat, Grad_Y_Sub, nabrAqDepth, AquiferDepth,
	                Deficit, elemSatn, satKfunc, effK, effKnabr, TotalY_Ele,
	                TotalY_Ele_down;
	double		Dif_Elev, Grad_Elev;
	double 		c1_square,c2_square,L,b;
	double		time_convert = 60*60*24/UNIT_C;
	double		h_regolith, Tau0, Critical_ShieldsStress,U_star, Renolds, ShieldsStress, q_star,R,q_b,Total_soil_out, Total_soil_in;
	double		Cf = 0.003;
	double		Cos_slope, Sec_slope;
	
	for (i=0;i<MD->NumEle;i++)
	{
        for (j = 0; j < 3; j++)
		{			
			DY[i] = DY[i] - MD->FluxSurf[i][j] / MD->Ele[i].area - MD->Seepage[i][j] / MD->Ele[i].area;
			DY[i + 2 * MD->NumEle] = DY[i + 2 * MD->NumEle] - MD->FluxSub[i][j] / MD->Ele[i].area;
            DY[i + 3 * MD->NumEle] = DY[i + 3 * MD->NumEle] - MD->FluxSalt[i][j] / MD->Ele[i].area;
		}
		DY[i] = DY[i] + MD->EleNetPrep[i] - MD->EleViR[i];
		DY[i] = DY[i] - MD->EleET[i][1];
		DY[i] = DY[i] / (UNIT_C);

		DY[i + MD->NumEle] = DY[i + MD->NumEle] + MD->EleViR[i] - MD->Recharge[i];
		DY[i + 2 * MD->NumEle] = DY[i + 2 * MD->NumEle] + MD->Recharge[i];
		DY[i + MD->NumEle] = DY[i + MD->NumEle] - MD->EleET[i][2];
		DY[i + 2 * MD->NumEle] = DY[i + 2 * MD->NumEle] - MD->EleET[i][3];
		DY[i + MD->NumEle] = DY[i + MD->NumEle] - MD->EleET[i][4];
		DY[i + 2 * MD->NumEle] = DY[i + 2 * MD->NumEle] - MD->EleET[i][5];
			
		DY[i + MD->NumEle] =  DY[i + MD->NumEle] / (MD->Ele[i].Porosity * UNIT_C);
		DY[i + 2 * MD->NumEle] =  DY[i + 2 * MD->NumEle] / (MD->Ele[i].Porosity * UNIT_C);
	}
	return DY;
}

void Adjust_water(Model_Data MD)
{
	int	i, j, inabr;
	
	for (i = 0; i < MD->NumEle; i++) 
	{

        //-----------Surface water------------------//
		if (MD->Total_Surfwater_out[i]>0&&MD->Total_Surfwater_out[i]/UNIT_C>MD->DummyY[i]&&MD->DummyY[i]>=0)
		{
			for (j=0;j<3;j++)
			{
				inabr = MD->Ele[i].nabr[j] - 1; 
				
				if (MD->FluxSurf[i][j]>0)
				{
					MD->FluxSurf[i][j] = ((MD->FluxSurf[i][j]/MD->Ele[i].area)/MD->Total_Surfwater_out[i])*MD->DummyY[i]*UNIT_C*MD->Ele[i].area;
					if (inabr>=0)
						MD->FluxSurf[inabr][MD->Ele[inabr].Nabr_index[i]] = -MD->FluxSurf[i][j];
				}
				if (MD->Seepage[i][j]>0)
				{
					MD->Seepage[i][j] = ((MD->Seepage[i][j]/MD->Ele[i].area)/MD->Total_Surfwater_out[i])*MD->DummyY[i]*UNIT_C*MD->Ele[i].area;
					if (inabr>=0)
						MD->FluxSub[inabr][MD->Ele[inabr].Nabr_index[i]] = -MD->Seepage[i][j];
				}	
			}
			if (MD->EleET[i][1]>0)
			{
				MD->EleET[i][1] = ((MD->EleET[i][1])/MD->Total_Surfwater_out[i])*MD->DummyY[i]*UNIT_C;
			}
			if (MD->EleViR[i]>0)
			{
				MD->EleViR[i] = ((MD->EleViR[i])/MD->Total_Surfwater_out[i])*MD->DummyY[i]*UNIT_C;
				if (MD->infil_mode[i]==-1)
				{
					MD->Recharge[i] = MD->EleViR[i];
				}
			}
			
		}
		else if(MD->DummyY[i]<0)
		{
			for (j=0;j<3;j++)
			{
				inabr = MD->Ele[i].nabr[j] - 1; 
				if (MD->FluxSurf[i][j]>0)
				{
					MD->FluxSurf[i][j] = 0;
					if (inabr>=0)
						MD->FluxSurf[inabr][MD->Ele[inabr].Nabr_index[i]] = -MD->FluxSurf[i][j];
				}
				if (MD->Seepage[i][j]>0)
				{
					MD->Seepage[i][j] = 0;
					if (inabr>=0)
						MD->FluxSub[inabr][MD->Ele[inabr].Nabr_index[i]] = -MD->Seepage[i][j];
				}
			}
			if (MD->EleET[i][1]>0)
			{
				MD->EleET[i][1] = 0;
			}
			if (MD->EleViR[i]>0)
			{
				MD->EleViR[i] = 0;
				if (MD->infil_mode[i]==-1)
				{
					MD->Recharge[i] = MD->EleViR[i];
				}
			}
		}
		//-----------Finish surface water------------------//
		
			//-----------Unsat water------------------//
			if (MD->Total_Unsatwater_out[i]>0 && MD->Total_Unsatwater_out[i]/UNIT_C>MD->DummyY[i + MD->NumEle]*MD->Ele[i].Porosity && MD->DummyY[i + MD->NumEle]>=0)
			{
				if (MD->EleET[i][2]>0)
				{
					MD->EleET[i][2] = ((MD->EleET[i][2])/MD->Total_Unsatwater_out[i])*MD->DummyY[i + MD->NumEle]*MD->Ele[i].Porosity*UNIT_C;
				}
				if (MD->EleET[i][4]>0)
				{
					MD->EleET[i][4] = ((MD->EleET[i][4])/MD->Total_Unsatwater_out[i])*MD->DummyY[i + MD->NumEle]*MD->Ele[i].Porosity*UNIT_C;
				}
				if (MD->Recharge[i]>0&&MD->infil_mode[i]==1)
				{
					MD->Recharge[i] = ((MD->Recharge[i])/MD->Total_Unsatwater_out[i])*MD->DummyY[i + MD->NumEle]*MD->Ele[i].Porosity*UNIT_C;
				}
						
			}
			else if (MD->DummyY[i + MD->NumEle]<0)
			{
				if (MD->EleET[i][2]>0)
				{
					MD->EleET[i][2] = 0;
				}
				if (MD->EleET[i][4]>0)
				{
					MD->EleET[i][4] = 0;
				}
				if (MD->Recharge[i]>0&&MD->infil_mode[i]==1)
				{
					MD->Recharge[i] = 0;
				}
			}
			//-----------Finish unsat water------------------//
			
			//-----------Saturated water------------------//
			if (MD->Total_Satwater_out[i]>0 && MD->Total_Satwater_out[i]/UNIT_C>MD->DummyY[i + 2*MD->NumEle]*MD->Ele[i].Porosity && MD->DummyY[i + 2*MD->NumEle]>=0)
			{
				for (j=0;j<3;j++)
				{
					inabr = MD->Ele[i].nabr[j] - 1; 
					
					if (MD->FluxSub[i][j]>0)
					{
						MD->FluxSub[i][j] = ((MD->FluxSub[i][j]/MD->Ele[i].area)/MD->Total_Satwater_out[i])*MD->DummyY[i + 2*MD->NumEle]*MD->Ele[i].Porosity*UNIT_C*MD->Ele[i].area;
						if (inabr>=0)
							MD->FluxSub[inabr][MD->Ele[inabr].Nabr_index[i]] = -MD->FluxSub[i][j];
					}	
				}				
				
				if (MD->EleET[i][3]>0)
				{
					MD->EleET[i][3] = ((MD->EleET[i][3])/MD->Total_Satwater_out[i])*MD->DummyY[i + 2*MD->NumEle]*MD->Ele[i].Porosity*UNIT_C;
				}
				if (MD->EleET[i][5]>0)
				{
					MD->EleET[i][5] = ((MD->EleET[i][5])/MD->Total_Satwater_out[i])*MD->DummyY[i + 2*MD->NumEle]*MD->Ele[i].Porosity*UNIT_C;
				}
				
				if (MD->EleViR[i]<0)
				{
					MD->EleViR[i] = ((MD->EleViR[i])/MD->Total_Satwater_out[i])*MD->DummyY[i + 2*MD->NumEle]*MD->Ele[i].Porosity*UNIT_C;
					MD->Recharge[i] = MD->EleViR[i];
				}
				if (MD->Recharge[i]<0&&MD->infil_mode[i]==1)
				{
					MD->Recharge[i] = ((MD->Recharge[i])/MD->Total_Satwater_out[i])*MD->DummyY[i + 2*MD->NumEle]*MD->Ele[i].Porosity*UNIT_C;
				}						
			}
			else if (MD->DummyY[i + 2*MD->NumEle]<0)
			{
				for (j=0;j<3;j++)
				{
					inabr = MD->Ele[i].nabr[j] - 1; 
					
					if (MD->FluxSub[i][j]>0)
					{
						MD->FluxSub[i][j] = 0;
						if (inabr>=0)
							MD->FluxSub[inabr][MD->Ele[inabr].Nabr_index[i]] = -MD->FluxSub[i][j];
					}	
				}				
				
				if (MD->EleET[i][3]>0)
				{
					MD->EleET[i][3] = 0;
				}
				if (MD->EleET[i][5]>0)
				{
					MD->EleET[i][5] = 0;
				}
				
				if (MD->EleViR[i]<0)
				{
					MD->EleViR[i] = 0;
					MD->Recharge[i] = MD->EleViR[i];
				}
				if (MD->Recharge[i]<0&&MD->infil_mode[i]==1)
				{
					MD->Recharge[i] = 0;
				}
			}
			//-----------Finish saturated water------------------//
        
        //-----------Saltwater------------------//
        if (MD->Total_Saltwater_out[i]>0 && MD->Total_Saltwater_out[i]/UNIT_C>MD->DummyY[i + 3*MD->NumEle]*MD->Ele[i].Porosity && MD->DummyY[i + 3*MD->NumEle]>=0)
        {
            for (j=0;j<3;j++)
            {
                inabr = MD->Ele[i].nabr[j] - 1;
                
                if (MD->FluxSalt[i][j]>0)
                {
                    MD->FluxSalt[i][j] = ((MD->FluxSalt[i][j]/MD->Ele[i].area)/MD->Total_Saltwater_out[i])*MD->DummyY[i + 3*MD->NumEle]*MD->Ele[i].Porosity*UNIT_C*MD->Ele[i].area;
                    if (inabr>=0)
                        MD->FluxSalt[inabr][MD->Ele[inabr].Nabr_index[i]] = -MD->FluxSalt[i][j];
                }
            }
        }
        else if (MD->DummyY[i + 3*MD->NumEle]<0)
        {
            for (j=0;j<3;j++)
            {
                inabr = MD->Ele[i].nabr[j] - 1;
                
                if (MD->FluxSalt[i][j]>0)
                {
                    MD->FluxSalt[i][j] = 0;
                    if (inabr>=0)
                        MD->FluxSalt[inabr][MD->Ele[inabr].Nabr_index[i]] = -MD->FluxSalt[i][j];
                }
            }
        }
        //-----------Finish saturated water------------------//
	}
}

void ET(Model_Data MD, int i,double t)
{
	double        Delta, Gamma;
	double        Rn, G, T, Vel, RH, VP, P, LAI, zero_dh, cnpy_h, rl,
	                r_a, r_s, alpha_r, f_r, eta_s, beta_s, gamma_s, Rmax,
	                Lambda, P_c, qv, qv_sat, ETp;
	double        ThetaRef, ThetaW;
	double        Avg_Y_Surf, Dif_Y_Surf, Grad_Y_Surf, Avg_Sf, Distance,Dif_Y_Surf1,Dif_Y_Surf2;
	double        Cwr, TotalY_Riv, TotalY_Riv_down, CrossA, CrossAdown, AvgCrossA,
	                Perem, Perem_down, Avg_Rough, Avg_Perem, Avg_Y_Riv,
	                Dif_Y_Riv, Grad_Y_Riv, Wid, Wid_down, Avg_Wid;
	double        Avg_Y_Sub, Dif_Y_Sub, Avg_Ksat, Grad_Y_Sub, nabrAqDepth, AquiferDepth,
	                Deficit, elemSatn, satKfunc, effK, effKnabr, TotalY_Ele,
	                TotalY_Ele_down;
	double		Dif_Elev, Grad_Elev;
	double 		c1_square,c2_square,L,b;
	double       *Y, *DY;
	double		time_convert = 60*60*24/UNIT_C;
	double		h_regolith, Tau0, Critical_ShieldsStress,U_star, Renolds, ShieldsStress, q_star,R,q_b,Total_soil_out, Total_soil_in;
	double		Cf = 0.003;
	double		Cos_slope, Sec_slope;
	double		ETa;
	double		tempET;
	
	
			//----------------Set up aquifer depth----------------------//
			AquiferDepth = (MD->Ele[i].zmax - MD->Ele[i].zmin);
			
			Deficit = AquiferDepth - MD->DummyY[i + 2 * MD->NumEle] - MD->DummyY[i + 3 * MD->NumEle];
			
			if (AquiferDepth < MD->Ele[i].RzD)
			{
				if (AquiferDepth<0)
				{
					MD->Ele[i].RzD = 0;
					AquiferDepth = 0;
				}
				else
				{
					MD->Ele[i].RzD = AquiferDepth;
				}
			}
			//----------------Finish Setting up aquifer depth----------------------//
			
			//----------------Calculation----------------------//
			Rn = Interpolation(&MD->TSD_Rn[MD->Ele[i].Rn - 1], t);
			//G = Interpolation(&MD->TSD_G[MD->Ele[i].G - 1], t);
			G = 0.1 * Rn;
			T = Interpolation(&MD->TSD_Temp[MD->Ele[i].temp - 1], t);
			Vel = Interpolation(&MD->TSD_WindVel[MD->Ele[i].WindVel - 1], t);
			RH = Interpolation(&MD->TSD_Humidity[MD->Ele[i].humidity - 1], t);
			VP = 611.2 * exp(17.67 * T / (T + 243.5)) * RH;
			P = 101.325 * pow(10, 3) * pow((293 - 0.0065 * MD->Ele[i].zmax) / 293, 5.26);
			
			qv = 0.622 * VP / P;
			qv_sat = 0.622 * (VP / RH) / P;
			//P = 101.325 * pow(10, 3) * pow((293 - 0.0065 * MD->Ele[i].zmax) / 293, 5.26);
			//Delta = 2503 * pow(10, 3) * exp(17.27 * T / (T + 237.3)) / (pow(237.3 + T, 2));
			//Gamma = P * 1.0035 * 0.92 / (0.622 * 2441);
			LAI = Interpolation(&MD->TSD_LAI[MD->Ele[i].LC - 1], t);
			/*
			* zero_dh=Interpolation(&MD->TSD_DH[MD->Ele[i].LC-1], t);
			* cnpy_h =
			* zero_dh/(1.1*(0.0000001+log(1+pow(0.007*LAI,0.25))));
			* if(LAI<2.85)	{ rl= 0.0002 + 0.3*cnpy_h*pow(0.07*LAI,0.5);
			* } else { rl= 0.3*cnpy_h*(1-(zero_dh/cnpy_h)); }
			*/
			rl = Interpolation(&MD->TSD_RL[MD->Ele[i].LC - 1], t);
			r_a = 12 * 4.72 * log(MD->Ele[i].windH / rl) / (0.54 * Vel / UNIT_C / 60 + 1) / UNIT_C / 60;

			Gamma = 4 * 0.7 * SIGMA * UNIT_C * R_dry / C_air * pow(T + 273.15, 4) / (P / r_a) + 1;
			Delta = Lv * Lv * 0.622 / R_v / C_air / pow(T + 273.15, 2) * qv_sat;
			ETp = (Rn * Delta + Gamma * (1.2 * Lv * (qv_sat - qv) / r_a)) / (1000.0 * Lv * (Delta + Gamma));
                        MD->ElePET[i] = ETp;
			//MD->EleET[i][2] = MD->pcCal.Et2 * (1 - MD->Ele[i].VegFrac) * (Rn * (1 - MD->Ele[i].Albedo) * Delta + (1.2 * 1003.5 * ((VP / RH) - VP) / r_a)) / (1000.0 * 2441000.0 * (Delta + Gamma));
			// BHATT: MAJOR BUG = AQUIFER DEPTH NOT CALCULATED EARLIER
			if (Deficit <= MD->Ele[i].RzD) 
			{
				elemSatn = 1.0;
			} 
			else 
			{
				elemSatn = ((MD->DummyY[i + MD->NumEle] / (AquiferDepth - MD->DummyY[i + 2 * MD->NumEle])) > 1) ? 1 : ((MD->DummyY[i + MD->NumEle] / (AquiferDepth - MD->DummyY[i + 2 * MD->NumEle])) < 0) ? 0 : 0.5 * (1 - cos(3.14 * (MD->DummyY[i + MD->NumEle] / (AquiferDepth - MD->DummyY[i + 2 * MD->NumEle]))));
			}
			ThetaRef = 0.7 * MD->Soil[(MD->Ele[i].soil - 1)].ThetaS;
			ThetaW = 1.05 * MD->Soil[(MD->Ele[i].soil - 1)].ThetaR;
			beta_s = (elemSatn * MD->Ele[i].Porosity + MD->Soil[(MD->Ele[i].soil - 1)].ThetaR - ThetaW) / (ThetaRef - ThetaW);
			beta_s = (beta_s < 0.0001) ? 0.0001 : (beta_s > 1 ? 1 : beta_s);
			ETa = MD->pcCal.Et2 * (1 - MD->Ele[i].VegFrac) * beta_s * ETp;
			if (LAI > 0.0) 
			{
				Rmax = 5000.0 / (60 * UNIT_C);	/* Unit day_per_m */
				f_r = 1.1 * 1.5 * Rn / (MD->Ele[i].Rs_ref * LAI);
				f_r = f_r < 0 ? 0 : f_r;
				alpha_r = (1 + f_r) / (f_r + (MD->Ele[i].Rmin / Rmax));
				alpha_r = alpha_r > 10000 ? 10000 : alpha_r;
				eta_s = 1 - 0.0016 * (pow((24.85 - T), 2));
				eta_s = eta_s < 0.0001 ? 0.0001 : eta_s;
				gamma_s = 1 / (1 + 0.00025 * (VP / RH - VP));
				gamma_s = (gamma_s < 0.01) ? 0.01 : gamma_s;
				r_s = ((MD->Ele[i].Rmin * alpha_r / (beta_s * LAI * eta_s * gamma_s)) > Rmax) ? Rmax : (MD->Ele[i].Rmin * alpha_r / (beta_s * LAI * eta_s * gamma_s));
				P_c = (1 + Delta / Gamma) / (1 + r_s / r_a + Delta / Gamma);
				tempET = MD->pcCal.Et1 * MD->Ele[i].VegFrac * P_c * (1 - pow(((MD->EleIS[i] + MD->EleSnowCanopy[i] < 0) ? 0 : (MD->EleIS[i] + MD->EleSnowCanopy[i])) / (MD->EleISmax[i] + MD->EleISsnowmax[i]), 1.0 / 2.0)) * ETp;
				tempET = tempET < 0 ? 0 : tempET;
			} 
			else 
			{
				tempET = 0.0;
			}
			
			if (AquiferDepth > 0)
			{
				if (MD->DummyY[i]>ETa)
				{
					MD->EleET[i][1] = ETa;
					MD->EleET[i][2] = 0;
					MD->EleET[i][3] = 0;
				}
				else
				{
					MD->EleET[i][1] = MD->DummyY[i];
					if (Deficit>MD->Ele[i].infD)
					{
						if (MD->DummyY[i + MD->NumEle]*MD->Ele[i].Porosity>ETa-MD->DummyY[i])
						{
							MD->EleET[i][2] = ETa-MD->DummyY[i];
							MD->EleET[i][3] = 0;
						}
						else
						{
							MD->EleET[i][2] = MD->DummyY[i + MD->NumEle]*MD->Ele[i].Porosity;
							MD->EleET[i][3] = 0;
						}
					}
					else
					{
						if (MD->DummyY[i + 2*MD->NumEle]*MD->Ele[i].Porosity>ETa-MD->DummyY[i])
						{
							MD->EleET[i][2] = 0;
							MD->EleET[i][3] = ETa-MD->DummyY[i];
						}
						else
						{
							MD->EleET[i][2] = 0;
							MD->EleET[i][3] = MD->DummyY[i + 2*MD->NumEle]*MD->Ele[i].Porosity;
						}
					}
				}
				
			
				if (Deficit>MD->Ele[i].RzD)
				{
					if (MD->DummyY[i + MD->NumEle]*MD->Ele[i].Porosity>tempET)
					{
						MD->EleET[i][4] = tempET;
						MD->EleET[i][5] = 0;
					}
					else
					{
						MD->EleET[i][4] = MD->DummyY[i + MD->NumEle]*MD->Ele[i].Porosity;
						MD->EleET[i][5] = 0;
					}
				}
				else
				{
					if (MD->DummyY[i + 2*MD->NumEle]*MD->Ele[i].Porosity>tempET)
					{
						MD->EleET[i][4] = 0;
						MD->EleET[i][5] = tempET;
					}
					else
					{
						MD->EleET[i][4] = 0;
						MD->EleET[i][5] = MD->DummyY[i + 2*MD->NumEle]*MD->Ele[i].Porosity;
					}
				}
			}
			else
			{
				if (MD->DummyY[i]>ETa)
				{
					MD->EleET[i][1] = ETa;
					MD->EleET[i][2] = 0;
					MD->EleET[i][3] = 0;
					MD->EleET[i][4] = 0;
					MD->EleET[i][5] = 0;
				}
				else
				{
					MD->EleET[i][1] = MD->DummyY[i];
					MD->EleET[i][2] = 0;
					MD->EleET[i][3] = 0;
					MD->EleET[i][4] = 0;
					MD->EleET[i][5] = 0;
				}
			}	
	MD->Total_Surfwater_out[i] = MD->Total_Surfwater_out[i] + MD->EleET[i][1];
	MD->Total_Unsatwater_out[i] = MD->Total_Unsatwater_out[i] + MD->EleET[i][2] + MD->EleET[i][4];
	MD->Total_Satwater_out[i] = MD->Total_Satwater_out[i] + MD->EleET[i][3] + MD->EleET[i][5];
}

void Infiltration(Model_Data MD, int i)
{
	double        Delta, Gamma;
	double        Rn, G, T, Vel, RH, VP, P, LAI, zero_dh, cnpy_h, rl,
	                r_a, r_s, alpha_r, f_r, eta_s, beta_s, gamma_s, Rmax,
	                Lambda, P_c, qv, qv_sat, ETp;
	double        ThetaRef, ThetaW;
	double        Avg_Y_Surf, Dif_Y_Surf, Grad_Y_Surf, Avg_Sf, Distance,Dif_Y_Surf1,Dif_Y_Surf2;
	double        Cwr, TotalY_Riv, TotalY_Riv_down, CrossA, CrossAdown, AvgCrossA,
	                Perem, Perem_down, Avg_Rough, Avg_Perem, Avg_Y_Riv,
	                Dif_Y_Riv, Grad_Y_Riv, Wid, Wid_down, Avg_Wid;
	double        Avg_Y_Sub, Dif_Y_Sub, Avg_Ksat, Grad_Y_Sub, nabrAqDepth, AquiferDepth,
	                Deficit, elemSatn, satKfunc, effK, effKnabr, TotalY_Ele,
	                TotalY_Ele_down;
	double		Dif_Elev, Grad_Elev;
	double 		c1_square,c2_square,L,b;
	double       *Y, *DY;
	double		time_convert = 60*60*24/UNIT_C;
	double		h_regolith, Tau0, Critical_ShieldsStress,U_star, Renolds, ShieldsStress, q_star,R,q_b,Total_soil_out, Total_soil_in;
	double		Cf = 0.003;
	double		Cos_slope, Sec_slope;
	double		ETa;
	double		temp;
	
	//----------------Set up aquifer depth----------------------//
	AquiferDepth = (MD->Ele[i].zmax - MD->Ele[i].zmin);
	Deficit = AquiferDepth - MD->DummyY[i + 2 * MD->NumEle] - MD->DummyY[i + 3 * MD->NumEle];
			
	if (AquiferDepth < MD->Ele[i].infD)
	{
		if (AquiferDepth<0)
		{
			MD->Ele[i].infD = 0;
			AquiferDepth = 0;
		}
		else
		{
			MD->Ele[i].infD = AquiferDepth;
		}
	}
	//----------------Finish Setting up aquifer depth----------------------//
		
	//----------------Calculation----------------------//
	if (AquiferDepth>0)
	{
		if (Deficit>MD->Ele[i].infD)
		{
			if ((MD->DummyY[i + MD->NumEle])<Deficit)
			{
				MD->infil_mode[i] = 1;
				elemSatn = MD->DummyY[i + MD->NumEle]/Deficit;
				elemSatn = (elemSatn < multF * EPS) ? multF * EPS : elemSatn;
				Avg_Y_Sub = (-(pow(pow(1 / elemSatn, MD->Ele[i].Beta / (MD->Ele[i].Beta - 1)) - 1, 1 / MD->Ele[i].Beta) / MD->Ele[i].Alpha) < MINpsi) ? MINpsi : -(pow(pow(1 / elemSatn, MD->Ele[i].Beta / (MD->Ele[i].Beta - 1)) - 1, 1 / MD->Ele[i].Beta) / MD->Ele[i].Alpha);
							
				TotalY_Ele = Avg_Y_Sub + MD->Ele[i].zmax - Deficit;
				Grad_Y_Sub = (MD->DummyY[i] + MD->Ele[i].zmax - TotalY_Ele) / Deficit;
					
				Grad_Y_Sub = ((MD->DummyY[i] < EPS / 100) && (Grad_Y_Sub > 0)) ? 0 : Grad_Y_Sub;
				satKfunc = pow(elemSatn, 0.5) * pow(-1 + pow(1 - pow(elemSatn, MD->Ele[i].Beta / (MD->Ele[i].Beta - 1)), (MD->Ele[i].Beta - 1) / MD->Ele[i].Beta), 2);
				effK = (MD->Ele[i].Macropore == 1) ? effKV(satKfunc, Grad_Y_Sub, MD->Ele[i].macKsatV, MD->Ele[i].infKsatV, MD->Ele[i].hAreaF) : MD->Ele[i].infKsatV;
				if (Avg_Y_Sub < 0)
				{
					MD->EleViR[i] = effK * Grad_Y_Sub;
				}
				else 
				{
					MD->EleViR[i] = 0;
				}
							
				effK = (MD->Ele[i].Macropore == 1) ? ((MD->DummyY[i + 2 * MD->NumEle] > AquiferDepth - MD->Ele[i].macD) ? effK : MD->Ele[i].KsatV * satKfunc) : MD->Ele[i].KsatV * satKfunc;
				MD->Recharge[i] = (elemSatn == 0.0) ? 0 : (Deficit <= 0) ? 0 : (MD->Ele[i].KsatV * MD->DummyY[i + 2 * MD->NumEle] + effK * Deficit) * (MD->Ele[i].Alpha * Deficit - 2 * pow(-1 + pow(elemSatn, MD->Ele[i].Beta / (-MD->Ele[i].Beta + 1)), 1 / MD->Ele[i].Beta)) / (MD->Ele[i].Alpha * pow(Deficit + MD->DummyY[i + 2 * MD->NumEle], 2));
				MD->Recharge[i] = (MD->Recharge[i] > 0 && MD->DummyY[i + MD->NumEle] <= 0) ? 0 : MD->Recharge[i];
				MD->Recharge[i] = (MD->Recharge[i] < 0 && MD->DummyY[i + 2 * MD->NumEle] <= 0) ? 0 : MD->Recharge[i];
			}
			else if(MD->DummyY[i + MD->NumEle]>=Deficit)
			{
				MD->infil_mode[i] = 1;
				elemSatn = 1;
				MD->EleViR[i] = 0;
							
				effK = (MD->Ele[i].Macropore == 1) ? ((MD->DummyY[i + 2 * MD->NumEle] > AquiferDepth - MD->Ele[i].macD) ? effK : MD->Ele[i].KsatV * satKfunc) : MD->Ele[i].KsatV * satKfunc;
				MD->Recharge[i] = (elemSatn == 0.0) ? 0 : (Deficit <= 0) ? 0 : (MD->Ele[i].KsatV * MD->DummyY[i + 2 * MD->NumEle] + effK * Deficit) * (MD->Ele[i].Alpha * Deficit - 2 * pow(-1 + pow(elemSatn, MD->Ele[i].Beta / (-MD->Ele[i].Beta + 1)), 1 / MD->Ele[i].Beta)) / (MD->Ele[i].Alpha * pow(Deficit + MD->DummyY[i + 2 * MD->NumEle], 2));
				MD->Recharge[i] = (MD->Recharge[i] > 0 && MD->DummyY[i + MD->NumEle] <= 0) ? 0 : MD->Recharge[i];
				MD->Recharge[i] = (MD->Recharge[i] < 0 && MD->DummyY[i + 2 * MD->NumEle] <= 0) ? 0 : MD->Recharge[i];
			}	
		}
		else
		{
			MD->infil_mode[i] = -1;
			Grad_Y_Sub = (MD->DummyY[i] + MD->Ele[i].zmax - (MD->DummyY[i + 2 * MD->NumEle] + MD->DummyY[i + 3 * MD->NumEle] + MD->Ele[i].zmin)) / AquiferDepth;
			Grad_Y_Sub = ((MD->DummyY[i] < EPS / 100) && (Grad_Y_Sub > 0)) ? 0 : Grad_Y_Sub;
			elemSatn = 1.0;
			satKfunc = pow(elemSatn, 0.5) * pow(-1 + pow(1 - pow(elemSatn, MD->Ele[i].Beta / (MD->Ele[i].Beta - 1)), (MD->Ele[i].Beta - 1) / MD->Ele[i].Beta), 2);
			effK = (MD->Ele[i].Macropore == 1) ? effKV(satKfunc, Grad_Y_Sub, MD->Ele[i].macKsatV, MD->Ele[i].infKsatV, MD->Ele[i].hAreaF) : MD->Ele[i].infKsatV;
						
			if(Grad_Y_Sub<0 && MD->DummyY[i + 2 * MD->NumEle]<=0)
            {
                MD->EleViR[i] = 0;
                MD->Recharge[i] = 0;
            }
            else
            {
                MD->EleViR[i] = effK * Grad_Y_Sub;
                MD->Recharge[i] = MD->EleViR[i];
            }
            
		}
	}
	else
	{
		MD->EleViR[i] = 0;
		MD->Recharge[i] = 0;
	}
	//----------------Calculation----------------------//
	if (MD->infil_mode[i] == -1)
	{
		if (MD->EleViR[i]>0)
		{
			MD->Total_Surfwater_out[i] = MD->Total_Surfwater_out[i] + MD->EleViR[i];
		}
		
		if (MD->EleViR[i]<0)
		{
			MD->Total_Satwater_out[i] = MD->Total_Satwater_out[i] - MD->EleViR[i];
		}
	}
	
	if (MD->infil_mode[i] == 1)
	{
		if (MD->EleViR[i]>0)
		{
			MD->Total_Surfwater_out[i] = MD->Total_Surfwater_out[i] + MD->EleViR[i];
		}
		
		if (MD->Recharge[i]>0)
		{
			MD->Total_Unsatwater_out[i] = MD->Total_Unsatwater_out[i] + MD->Recharge[i];
		}
		
		if (MD->Recharge[i]<0)
		{
			MD->Total_Satwater_out[i] = MD->Total_Satwater_out[i] - MD->Recharge[i];
		}
	}
}

void HydroSurfLateral(Model_Data MD, int i,int j,int inabr, double t)
{
	//int             i, j, k, inabr,boundary_index,q,is_outlet,out_i,out_j;
	double        Delta, Gamma;
	double        Rn, G, T, Vel, RH, VP, P, LAI, zero_dh, cnpy_h, rl,
	                r_a, r_s, alpha_r, f_r, eta_s, beta_s, gamma_s, Rmax,
	                Lambda, P_c, qv, qv_sat, ETp;
	double        ThetaRef, ThetaW;
	double        Avg_Y_Surf, Dif_Y_Surf, Grad_Y_Surf, Avg_Sf, Distance,Dif_Y_Surf1,Dif_Y_Surf2;
	double        Cwr, TotalY_Riv, TotalY_Riv_down, CrossA, CrossAdown, AvgCrossA,
	                Perem, Perem_down, Avg_Rough, Avg_Perem, Avg_Y_Riv,
	                Dif_Y_Riv, Grad_Y_Riv, Wid, Wid_down, Avg_Wid;
	double        Avg_Y_Sub, Dif_Y_Sub, Avg_Ksat, Grad_Y_Sub, nabrAqDepth, AquiferDepth,
	                Deficit, elemSatn, satKfunc, effK, effKnabr, TotalY_Ele,
	                TotalY_Ele_down;
	double		Dif_Elev, Grad_Elev;
	double 		c1_square,c2_square,L,b;
	double       *Y, *DY;
	double		time_convert = 60*60*24/UNIT_C;
	double		h_regolith, Tau0, Critical_ShieldsStress,U_star, Renolds, ShieldsStress, q_star,R,q_b,Total_soil_out, Total_soil_in;
	double		Cf = 0.003;
	double		Cos_slope, Sec_slope;
	double 		lake_level_adjust=0;
	double		spill_flow = 0;
	double		Orifice_flow = 0;
	double		total_sea_level = 0;
	
	//----------------Calculation----------------------//
	if (inabr>=0) 
	{
		Dif_Y_Surf = (MD->SurfMode == 1) ? (MD->Ele[i].zmax - MD->Ele[inabr].zmax) : (MD->DummyY[i] + MD->Ele[i].zmax) - (MD->DummyY[inabr] + MD->Ele[inabr].zmax);
		
		Avg_Y_Surf = avgY(Dif_Y_Surf, MD->DummyY[i], MD->DummyY[inabr]);
		Distance = sqrt(pow((MD->Ele[i].x - MD->Ele[inabr].x), 2) + pow((MD->Ele[i].y - MD->Ele[inabr].y), 2));
        MD->Distance[i][j] = Distance;
		Grad_Y_Surf = Dif_Y_Surf / Distance;
		Avg_Sf = 0.5*(sqrt(pow(MD->Ele[i].dhBYdx, 2) + pow(MD->Ele[i].dhBYdy, 2))+sqrt(pow(MD->Ele[inabr].dhBYdx, 2) + pow(MD->Ele[inabr].dhBYdy, 2)));//?? Xuan Weighting needed
		//updated in 2.2
		Avg_Sf = (MD->SurfMode == 1) ? (Grad_Y_Surf > 0 ? Grad_Y_Surf : EPS / pow(10.0, 6)) : (Avg_Sf > EPS / pow(10.0, 6)) ? Avg_Sf : EPS / pow(10.0, 6);
		/* Weighting needed */
		Avg_Rough = 0.5 * (MD->Ele[i].Rough + MD->Ele[inabr].Rough);
		CrossA = Avg_Y_Surf * MD->Ele[i].edge[j];
        MD->CrossA_surf[i][j] = CrossA;
		MD->FluxSurf[i][j] = CrossA * pow(Avg_Y_Surf, 2.0 / 3.0) * Grad_Y_Surf / (sqrt(fabs(Avg_Sf)) * Avg_Rough);
	}
	else
	{
		if (MD->Ele[i].BC[j] == 0) 
		{
			MD->FluxSurf[i][j] = 0;
            MD->Distance[i][j] = -99999;
		}
		else if (MD->Ele[i].BC[j] == 1) {	
			/* Note: ideally
			* different boundary
			* conditions need to be
			* incorporated	for surf
			* and subsurf
			* respectively */
			/*
			* Note: the formulation assumes only
			* dirichlet TS right now
			*/
			MD->FluxSurf[i][j] = 0;	/* Note the assumption
									 * here is no flow for
									 * surface */
			Dif_Y_Sub = (MD->DummyY[i + 2 * MD->NumEle] + MD->Ele[i].zmin) - Interpolation(&MD->TSD_EleBC[(MD->Ele[i].BC[j]) - 1], t);
			Avg_Y_Sub = avgY(Dif_Y_Sub, MD->DummyY[i + 2 * MD->NumEle], (Interpolation(&MD->TSD_EleBC[(MD->Ele[i].BC[j]) - 1], t) - MD->Ele[i].zmin));
			//Avg_Y_Sub = (MD->DummyY[i + 2 * MD->NumEle] + (Interpolation(&MD->TSD_EleBC[(MD->Ele[i].BC[j]) - 1], t) - MD->Ele[i].zmin)) / 2;
			/*
			 * Minimum Distance from circumcenter
			 * to the edge of the triangle on
			 * which BDD. condition is defined
			 */
			Distance = sqrt(pow(MD->Ele[i].edge[0] * MD->Ele[i].edge[1] * MD->Ele[i].edge[2] / (4 * MD->Ele[i].area), 2) - pow(MD->Ele[i].edge[j] / 2, 2));
			effK = effKH(MD->Ele[i].Macropore, MD->DummyY[i + 2 * MD->NumEle], AquiferDepth, MD->Ele[i].macD, MD->Ele[i].macKsatH, MD->Ele[i].vAreaF, MD->Ele[i].KsatH);
			Avg_Ksat = effK;
			Grad_Y_Sub = Dif_Y_Sub / Distance;
			MD->FluxSub[i][j] = Avg_Ksat * Grad_Y_Sub * Avg_Y_Sub * MD->Ele[i].edge[j];
		} 
		else if (MD->Ele[i].BC[j] == 2)
		{	
			/* Neumann BC (Note:
			 * MD->Ele[i].BC[j] value
			 * have to be = 2+(index of
			 * neumann boundary TS) 
			 */
			MD->FluxSurf[i][j] = Interpolation(&MD->TSD_EleBC[(MD->Ele[i].BC[j]) - 1], t);
			MD->FluxSub[i][j] = Interpolation(&MD->TSD_EleBC[(-MD->Ele[i].BC[j]) - 1], t);
		}
		else if (MD->Ele[i].BC[j] == 3)
		{
			Dif_Y_Surf = MD->DummyY[i];
			Avg_Y_Surf = MD->DummyY[i];
			Distance = sqrt(pow(MD->Ele[i].edge[0] * MD->Ele[i].edge[1] * MD->Ele[i].edge[2] / (4 * MD->Ele[i].area), 2) - pow(MD->Ele[i].edge[j] / 2, 2));
            MD->Distance[i][j] = Distance;
			Grad_Y_Surf = Dif_Y_Surf/Distance;
			Avg_Rough = MD->Ele[i].Rough;
			CrossA = MD->DummyY[i]*MD->Ele[i].edge[j];
            MD->CrossA_surf[i][j] = CrossA;
			MD->FluxSurf[i][j] = sqrt((Grad_Y_Surf>0)?Grad_Y_Surf:0)*CrossA*((Avg_Y_Surf>EPS/ pow(10.0, 6))?pow(Avg_Y_Surf,2.0/3.0):0)/Avg_Rough; 														
		}
		else if (MD->Ele[i].BC[j] == 4)
		{
			lake_level_adjust  = MD->LakeSurfElev[MD->Land_lake_index[i]-1]-MD->LakeInitialElev[MD->Land_lake_index[i]-1];
			Dif_Y_Surf = MD->DummyY[i]+MD->Ele[i].zmax+MD->BankElevAdjust[MD->Land_lake_index[i]-1] - (MD->Ele[i].zmax + lake_level_adjust);
			if (lake_level_adjust>=0)
			{
				Avg_Y_Surf = avgY(Dif_Y_Surf, MD->DummyY[i], lake_level_adjust);
			}
			else
			{
				Avg_Y_Surf = avgY(Dif_Y_Surf, MD->DummyY[i], MD->DummyY[i]);
			}
			
			Distance = sqrt(pow(MD->Ele[i].edge[0] * MD->Ele[i].edge[1] * MD->Ele[i].edge[2] / (4 * MD->Ele[i].area), 2) - pow(MD->Ele[i].edge[j] / 2, 2));
            MD->Distance[i][j] = Distance;
			Grad_Y_Surf = Dif_Y_Surf/Distance;
			Avg_Rough = MD->Ele[i].Rough;
			CrossA = Avg_Y_Surf*MD->Ele[i].edge[j];
            MD->CrossA_surf[i][j] = CrossA;
			MD->FluxSurf[i][j] = sqrt(fabs(Grad_Y_Surf))*CrossA*((Avg_Y_Surf>EPS/ pow(10.0, 6))?pow(Avg_Y_Surf,2.0/3.0):0)/Avg_Rough;
            if (Dif_Y_Surf<0)
            {
                MD->FluxSurf[i][j] = -MD->FluxSurf[i][j];
            }
		}
        else if (MD->Ele[i].BC[j] == 5)
        {
            lake_level_adjust  = MD->LakeSurfElev[MD->Land_lake_index[i]-1] - MD->LakeInitialElev[MD->Land_lake_index[i]-1] + MD->LakeOutletDecline[MD->Land_lake_index[i]-1];
            Dif_Y_Surf = MD->DummyY[i]+MD->Ele[i].zmax - (MD->Ele[i].zmax + lake_level_adjust);
            if (lake_level_adjust>=0)
            {
                Avg_Y_Surf = avgY(Dif_Y_Surf, MD->DummyY[i], lake_level_adjust);
            }
            else
            {
                Avg_Y_Surf = avgY(Dif_Y_Surf, MD->DummyY[i], MD->DummyY[i]);
            }
            
            Distance = sqrt(pow(MD->Ele[i].edge[0] * MD->Ele[i].edge[1] * MD->Ele[i].edge[2] / (4 * MD->Ele[i].area), 2) - pow(MD->Ele[i].edge[j] / 2, 2));
            MD->Distance[i][j] = Distance;
            Grad_Y_Surf = Dif_Y_Surf/Distance;
            Avg_Rough = MD->Ele[i].Rough;
            CrossA = Avg_Y_Surf*MD->Ele[i].edge[j];
            MD->CrossA_surf[i][j] = CrossA;
            MD->FluxSurf[i][j] = sqrt(fabs(Grad_Y_Surf))*CrossA*((Avg_Y_Surf>EPS/ pow(10.0, 6))?pow(Avg_Y_Surf,2.0/3.0):0)/Avg_Rough;
            if (Dif_Y_Surf<=0)
            {
                MD->FluxSurf[i][j] = -MD->FluxSurf[i][j];
            }
        }
		else if (MD->Ele[i].BC[j] == 6)
        {
            Distance = sqrt(pow(MD->Ele[i].edge[0] * MD->Ele[i].edge[1] * MD->Ele[i].edge[2] / (4 * MD->Ele[i].area), 2) - pow(MD->Ele[i].edge[j] / 2, 2));
            MD->Distance[i][j] = Distance;
			lake_level_adjust  = MD->LakeSurfElev[MD->Land_lake_index[i]-1] - MD->LakeInitialElev[MD->Land_lake_index[i]-1];
            Dif_Y_Surf = MD->Ele[i].zmax - (MD->Ele[i].zmax + lake_level_adjust);
            if (Dif_Y_Surf<=0)
            {
                Avg_Y_Surf = MD->Ele[i].zmax + lake_level_adjust-MD->LakeBed[MD->Land_lake_index[i]-1];
				spill_flow = 2/3*0.62*sqrt(2*GRAV)*MD->Weir_length[MD->Land_lake_index[i]-1]*pow(lake_level_adjust,3/2);
				Orifice_flow = 0.61*MD->Orifice_height[MD->Land_lake_index[i]-1]*MD->Weir_length[MD->Land_lake_index[i]-1]*sqrt(2*GRAV*(Avg_Y_Surf-MD->Orifice_height[MD->Land_lake_index[i]-1]));
				MD->FluxSurf[i][j] = - spill_flow*UNIT_C - Orifice_flow*UNIT_C;
				if(isnan(MD->FluxSurf[i][j])==1)
				{
					printf("is nan");
				}
            }
            else
            {
                MD->FluxSurf[i][j] = 0;
            }
        }
		
		else if (MD->Ele[i].BC[j] == 7) {	
			/* Boundary condition for Tide and sea level rise*/
			/*
			* Note: the formulation assumes only
			* dirichlet TS right now
			*/
			if (MD->Tide_mode == 1 && MD->Sea_level_rise_mode == 1)
			{
				total_sea_level = Interpolation(&MD->TSD_Tide[0], t) + 0.0003 * t / 43200 - 0.0243;
				Dif_Y_Surf = (MD->DummyY[i] + MD->Ele[i].zmax) - total_sea_level;
				
				Avg_Y_Surf = avgY(Dif_Y_Surf, MD->DummyY[i], total_sea_level-MD->Ele[i].zmax);
				//Avg_Y_Sub = (MD->DummyY[i + 2 * MD->NumEle] + (Interpolation(&MD->TSD_EleBC[(MD->Ele[i].BC[j]) - 1], t) - MD->Ele[i].zmin)) / 2;
				/*
				 * Minimum Distance from circumcenter
				 * to the edge of the triangle on
				 * which BDD. condition is defined
				 */
				Distance = sqrt(pow(MD->Ele[i].edge[0] * MD->Ele[i].edge[1] * MD->Ele[i].edge[2] / (4 * MD->Ele[i].area), 2) - pow(MD->Ele[i].edge[j] / 2, 2));
				MD->Distance[i][j] = Distance;
				Grad_Y_Surf = Dif_Y_Surf/Distance;
				Avg_Rough = MD->Ele[i].Rough;
				CrossA = Avg_Y_Surf*MD->Ele[i].edge[j];
				MD->CrossA_surf[i][j] = CrossA;
				MD->FluxSurf[i][j] = sqrt(fabs(Grad_Y_Surf))*CrossA*((Avg_Y_Surf>EPS/ pow(10.0, 6))?pow(Avg_Y_Surf,2.0/3.0):0)/Avg_Rough;
				if (Dif_Y_Surf<0)
				{
					MD->FluxSurf[i][j] = -MD->FluxSurf[i][j];
				}
			}
			else if (MD->Tide_mode == 1 && MD->Sea_level_rise_mode == 0)
			{
				total_sea_level = Interpolation(&MD->TSD_Tide[0], t);
				Dif_Y_Surf = (MD->DummyY[i] + MD->Ele[i].zmax) - total_sea_level;
				
				Avg_Y_Surf = avgY(Dif_Y_Surf, MD->DummyY[i], total_sea_level-MD->Ele[i].zmax);
				//Avg_Y_Sub = (MD->DummyY[i + 2 * MD->NumEle] + (Interpolation(&MD->TSD_EleBC[(MD->Ele[i].BC[j]) - 1], t) - MD->Ele[i].zmin)) / 2;
				/*
				 * Minimum Distance from circumcenter
				 * to the edge of the triangle on
				 * which BDD. condition is defined
				 */
				Distance = sqrt(pow(MD->Ele[i].edge[0] * MD->Ele[i].edge[1] * MD->Ele[i].edge[2] / (4 * MD->Ele[i].area), 2) - pow(MD->Ele[i].edge[j] / 2, 2));
				MD->Distance[i][j] = Distance;
				Grad_Y_Surf = Dif_Y_Surf/Distance;
				Avg_Rough = MD->Ele[i].Rough;
				CrossA = Avg_Y_Surf*MD->Ele[i].edge[j];
				MD->CrossA_surf[i][j] = CrossA;
				MD->FluxSurf[i][j] = sqrt(fabs(Grad_Y_Surf))*CrossA*((Avg_Y_Surf>EPS/ pow(10.0, 6))?pow(Avg_Y_Surf,2.0/3.0):0)/Avg_Rough;
				if (Dif_Y_Surf<0)
				{
					MD->FluxSurf[i][j] = -MD->FluxSurf[i][j];
				}
			}
			else if (MD->Tide_mode == 0 && MD->Sea_level_rise_mode == 1)
			{
				total_sea_level = 0.0003 * t / 43200 - 0.0243;
				Dif_Y_Surf = (MD->DummyY[i] + MD->Ele[i].zmax) - total_sea_level;
				
				Avg_Y_Surf = avgY(Dif_Y_Surf, MD->DummyY[i], total_sea_level-MD->Ele[i].zmax);
				//Avg_Y_Sub = (MD->DummyY[i + 2 * MD->NumEle] + (Interpolation(&MD->TSD_EleBC[(MD->Ele[i].BC[j]) - 1], t) - MD->Ele[i].zmin)) / 2;
				/*
				 * Minimum Distance from circumcenter
				 * to the edge of the triangle on
				 * which BDD. condition is defined
				 */
				Distance = sqrt(pow(MD->Ele[i].edge[0] * MD->Ele[i].edge[1] * MD->Ele[i].edge[2] / (4 * MD->Ele[i].area), 2) - pow(MD->Ele[i].edge[j] / 2, 2));
				MD->Distance[i][j] = Distance;
				Grad_Y_Surf = Dif_Y_Surf/Distance;
				Avg_Rough = MD->Ele[i].Rough;
				CrossA = Avg_Y_Surf*MD->Ele[i].edge[j];
				MD->CrossA_surf[i][j] = CrossA;
				MD->FluxSurf[i][j] = sqrt(fabs(Grad_Y_Surf))*CrossA*((Avg_Y_Surf>EPS/ pow(10.0, 6))?pow(Avg_Y_Surf,2.0/3.0):0)/Avg_Rough;
				if (Dif_Y_Surf<0)
				{
					MD->FluxSurf[i][j] = -MD->FluxSurf[i][j];
				}
			}
			else
			{
				MD->FluxSurf[i][j] = 0;
				MD->Distance[i][j] = -99999;
			}	
		} 
		else
		{
			MD->FluxSurf[i][j] = 0;
            MD->Distance[i][j] = -99999;
		}
	}
	//----------------Finish Calculation----------------------//
	if (MD->FluxSurf[i][j]>0)
	{
		MD->Total_Surfwater_out[i] = MD->Total_Surfwater_out[i] + MD->FluxSurf[i][j]/MD->Ele[i].area;
	}
}


void HydroSubLateral(Model_Data MD, int i,int j,int inabr, double t)
{
	//int             i, j, k, inabr,boundary_index,q,is_outlet,out_i,out_j;
	double        Avg_Y_Surf, Dif_Y_Surf, Grad_Y_Surf, Avg_Sf, Distance,Dif_Y_Surf1,Dif_Y_Surf2;
	double        Cwr, TotalY_Riv, TotalY_Riv_down, CrossA, CrossAdown, AvgCrossA,
	                Perem, Perem_down, Avg_Rough, Avg_Perem, Avg_Y_Riv,
	                Dif_Y_Riv, Grad_Y_Riv, Wid, Wid_down, Avg_Wid;
	double        Avg_Y_Sub, Dif_Y_Sub, Avg_Ksat, Grad_Y_Sub, nabrAqDepth, AquiferDepth,
	                Deficit, elemSatn, satKfunc, effK, effKnabr, TotalY_Ele,
	                TotalY_Ele_down;
	double		Dif_Elev, Grad_Elev;
	double 		c1_square,c2_square,L,b;
	double       *Y, *DY;
	double		time_convert = 60*60*24/UNIT_C;
	double		h_regolith, Tau0, Critical_ShieldsStress,U_star, Renolds, ShieldsStress, q_star,R,q_b,Total_soil_out, Total_soil_in;
	double		Cf = 0.003;
	double		Cos_slope, Sec_slope;
	double		lake_GW_level_adjust=0;
	double		total_sea_level = 0;
	

	if (inabr>=0) 
	{
		//----------------Set up aquifer depth----------------------//
		AquiferDepth = (MD->Ele[i].zmax - MD->Ele[i].zmin);
		nabrAqDepth = (MD->Ele[inabr].zmax - MD->Ele[inabr].zmin);
		if (AquiferDepth < MD->Ele[i].macD)
		{
			if (AquiferDepth<0)
			{
				MD->Ele[i].macD = 0;
				AquiferDepth = 0;
			}
			else
			{
				MD->Ele[i].macD = AquiferDepth;
			}
		}
			
		if (nabrAqDepth < MD->Ele[inabr].macD)
		{
			if (nabrAqDepth<0)
			{
				MD->Ele[inabr].macD = 0;
				nabrAqDepth = 0;
			}
			else
			{
				MD->Ele[inabr].macD = nabrAqDepth;
			}
		}
		//----------------Finish Setting up aquifer depth----------------------//
			
		//----------------Calculation of sub fresh water----------------------//
			Dif_Y_Sub = (MD->DummyY[i + 2 * MD->NumEle] + MD->DummyY[i + 3 * MD->NumEle] + MD->Ele[i].zmin) - (MD->DummyY[inabr + 2 * MD->NumEle] + MD->DummyY[inabr + 3 * MD->NumEle] + MD->Ele[inabr].zmin);
			Avg_Y_Sub = avgY(Dif_Y_Sub, MD->DummyY[i + 2 * MD->NumEle], MD->DummyY[inabr + 2 * MD->NumEle]);
			Distance = sqrt(pow((MD->Ele[i].x - MD->Ele[inabr].x), 2) + pow((MD->Ele[i].y - MD->Ele[inabr].y), 2));
			Grad_Y_Sub = Dif_Y_Sub / Distance;
			/* take care of macropore effect */
			effK = effKH(MD->Ele[i].Macropore, MD->DummyY[i + 2 * MD->NumEle], AquiferDepth, MD->Ele[i].macD, MD->Ele[i].macKsatH, MD->Ele[i].vAreaF, MD->Ele[i].KsatH);
			effKnabr = effKH(MD->Ele[inabr].Macropore, MD->DummyY[inabr + 2 * MD->NumEle], nabrAqDepth, MD->Ele[inabr].macD, MD->Ele[inabr].macKsatH, MD->Ele[inabr].vAreaF, MD->Ele[inabr].KsatH);
			Avg_Ksat = 0.5 * (effK + effKnabr);
            MD->CrossA_sub[i][j] = Avg_Y_Sub * MD->Ele[i].edge[j];
			/* groundwater flow modeled by Darcy's law */
			MD->FluxSub[i][j] = Avg_Ksat * Grad_Y_Sub * Avg_Y_Sub * MD->Ele[i].edge[j];
			MD->Seepage[i][j] = 0;
			//----------------Finish Calculation----------------------//
        
            //----------------Calculation of sub saltwater----------------------//
            Dif_Y_Sub = (1000/1025*MD->DummyY[i + 2 * MD->NumEle] + MD->DummyY[i + 3 * MD->NumEle] + MD->Ele[i].zmin) - (1000/1025*MD->DummyY[inabr + 2 * MD->NumEle] + MD->DummyY[inabr + 3 * MD->NumEle] + MD->Ele[inabr].zmin);
            Avg_Y_Sub = avgY(Dif_Y_Sub, MD->DummyY[i + 3 * MD->NumEle], MD->DummyY[inabr + 3 * MD->NumEle]);
            Grad_Y_Sub = Dif_Y_Sub / Distance;
            /* Saltwater flow */
            MD->FluxSalt[i][j] = Avg_Ksat * Grad_Y_Sub * Avg_Y_Sub * MD->Ele[i].edge[j];
            //----------------Finish Calculation----------------------//
	}
	else if (MD->Ele[i].BC[j] == 0)
	{
		MD->FluxSub[i][j] = 0;
		MD->Seepage[i][j] = 0;
	}
	else if (MD->Ele[i].BC[j] == 1) {	
		/* Note: ideally
		* different boundary
		* conditions need to be
		* incorporated	for surf
		* and subsurf
		* respectively */
		/*
		* Note: the formulation assumes only
		* dirichlet TS right now
		*/
		MD->FluxSurf[i][j] = 0;	/* Note the assumption
								 * here is no flow for
								 * surface */
		Dif_Y_Sub = (MD->DummyY[i + 2 * MD->NumEle] + MD->Ele[i].zmin) - Interpolation(&MD->TSD_EleBC[(MD->Ele[i].BC[j]) - 1], t);
		Avg_Y_Sub = avgY(Dif_Y_Sub, MD->DummyY[i + 2 * MD->NumEle], (Interpolation(&MD->TSD_EleBC[(MD->Ele[i].BC[j]) - 1], t) - MD->Ele[i].zmin));
		//Avg_Y_Sub = (MD->DummyY[i + 2 * MD->NumEle] + (Interpolation(&MD->TSD_EleBC[(MD->Ele[i].BC[j]) - 1], t) - MD->Ele[i].zmin)) / 2;
		/*
		 * Minimum Distance from circumcenter
		 * to the edge of the triangle on
		 * which BDD. condition is defined
		 */
		Distance = sqrt(pow(MD->Ele[i].edge[0] * MD->Ele[i].edge[1] * MD->Ele[i].edge[2] / (4 * MD->Ele[i].area), 2) - pow(MD->Ele[i].edge[j] / 2, 2));
		effK = effKH(MD->Ele[i].Macropore, MD->DummyY[i + 2 * MD->NumEle], AquiferDepth, MD->Ele[i].macD, MD->Ele[i].macKsatH, MD->Ele[i].vAreaF, MD->Ele[i].KsatH);
		Avg_Ksat = effK;
		Grad_Y_Sub = Dif_Y_Sub / Distance;
		MD->FluxSub[i][j] = Avg_Ksat * Grad_Y_Sub * Avg_Y_Sub * MD->Ele[i].edge[j];
	} 
	else if (MD->Ele[i].BC[j] == 2)
	{	
		/* Neumann BC (Note:
		 * MD->Ele[i].BC[j] value
		 * have to be = 2+(index of
		 * neumann boundary TS) 
		 */
		MD->FluxSurf[i][j] = Interpolation(&MD->TSD_EleBC[(MD->Ele[i].BC[j]) - 1], t);
		MD->FluxSub[i][j] = Interpolation(&MD->TSD_EleBC[(-MD->Ele[i].BC[j]) - 1], t);
	}
	else if (MD->Ele[i].BC[j] == 3)
	{
		//----------------Calculation----------------------//
		Dif_Y_Sub = (MD->DummyY[i + 2 * MD->NumEle]  + MD->DummyY[i  + 3 * MD->NumEle] + MD->Ele[i].zmin) - MD->Ele[i].zmax;									
		Avg_Y_Sub = avgY(Dif_Y_Sub, MD->DummyY[i + 2 * MD->NumEle], MD->DummyY[i + 2 * MD->NumEle]);
		Distance = sqrt(pow(MD->Ele[i].edge[0] * MD->Ele[i].edge[1] * MD->Ele[i].edge[2] / (4 * MD->Ele[i].area), 2) - pow(MD->Ele[i].edge[j] / 2, 2));
		Grad_Y_Sub = Dif_Y_Sub / Distance;
		/* take care of macropore effect */
		effK = effKH(MD->Ele[i].Macropore, MD->DummyY[i + 2 * MD->NumEle], AquiferDepth, MD->Ele[i].macD, MD->Ele[i].macKsatH, MD->Ele[i].vAreaF, MD->Ele[i].KsatH);
		Avg_Ksat = effK;
        MD->CrossA_sub[i][j] = Avg_Y_Sub * MD->Ele[i].edge[j];
		/* groundwater flow modeled by Darcy's law */
		MD->FluxSub[i][j] = Avg_Ksat * Grad_Y_Sub * Avg_Y_Sub * MD->Ele[i].edge[j];
		MD->Seepage[i][j] = 0;
	}
	else if (MD->Ele[i].BC[j] == 4)
	{
		lake_GW_level_adjust = MD->LakeGWElev[MD->Land_lake_index[i]-1]-MD->LakeInitialElev[MD->Land_lake_index[i]-1];
		Dif_Y_Sub = (MD->DummyY[i + 2 * MD->NumEle] + MD->Ele[i].zmin) - (MD->Ele[i].zmax + lake_GW_level_adjust);
		Avg_Y_Sub = avgY(Dif_Y_Sub, MD->DummyY[i + 2 * MD->NumEle], (MD->Ele[i].zmax + lake_GW_level_adjust - MD->Ele[i].zmin));
		Grad_Y_Sub = Dif_Y_Sub / Distance;
		/* take care of macropore effect */
		// effK = effKH(MD->Ele[i].Macropore, MD->DummyY[i + 2 * MD->NumEle], AquiferDepth, MD->Ele[i].macD, MD->Ele[i].macKsatH, MD->Ele[i].vAreaF, MD->Ele[i].KsatH);
		// effKnabr = effKH(MD->Ele[inabr].Macropore, MD->DummyY[inabr + 2 * MD->NumEle], nabrAqDepth, MD->Ele[inabr].macD, MD->Ele[inabr].macKsatH, MD->Ele[inabr].vAreaF, MD->Ele[inabr].KsatH);
		/*
			* It should be weighted average. However,
			* there is an ambiguity about distance used
		*/
		// Avg_Ksat = 0.5 * (effK + effKnabr);
		Avg_Ksat = MD->Ele[i].KsatH;
        MD->CrossA_sub[i][j] = Avg_Y_Sub * MD->Ele[i].edge[j];
		/* groundwater flow modeled by Darcy's law */
		MD->FluxSub[i][j] = Avg_Ksat * Grad_Y_Sub * Avg_Y_Sub * MD->Ele[i].edge[j];
		MD->Seepage[i][j] = 0;
	}
    else if (MD->Ele[i].BC[j] == 5)
    {
        lake_GW_level_adjust = MD->LakeGWElev[MD->Land_lake_index[i]-1]-MD->LakeInitialElev[MD->Land_lake_index[i]-1]+MD->LakeOutletDecline[MD->Land_lake_index[i]-1];
        Dif_Y_Sub = (MD->DummyY[i + 2 * MD->NumEle] + MD->Ele[i].zmin) - (MD->Ele[i].zmax + lake_GW_level_adjust);
        Avg_Y_Sub = avgY(Dif_Y_Sub, MD->DummyY[i + 2 * MD->NumEle], (MD->Ele[i].zmax + lake_GW_level_adjust - MD->Ele[i].zmin));
        Grad_Y_Sub = Dif_Y_Sub / Distance;
        /* take care of macropore effect */
        // effK = effKH(MD->Ele[i].Macropore, MD->DummyY[i + 2 * MD->NumEle], AquiferDepth, MD->Ele[i].macD, MD->Ele[i].macKsatH, MD->Ele[i].vAreaF, MD->Ele[i].KsatH);
        // effKnabr = effKH(MD->Ele[inabr].Macropore, MD->DummyY[inabr + 2 * MD->NumEle], nabrAqDepth, MD->Ele[inabr].macD, MD->Ele[inabr].macKsatH, MD->Ele[inabr].vAreaF, MD->Ele[inabr].KsatH);
        /*
         * It should be weighted average. However,
         * there is an ambiguity about distance used
         */
        // Avg_Ksat = 0.5 * (effK + effKnabr);
        Avg_Ksat = MD->Ele[i].KsatH;
        MD->CrossA_sub[i][j] = Avg_Y_Sub * MD->Ele[i].edge[j];
        /* groundwater flow modeled by Darcy's law */
        MD->FluxSub[i][j] = Avg_Ksat * Grad_Y_Sub * Avg_Y_Sub * MD->Ele[i].edge[j];
		MD->Seepage[i][j] = 0;
    }
	else if (MD->Ele[i].BC[j] == 6)
    {
        lake_GW_level_adjust = MD->LakeGWElev[MD->Land_lake_index[i]-1]-MD->LakeInitialElev[MD->Land_lake_index[i]-1]+MD->LakeOutletDecline[MD->Land_lake_index[i]-1];
        Dif_Y_Sub = (MD->DummyY[i + 2 * MD->NumEle] + MD->Ele[i].zmin) - (MD->Ele[i].zmax + lake_GW_level_adjust);
        Avg_Y_Sub = avgY(Dif_Y_Sub, MD->DummyY[i + 2 * MD->NumEle], (MD->Ele[i].zmax + lake_GW_level_adjust - MD->Ele[i].zmin));
        Grad_Y_Sub = Dif_Y_Sub / Distance;
        /* take care of macropore effect */
        // effK = effKH(MD->Ele[i].Macropore, MD->DummyY[i + 2 * MD->NumEle], AquiferDepth, MD->Ele[i].macD, MD->Ele[i].macKsatH, MD->Ele[i].vAreaF, MD->Ele[i].KsatH);
        // effKnabr = effKH(MD->Ele[inabr].Macropore, MD->DummyY[inabr + 2 * MD->NumEle], nabrAqDepth, MD->Ele[inabr].macD, MD->Ele[inabr].macKsatH, MD->Ele[inabr].vAreaF, MD->Ele[inabr].KsatH);
        /*
         * It should be weighted average. However,
         * there is an ambiguity about distance used
         */
        // Avg_Ksat = 0.5 * (effK + effKnabr);
        Avg_Ksat = MD->Ele[i].KsatH;
        MD->CrossA_sub[i][j] = Avg_Y_Sub * MD->Ele[i].edge[j];
        /* groundwater flow modeled by Darcy's law */
        MD->FluxSub[i][j] = Avg_Ksat * Grad_Y_Sub * Avg_Y_Sub * MD->Ele[i].edge[j];
		MD->Seepage[i][j] = 0;
    }
	else if (MD->Ele[i].BC[j] == 7) 
	{	
		/* Boundary condition for Tide and sea level rise*/
		/*
		* Note: the formulation assumes only
		* dirichlet TS right now
		*/
		if (MD->Tide_mode == 1 && MD->Sea_level_rise_mode == 1)
		{
			total_sea_level = Interpolation(&MD->TSD_Tide[0], t) + 0.0003 * t / 43200 - 0.0243;
			Dif_Y_Sub = (MD->DummyY[i  + 2 * MD->NumEle] + MD->DummyY[i  + 3 * MD->NumEle] + MD->Ele[i].zmin) - total_sea_level;
					
			Avg_Y_Sub = avgY(Dif_Y_Sub, MD->DummyY[i  + 2 * MD->NumEle], MD->DummyY[i  + 2 * MD->NumEle]);
			//Avg_Y_Sub = (MD->DummyY[i + 2 * MD->NumEle] + (Interpolation(&MD->TSD_EleBC[(MD->Ele[i].BC[j]) - 1], t) - MD->Ele[i].zmin)) / 2;
			/*
			 * Minimum Distance from circumcenter
			 * to the edge of the triangle on
			 * which BDD. condition is defined
			 */
			Distance = sqrt(pow(MD->Ele[i].edge[0] * MD->Ele[i].edge[1] * MD->Ele[i].edge[2] / (4 * MD->Ele[i].area), 2) - pow(MD->Ele[i].edge[j] / 2, 2));
			MD->Distance[i][j] = Distance;
			Grad_Y_Sub = Dif_Y_Sub/Distance;
			Avg_Ksat = MD->Ele[i].KsatH;
			MD->CrossA_sub[i][j] = Avg_Y_Sub * MD->Ele[i].edge[j];
			/* groundwater flow modeled by Darcy's law */
			MD->FluxSub[i][j] = Avg_Ksat * Grad_Y_Sub * Avg_Y_Sub * MD->Ele[i].edge[j];
            if(MD->FluxSub[i][j]<0)
            {
                MD->FluxSub[i][j]=0;
            }
            
			MD->Seepage[i][j] = 0;
            /*----------Saltwater intrusion----------------*/
            
            Dif_Y_Sub = (1000/1025*MD->DummyY[i  + 2 * MD->NumEle] + MD->DummyY[i  + 3 * MD->NumEle] + MD->Ele[i].zmin) - total_sea_level;
            
            Avg_Y_Sub = avgY(Dif_Y_Sub, MD->DummyY[i  + 3 * MD->NumEle], total_sea_level-MD->Ele[i].zmin);
            //Avg_Y_Sub = (MD->DummyY[i + 2 * MD->NumEle] + (Interpolation(&MD->TSD_EleBC[(MD->Ele[i].BC[j]) - 1], t) - MD->Ele[i].zmin)) / 2;
            /*
             * Minimum Distance from circumcenter
             * to the edge of the triangle on
             * which BDD. condition is defined
             */
            Grad_Y_Sub = Dif_Y_Sub/Distance;
            Avg_Ksat = MD->Ele[i].KsatH;
            /* groundwater flow modeled by Darcy's law */
            MD->FluxSalt[i][j] = Avg_Ksat * Grad_Y_Sub * Avg_Y_Sub * MD->Ele[i].edge[j];
		}
		else if (MD->Tide_mode == 1 && MD->Sea_level_rise_mode == 0)
		{
			total_sea_level = Interpolation(&MD->TSD_Tide[0], t);
			Dif_Y_Sub = (MD->DummyY[i  + 2 * MD->NumEle] + MD->Ele[i].zmin) - total_sea_level;
					
			Avg_Y_Sub = avgY(Dif_Y_Sub, MD->DummyY[i  + 2 * MD->NumEle], MD->DummyY[i  + 2 * MD->NumEle]);
			//Avg_Y_Sub = (MD->DummyY[i + 2 * MD->NumEle] + (Interpolation(&MD->TSD_EleBC[(MD->Ele[i].BC[j]) - 1], t) - MD->Ele[i].zmin)) / 2;
			/*
			 * Minimum Distance from circumcenter
			 * to the edge of the triangle on
			 * which BDD. condition is defined
			 */
			Distance = sqrt(pow(MD->Ele[i].edge[0] * MD->Ele[i].edge[1] * MD->Ele[i].edge[2] / (4 * MD->Ele[i].area), 2) - pow(MD->Ele[i].edge[j] / 2, 2));
			MD->Distance[i][j] = Distance;
			Grad_Y_Sub = Dif_Y_Sub/Distance;
			Avg_Ksat = MD->Ele[i].KsatH;
			MD->CrossA_sub[i][j] = Avg_Y_Sub * MD->Ele[i].edge[j];
			/* groundwater flow modeled by Darcy's law */
			MD->FluxSub[i][j] = Avg_Ksat * Grad_Y_Sub * Avg_Y_Sub * MD->Ele[i].edge[j];
			MD->Seepage[i][j] = 0;
		}
		else if (MD->Tide_mode == 0 && MD->Sea_level_rise_mode == 1)
		{
			total_sea_level = 0.0003 * t / 43200 - 0.0243;
			Dif_Y_Sub = (MD->DummyY[i  + 2 * MD->NumEle] + MD->Ele[i].zmin) - total_sea_level;
					
			Avg_Y_Sub = avgY(Dif_Y_Sub, MD->DummyY[i  + 2 * MD->NumEle], MD->DummyY[i  + 2 * MD->NumEle]);
			//Avg_Y_Sub = (MD->DummyY[i + 2 * MD->NumEle] + (Interpolation(&MD->TSD_EleBC[(MD->Ele[i].BC[j]) - 1], t) - MD->Ele[i].zmin)) / 2;
			/*
			 * Minimum Distance from circumcenter
			 * to the edge of the triangle on
			 * which BDD. condition is defined
			 */
			Distance = sqrt(pow(MD->Ele[i].edge[0] * MD->Ele[i].edge[1] * MD->Ele[i].edge[2] / (4 * MD->Ele[i].area), 2) - pow(MD->Ele[i].edge[j] / 2, 2));
			MD->Distance[i][j] = Distance;
			Grad_Y_Sub = Dif_Y_Sub/Distance;
			Avg_Ksat = MD->Ele[i].KsatH;
			MD->CrossA_sub[i][j] = Avg_Y_Sub * MD->Ele[i].edge[j];
			/* groundwater flow modeled by Darcy's law */
			MD->FluxSub[i][j] = Avg_Ksat * Grad_Y_Sub * Avg_Y_Sub * MD->Ele[i].edge[j];
			MD->Seepage[i][j] = 0;
		}
		else
		{
			MD->FluxSub[i][j] = 0;
			MD->Seepage[i][j] = 0;
		}
	} 
	else
	{
		MD->FluxSub[i][j] = 0;
		MD->Seepage[i][j] = 0;
	}
		
	if (MD->FluxSub[i][j]>0)
		MD->Total_Satwater_out[i] = MD->Total_Satwater_out[i] + MD->FluxSub[i][j]/MD->Ele[i].area;
	if (MD->Seepage[i][j]>0)
		MD->Total_Surfwater_out[i] = MD->Total_Surfwater_out[i] + MD->Seepage[i][j]/MD->Ele[i].area;
    if (MD->FluxSalt[i][j]>0)
        MD->Total_Saltwater_out[i] = MD->Total_Saltwater_out[i] + MD->FluxSalt[i][j]/MD->Ele[i].area;
}

double* initiation(Model_Data MD,double *DY,double t)
{
	int             i, j, k, inabr;
	for (i = 0; i < MD->NumEle; i++) 
	{
		//MD->DummyY[i] = (Y[i] >= 0) ? Y[i] : 0;
		for(j=0;j<3;j++)
		{
			MD->FluxSub[i][j] = 0;
			MD->FluxSurf[i][j] = 0;
            MD->FluxSalt[i][j] = 0;
			MD->Seepage[i][j] = 0;
            MD->CrossA_surf[i][j] = 0;
            MD->CrossA_sub[i][j] = 0;
            MD->Distance[i][j] = 0;
			
		}
		MD->EleViR[i] = 0;
        MD->EleViR_Salt = 0;
		MD->Recharge[i] = 0;
		MD->EleET[i][1] = 0;
		MD->EleET[i][2] = 0;
		MD->EleET[i][3] = 0;
		MD->EleET[i][4] = 0;
		MD->EleET[i][5] = 0;
		
		
		
        DY[i] = 0;
        DY[i + MD->NumEle] = 0;
        DY[i + 2*MD->NumEle] = 0;
        DY[i + 3*MD->NumEle] = 0;
        
        MD->Total_Surfwater_out[i] = 0;
        MD->Total_Unsatwater_out[i] = 0;
        MD->Total_Satwater_out[i] = 0;
        MD->Total_Saltwater_out[i] = 0;
        MD->infil_mode[i] = 0;
	}
	
	for (i = 0; i < MD->NumEle; i++) 
	{
		if ((MD->SurfMode == 2) && (i < MD->NumEle)) 
		{
			for (j = 0; j < 3; j++) 
			{
				inabr = MD->Ele[i].nabr[j] - 1;
				MD->Ele[i].surfH[j] = (inabr >=0) ? (MD->Ele[inabr].zmax + MD->DummyY[inabr]) : ((MD->Ele[i].BC[j] != 1) ? (MD->Ele[i].zmax + MD->DummyY[i]) : Interpolation(&MD->TSD_EleBC[(MD->Ele[i].BC[j]) - 1], t));
			}
			MD->Ele[i].dhBYdx = -1 * (MD->Ele[i].surfY[2] * (MD->Ele[i].surfH[1] - MD->Ele[i].surfH[0]) + MD->Ele[i].surfY[1] * (MD->Ele[i].surfH[0] - MD->Ele[i].surfH[2]) + MD->Ele[i].surfY[0] * (MD->Ele[i].surfH[2] - MD->Ele[i].surfH[1])) / (MD->Ele[i].surfX[2] * (MD->Ele[i].surfY[1] - MD->Ele[i].surfY[0]) + MD->Ele[i].surfX[1] * (MD->Ele[i].surfY[0] - MD->Ele[i].surfY[2]) + MD->Ele[i].surfX[0] * (MD->Ele[i].surfY[2] - MD->Ele[i].surfY[1]));
			MD->Ele[i].dhBYdy = -1 * (MD->Ele[i].surfX[2] * (MD->Ele[i].surfH[1] - MD->Ele[i].surfH[0]) + MD->Ele[i].surfX[1] * (MD->Ele[i].surfH[0] - MD->Ele[i].surfH[2]) + MD->Ele[i].surfX[0] * (MD->Ele[i].surfH[2] - MD->Ele[i].surfH[1])) / (MD->Ele[i].surfY[2] * (MD->Ele[i].surfX[1] - MD->Ele[i].surfX[0]) + MD->Ele[i].surfY[1] * (MD->Ele[i].surfX[0] - MD->Ele[i].surfX[2]) + MD->Ele[i].surfY[0] * (MD->Ele[i].surfX[2] - MD->Ele[i].surfX[1]));
		}
	}
	return DY;
}

int
f_explicit(double t, Model_Data MD)
{
	int             i, j, k, inabr,boundary_index,q,is_outlet,out_i,out_j;
	double        Delta, Gamma;
	double        Rn, G, T, Vel, RH, VP, P, LAI, zero_dh, cnpy_h, rl,
	                r_a, r_s, alpha_r, f_r, eta_s, beta_s, gamma_s, Rmax,
	                Lambda, P_c, qv, qv_sat, ETp;
	double        ThetaRef, ThetaW;
	double        Avg_Y_Surf, Dif_Y_Surf, Grad_Y_Surf, Avg_Sf, Distance,Dif_Y_Surf1,Dif_Y_Surf2;
	double        Cwr, TotalY_Riv, TotalY_Riv_down, CrossA, CrossAdown, AvgCrossA,
	                Perem, Perem_down, Avg_Rough, Avg_Perem, Avg_Y_Riv,
	                Dif_Y_Riv, Grad_Y_Riv, Wid, Wid_down, Avg_Wid;
	double        Avg_Y_Sub, Dif_Y_Sub, Avg_Ksat, Grad_Y_Sub, nabrAqDepth, AquiferDepth,
	                Deficit, elemSatn, satKfunc, effK, effKnabr, TotalY_Ele,
	                TotalY_Ele_down;
	double		Dif_Elev, Grad_Elev;
	double 		c1_square,c2_square,L,b;
	double       *Y, *DY;
	double		time_convert = 60*60*24/UNIT_C;
	double		h_regolith, Tau0, Critical_ShieldsStress,U_star, Renolds, ShieldsStress, q_star,R,q_b,Total_soil_out, Total_soil_in;
	double		Cf = 0.003;
	double		Cos_slope, Sec_slope;
	
	
	
	/*********initialize some parameters******************/
	MD->DY = initiation(MD,MD->DY,t); //All DY[*] = 0
	/********************************************************/

	/*********Calculations among different Processes******************/
	for (i = 0; i < MD->NumEle; i++) 
	{		
		for (j = 0; j < 3; j++) 
		{
            inabr = MD->Ele[i].nabr[j] - 1;
			HydroSubLateral(MD,i,j,inabr,t);
			HydroSurfLateral(MD,i,j,inabr,t);
		}
       
		Infiltration(MD,i);
		ET(MD,i,t);
	}
	Adjust_water(MD);
	//////////////////////calculate new state variables////////////////////
	MD->DY = Calculate_DY(MD,MD->DY);
    
    for (i = 0; i < MD->NumEle; i++)
    {
        MD->DummyY[i] = MD->DummyY[i]+MD->DY[i];
        MD->DummyY[i + MD->NumEle] = MD->DummyY[i + MD->NumEle]+MD->DY[i + MD->NumEle];
        MD->DummyY[i + 2 * MD->NumEle] = MD->DummyY[i + 2 * MD->NumEle]+MD->DY[i + 2 * MD->NumEle];
		MD->DummyY[i + 3 * MD->NumEle] = MD->DummyY[i + 3 * MD->NumEle]+MD->DY[i + 3 * MD->NumEle];
    }
    return 0;
}

double
Interpolation(TSD * Data, double t)
{
	int             i, success;
	double        result;
	i = Data->iCounter;
	success = 0;
	t = t / (UNIT_C);
	while (i < Data->length && t > Data->TS[i][0]) {
		i++;
	}
	if (i == 0) {
		/* t is smaller than the 1st node */
		result = Data->TS[i][1];
	} else if (i >= Data->length) {
		result = Data->TS[i - 1][1];
	} else {
		result = ((Data->TS[i][0] - t) * Data->TS[i - 1][1] + (t - Data->TS[i - 1][0]) * Data->TS[i][1]) / (Data->TS[i][0] - Data->TS[i - 1][0]);
		success = 1;
	}
	if (success == 0) {
		/*
    	printf("\nWarning:  Extrapolation is used ...\n");
    	*/
	}
	return result;
}
