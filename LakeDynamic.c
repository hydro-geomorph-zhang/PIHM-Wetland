//
//  Lake_land_exchange.c
//  PIHM-Lake
//
//  Created by Yu Zhang on 2/21/16.
//  Copyright © 2016 Yu Zhang. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>



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


double 		avgY(double diff, double yi, double yinabr);
double        Interpolation(TSD * Data, double t);


void LakeDynamic(double t, LakeData LD, Model_Data MD)
{
	int             i, j, k, boundary_index,is_outlet,out_i,out_j;
    int             EleID;
	double        Delta, Gamma;
	double        Rn, G, T, Vel, RH, VP, P, LAI, zero_dh, cnpy_h, rl,
	                r_a, r_s, alpha_r, f_r, eta_s, beta_s, gamma_s, Rmax,
	                Lambda, P_c, qv, qv_sat, ETp;
	double        ThetaRef, ThetaW;
	double        Avg_Y_Surf, Dif_Y_Surf, Grad_Y_Surf, Avg_Sf, Distance,Dif_Y_Surf1,Dif_Y_Surf2;
	double        Cwr, TotalY_Riv, TotalY_Riv_down, CrossA, CrossAdown,
	                Perem, Perem_down, Avg_Rough, Avg_Perem, Avg_Y_Riv,
	                Dif_Y_Riv, Grad_Y_Riv, Wid, Wid_down, Avg_Wid;
	double        Avg_Y_Sub, Dif_Y_Sub, Avg_Ksat, Grad_Y_Sub, AquiferDepth,
	                Deficit, elemSatn, satKfunc, effK, effKnabr, TotalY_Ele,
	                TotalY_Ele_down;
	double		Dif_Elev, Grad_Elev;
	double 		c1_square,c2_square,L,b;
	double		h_regolith, Tau0, Critical_ShieldsStress,U_star, Renolds, ShieldsStress, q_star,R,q_b,Total_soil_out, Total_soil_in;
	double		Cf = 0.003;
    double        D_surf=0,D_GW=0;
	

	
    //*****/**Calculate Lake Water**/**************/
	
	for(i=0;i<LD->NumLake;i++)
	{    
        for(j=0;j<LD->Lake[i].NumBankEle;j++)
		{
            EleID = LD->Lake[i].BankEle[j]-1;
            
            for (k=0;k<3;k++)
            {
                if (MD->Ele[EleID].BC[k] == 4||MD->Ele[EleID].BC[k] == 5||MD->Ele[EleID].BC[k] == 6)
                {
                    
                    LD->FluxSurf[i][j] = MD->FluxSurf[EleID][k];
                    D_surf = D_surf + LD->FluxSurf[i][j]/LD->Lake[i].SurfArea;
                    
                    
                    LD->FluxSub[i][j] = MD->FluxSub[EleID][k];
                    D_GW = D_GW + LD->FluxSub[i][j]/LD->Lake[i].SurfArea;
                }
            }
		}
		
        /*****Calculate water exchange between lake surface water and groundwater*****/
        Grad_Y_Sub = (LD->LakeSurfDepth[i] + LD->Lake[i].BedElev - (LD->LakeGW[i] + LD->Lake[i].BaseElev)) / (LD->Lake[i].BedElev-LD->Lake[i].BaseElev);
        Grad_Y_Sub = ((LD->LakeSurfDepth[i] < EPS / 100) && (Grad_Y_Sub > 0)) ? 0 : Grad_Y_Sub;
        
        LD->Infil[i] = LD->Lake[i].Kv * Grad_Y_Sub;
        D_surf = D_surf - LD->Infil[i];
        D_GW = D_GW + LD->Infil[i];
        
		/*****Calculate Evaporation*****/
        LD->Lake[i].Precip = Interpolation(&(MD->TSD_Prep[LD->Lake[i].Meteo - 1]), t);
        Rn = Interpolation(&(MD->TSD_Rn[LD->Lake[i].Meteo - 1]), t);
        //G = Interpolation(&MD->TSD_G[MD->Ele[i].G - 1], t);
        G = 0.1 * Rn;
        LD->Lake[i].Temp = Interpolation(&(MD->TSD_Temp[LD->Lake[i].Meteo - 1]), t);
        Vel = Interpolation(&(MD->TSD_WindVel[LD->Lake[i].Meteo - 1]), t);
        RH = Interpolation(&(MD->TSD_Humidity[LD->Lake[i].Meteo - 1]), t);
        VP = 611.2 * exp(17.67 * LD->Lake[i].Temp / (LD->Lake[i].Temp + 243.5)) * RH;
        P = 101.325 * pow(10, 3) * pow((293 - 0.0065 * MD->Ele[i].zmax) / 293, 5.26);
        qv = 0.622 * VP / P;
        qv_sat = 0.622 * (VP / RH) / P;
        //P = 101.325 * pow(10, 3) * pow((293 - 0.0065 * MD->Ele[i].zmax) / 293, 5.26);
        //Delta = 2503 * pow(10, 3) * exp(17.27 * T / (T + 237.3)) / (pow(237.3 + T, 2));
        //Gamma = P * 1.0035 * 0.92 / (0.622 * 2441);
        /*
            * zero_dh=Interpolation(&MD->TSD_DH[MD->Ele[i].LC-1], t);
            * cnpy_h =
            * zero_dh/(1.1*(0.0000001+log(1+pow(0.007*LAI,0.25))));
            * if(LAI<2.85)	{ rl= 0.0002 + 0.3*cnpy_h*pow(0.07*LAI,0.5);
            * } else { rl= 0.3*cnpy_h*(1-(zero_dh/cnpy_h)); }
        */
        rl = Interpolation(&(MD->TSD_RL[LD->Lake[i].Meteo - 1]), t);
        r_a = 12 * 4.72 * log(MD->Ele[i].windH / rl) / (0.54 * Vel / UNIT_C / 60 + 1) / UNIT_C / 60;

        Gamma = 4 * 0.7 * SIGMA * UNIT_C * R_dry / C_air * pow(LD->Lake[i].Temp + 273.15, 4) / (P / r_a) + 1;
        Delta = Lv * Lv * 0.622 / R_v / C_air / pow(LD->Lake[i].Temp + 273.15, 2) * qv_sat;
        ETp = (Rn * Delta + Gamma * (1.2 * Lv * (qv_sat - qv) / r_a)) / (1000.0 * Lv * (Delta + Gamma));
        LD->ET[i] = ETp;
        D_surf = D_surf - LD->ET[i];
        
        /*****Calculate New State Variables*****/
        LD->LakeSurfDepth[i] = LD->LakeSurfDepth[i] + D_surf/UNIT_C + LD->Lake[i].Precip/UNIT_C;
        LD->LakeGW[i] = LD->LakeGW[i] + D_GW/(LD->Lake[i].ThetaS * UNIT_C);
		D_surf = 0;
		D_GW = 0;
	}
}

























