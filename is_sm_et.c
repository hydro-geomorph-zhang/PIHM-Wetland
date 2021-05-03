/********************************************************************************
 * File        : is_sm_et.c                                                     *
 * Function    : for calculation of canopy ET, IS and snow melt       	        *
 * Version     : Nov, 2007 (2.0)                                                *
 * Developer of PIHM2.0:        Mukesh Kumar (muk139@psu.edu)                   *
 * Developer of PIHM1.0:        Yizhong Qu   (quyizhong@gmail.com)              *
 *------------------------------------------------------------------------------*
 *                                                                              *
 *..............MODIFICATIONS/ADDITIONS in PIHM 2.0.............................*
 * a) Modification of ET components from canopy, ground and transpiration       *
 * b) Addition of snow melt and throughflow processes                           *
 * c) Fully Coupled Inclusion of ET from transpiration and OVL 		        *
 * d) Rain/Snow Fraction calculation according to USACE (1956)		        *
 * e) Change in maximum interception storage due to snow accretion has been     *
 *	accounted for								*
 * f) Incorporation of Interception storage for rainfall as well as snow	*
 ********************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "nvector_serial.h"
#include "sundials_types.h"

#include "pihm.h"
#include "lake.h"

#define EPSILON 0.05
#define UNIT_C 1440  /* 60*24 for calculation of yDot in m/min
				 * units while forcing is in m/day. */
#define multF1 1
#define multF2 5//??
#define multF3 1
#define multF4 1.0
#define C_air 1004.0
#define Lv (2.503*pow(10,6))
#define SIGMA (5.67*pow(10,-8)*60)
#define R_dry 287.04
#define R_v 461.5




double Interpolation(TSD * Data, double t);

void
is_sm_et(double t, double stepsize, Model_Data DS)
{
    int i;
    double totEvap;
    double Delta, Gamma, Lambda;
    double Rn, G, T, Vel, RH, VP, P, LAI = 0.000, zero_dh, cnpy_h, rl,
            r_a, qv_sat, qv, ETp;
    double isval = 0, etval = 0;
    double fracSnow, snowRate, MeltRateGrnd, MeltRateCanopy, eltRate,
            MF, Ts = -3.0, Tr = 1.0, To = 0.0, ret = 0;
    double snow_last, is_last, balance;

    Model_Data MD;

    MD = (Model_Data) DS;

    stepsize = stepsize / UNIT_C;


    for (i = 0; i < MD->NumEle; i++)
    {
        snow_last = MD->EleSnowGrnd[i] + MD->EleSnowCanopy[i];
        is_last = MD->EleIS[i];
        /* Note the dependence on physical units */
        MD->ElePrep[i] = (double) Interpolation(&MD->TSD_Prep[MD->Ele[i].prep - 1], t);
        Rn = Interpolation(&MD->TSD_Rn[MD->Ele[i].Rn - 1], t);
        //G = Interpolation(&MD->TSD_G[MD->Ele[i].G - 1], t);
        G = 0.1 * Rn;
        T = Interpolation(&MD->TSD_Temp[MD->Ele[i].temp - 1], t);
        MD->EleTemp[i] = T;
        Vel = Interpolation(&MD->TSD_WindVel[MD->Ele[i].WindVel - 1], t);
        RH = Interpolation(&MD->TSD_Humidity[MD->Ele[i].humidity - 1], t);
        //VP = Interpolation(&MD->TSD_Pressure[MD->Ele[i].pressure - 1], t);
        VP = 611.2 * exp(17.67 * T / (T + 243.5)) * RH;
        P = 101.325 * pow(10, 3) * pow((293 - 0.0065 * MD->Ele[i].zmax) / 293, 5.26);
        qv = 0.622 * VP / P;
        qv_sat = 0.622 * (VP / RH) / P;
        LAI = Interpolation(&MD->TSD_LAI[MD->Ele[i].LC - 1], t);
        MF = multF2 * Interpolation(&MD->TSD_MeltF[MD->Ele[i].meltF - 1], t);
        /******************************************************************************************/
        /* Snow Accumulation/Melt Calculation				  */
        /******************************************************************************************/
        fracSnow = T < Ts ? 1.0 : T > Tr ? 0 : (Tr - T) / (Tr - Ts);
        snowRate = fracSnow * MD->ElePrep[i];
        /*
         * EleSnowGrnd, EleSnowCanopy, EleISsnowmax,
         * MeltRateGrnd,MeltRateCanopy are the average value prorated
         * over the whole elemental area
         */
        MD->EleSnowGrnd[i] = MD->EleSnowGrnd[i] + (1 - MD->Ele[i].VegFrac) * snowRate * stepsize;
        MD->EleSnowCanopy[i] = MD->EleSnowCanopy[i] + MD->Ele[i].VegFrac * snowRate * stepsize;
        MD->EleISsnowmax[i] = MD->EleSnowCanopy[i] > 0 ? 0.003 * LAI * MD->Ele[i].VegFrac : 0;
        MD->EleISsnowmax[i] = multF1 * MD->EleISsnowmax[i];
        if (MD->EleSnowCanopy[i] > MD->EleISsnowmax[i])
        {
            MD->EleSnowGrnd[i] = MD->EleSnowGrnd[i] + MD->EleSnowCanopy[i] - MD->EleISsnowmax[i];
            MD->EleSnowCanopy[i] = MD->EleISsnowmax[i]; /* update in 2.2  */
        }
        MeltRateGrnd = MeltRateCanopy = (T > To ? (T - To) * MF : 0); /* Note the units for
											 * MF. */
        if (MD->EleSnowGrnd[i] > MeltRateGrnd * stepsize)
        {
            MD->EleSnowGrnd[i] = MD->EleSnowGrnd[i] - MeltRateGrnd * stepsize;
        }
        else
        {
            MeltRateGrnd = MD->EleSnowGrnd[i] / stepsize;
            MD->EleSnowGrnd[i] = 0;
        }
        if (MD->EleSnowCanopy[i] > MeltRateCanopy * stepsize)
        {
            MD->EleSnowCanopy[i] = MD->EleSnowCanopy[i] - MeltRateCanopy * stepsize;
        }
        else
        {
            MeltRateCanopy = MD->EleSnowCanopy[i] / stepsize;
            MD->EleSnowCanopy[i] = 0;
        }
        /************************************************************************/
        /* ThroughFall and Evaporation from canopy			 */
        /************************************************************************/
        /*
         * EleIS, EleET[0] and ret are prorated for the whole
         * element. Logistics are simpler if assumed in volumetric
         * form by multiplication of Area on either side of equation
         */
        MD->EleISmax[i] = multF1 * MD->ISFactor[MD->Ele[i].LC - 1] * LAI * MD->Ele[i].VegFrac;
        /* Note the dependence on physical units */
        if (LAI > 0.0)
        {

            /*
             * zero_dh=Interpolation(&MD->TSD_DH[MD->Ele[i].LC-1],
             *  t); cnpy_h =
             * zero_dh/(1.1*(0.0000001+log(1+pow(0.007*LAI,0.25)))
             * ); if(LAI<2.85)	{ rl= 0.0002 +
             * 0.3*cnpy_h*pow(0.07*LAI,0.5); } else { rl=
             * 0.3*cnpy_h*(1-(zero_dh/cnpy_h)); }
             */ rl = Interpolation(&MD->TSD_RL[MD->Ele[i].LC - 1], t);
            //r_a = log(MD->Ele[i].windH / rl) * log(10 * MD->Ele[i].windH / rl) / (Vel * 0.16);
            r_a = 12 * 4.72 * log(MD->Ele[i].windH / rl) / (0.54 * Vel / UNIT_C / 60 + 1) / UNIT_C / 60;

            Gamma = 4 * 0.7 * SIGMA * UNIT_C * R_dry / C_air * pow(T + 273.15, 4) / (P / r_a) + 1;
            Delta = Lv * Lv * 0.622 / R_v / C_air / pow(T + 273.15, 2) * qv_sat;

            ETp = (Rn * Delta + Gamma * (1.2 * Lv * (qv_sat - qv) / r_a)) / (1000.0 * Lv * (Delta + Gamma));

            MD->EleET[i][0] = MD->pcCal.Et0 * MD->Ele[i].VegFrac * (pow((MD->EleIS[i] < 0 ? 0 : (MD->EleIS[i] > MD->EleISmax[i] ? MD->EleISmax[i] : MD->EleIS[i])) / MD->EleISmax[i], 1.0 / 2.0)) * ETp;
            MD->EleET[i][0] = MD->EleET[i][0] < 0 ? 0 : MD->EleET[i][0];

            //MD->EleET[i][0] = MD->pcCal.Et0 * MD->Ele[i].VegFrac * (LAI / MD->Ele[i].LAImax) * (pow((MD->EleIS[i] < 0 ? 0 : MD->EleIS[i]) / MD->EleISmax[i], 2.0 / 3.0)) * (Rn * (1 - MD->Ele[i].Albedo) * Delta + (1.2 * 1003.5 * ((VP / RH) - VP) / r_a)) / (1000 * 2441000.0 * (Delta + Gamma));
            MD->EleTF[i] = MD->EleIS[i] <= 0 ? 0 : 5.65 * pow(10, -2) * MD->EleISmax[i] * exp(3.89 * (MD->EleIS[i] < 0 ? 0 : MD->EleIS[i]) / MD->EleISmax[i]); /* Note the dependece on
																							 * physical units */
            MD->EleTF[i] = multF3 * MD->EleTF[i];
        }
        else
        {
            MD->EleET[i][0] = 0.0;
            MD->EleTF[i] = 0.0;
        }
        if (MD->EleTF[i] < 0)MD->EleTF[i] = 0.0;
        if (MD->EleTF[i] * stepsize > MD->EleIS[i])MD->EleTF[i] = MD->EleIS[i] / stepsize;
        if (MD->EleIS[i] >= MD->EleISmax[i])
        {
            if (((1 - fracSnow) * MD->ElePrep[i] * MD->Ele[i].VegFrac + MeltRateCanopy) >= MD->EleET[i][0] + MD->EleTF[i])
            {
                MD->EleETloss[i] = MD->EleET[i][0];
                ret = MD->EleTF[i] + (MD->EleIS[i] - MD->EleISmax[i]) / stepsize + (((1 - fracSnow) * MD->ElePrep[i] * MD->Ele[i].VegFrac + MeltRateCanopy) - (MD->EleET[i][0] + MD->EleTF[i]));
                isval = MD->EleISmax[i];
                //MD->EleIS[i] = MD->EleISmax[i];
            }
            else if ((((1 - fracSnow) * MD->ElePrep[i] * MD->Ele[i].VegFrac + MeltRateCanopy) < MD->EleET[i][0] + MD->EleTF[i]) && (MD->EleIS[i] + stepsize * ((1 - fracSnow) * MD->ElePrep[i] * MD->Ele[i].VegFrac + MeltRateCanopy - MD->EleET[i][0] - MD->EleTF[i]) <= 0))
            {
                MD->EleETloss[i] = (MD->EleET[i][0] / (MD->EleET[i][0] + MD->EleTF[i])) * (MD->EleIS[i] / stepsize + ((1 - fracSnow) * MD->ElePrep[i] * MD->Ele[i].VegFrac + MeltRateCanopy));
                ret = (MD->EleTF[i] / (MD->EleET[i][0] + MD->EleTF[i])) * (MD->EleIS[i] / stepsize + ((1 - fracSnow) * MD->ElePrep[i] * MD->Ele[i].VegFrac + MeltRateCanopy));
                MD->EleET[i][0] = MD->EleETloss[i];
                //MD->EleIS[i] = 0;
                isval = 0;
                MD->EleETloss[i] = MD->EleET[i][0];
            }
            else
            {
                if (MD->EleIS[i] + stepsize * ((1 - fracSnow) * MD->ElePrep[i] * MD->Ele[i].VegFrac + MeltRateCanopy - MD->EleET[i][0] - MD->EleTF[i]) >= MD->EleISmax[i])
                {
                    isval = MD->EleISmax[i];
                    ret = MD->EleTF[i] + (MD->EleIS[i] - MD->EleISmax[i]) / stepsize + (((1 - fracSnow) * MD->ElePrep[i] * MD->Ele[i].VegFrac + MeltRateCanopy) - (MD->EleET[i][0] + MD->EleTF[i]));
                    MD->EleETloss[i] = MD->EleET[i][0];
                }
                else
                {
                    isval = MD->EleIS[i] + stepsize * (((1 - fracSnow) * MD->ElePrep[i] * MD->Ele[i].VegFrac + MeltRateCanopy) - MD->EleET[i][0] - MD->EleTF[i]);
                    //MD->EleIS[i] = MD->EleIS[i] + stepsize * (((1 - fracSnow) * MD->ElePrep[i] * MD->Ele[i].VegFrac + MeltRateCanopy) - MD->EleET[i][0] - MD->EleTF[i]);
                    ret = MD->EleTF[i];
                    MD->EleETloss[i] = MD->EleET[i][0];
                }
            }
        }
        else if ((MD->EleIS[i] < MD->EleISmax[i]) && ((MD->EleIS[i] + (((1 - fracSnow) * MD->ElePrep[i] * MD->Ele[i].VegFrac + MeltRateCanopy) - MD->EleET[i][0] - MD->EleTF[i]) * stepsize) >= MD->EleISmax[i]))
        {
            MD->EleETloss[i] = MD->EleET[i][0];
            isval = MD->EleISmax[i];
            ret = MD->EleTF[i] + (((MD->EleIS[i] + (((1 - fracSnow) * MD->ElePrep[i] * MD->Ele[i].VegFrac + MeltRateCanopy) - MD->EleET[i][0] - MD->EleTF[i]) * stepsize) - MD->EleISmax[i])) / stepsize;
        }
        else if ((MD->EleIS[i] < MD->EleISmax[i]) && ((MD->EleIS[i] + (((1 - fracSnow) * MD->ElePrep[i] * MD->Ele[i].VegFrac + MeltRateCanopy) - MD->EleET[i][0] - MD->EleTF[i]) * stepsize) <= 0))
        {
            if ((MD->EleET[i][0] > 0) || (MD->EleTF[i] > 0))
            {
                MD->EleETloss[i] = (MD->EleET[i][0] / (MD->EleET[i][0] + MD->EleTF[i])) * (MD->EleIS[i] / stepsize + ((1 - fracSnow) * MD->ElePrep[i] * MD->Ele[i].VegFrac + MeltRateCanopy));
                ret = (MD->EleTF[i] / (MD->EleET[i][0] + MD->EleTF[i])) * (MD->EleIS[i] / stepsize + ((1 - fracSnow) * MD->ElePrep[i] * MD->Ele[i].VegFrac + MeltRateCanopy));
                MD->EleET[i][0] = MD->EleETloss[i];
            }
            else
            {
                MD->EleET[i][0] = 0;
                ret = 0;
            }
            MD->EleETloss[i] = MD->EleET[i][0];
            isval = 0;
        }
        else
        {
            isval = MD->EleIS[i] + (((1 - fracSnow) * MD->ElePrep[i] * MD->Ele[i].VegFrac + MeltRateCanopy) - MD->EleET[i][0] - MD->EleTF[i]) * stepsize;
            MD->EleETloss[i] = MD->EleET[i][0];
            ret = MD->EleTF[i];
        }
        MD->EleNetPrep[i] = (1 - MD->Ele[i].VegFrac) * (1 - fracSnow) * MD->ElePrep[i] + ret + MeltRateGrnd;
        MD->EleTF[i] = ret;
        MD->EleIS[i] = isval;
    }
}
