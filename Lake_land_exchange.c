//
//  Lake_land_exchange.c
//  PIHM-Lake
//
//  Created by Yu Zhang on 2/21/16.
//  Copyright © 2016 Yu Zhang. All rights reserved.
//

#include <stdio.h>

#include "pihm.h"
#include "lake.h"

#include "sundials_types.h"
#include "nvector_serial.h"

void Lake_land_exchange(Model_Data DS, LakeData LD)
{
    int k,w;
    
    for(k=0;k<LD->NumLake;k++)
    {    
        for(w=0;w<LD->Lake[k].NumBankEle;w++)
        {
            LD->BankSurf[k][w] = DS->DummyY[LD->Lake[k].BankEle[w]-1];
            LD->BankGW[k][w] = DS->DummyY[LD->Lake[k].BankEle[w]-1+2*DS->NumEle];
        }
        
        DS->LakeSurfElev[k] = LD->LakeSurfDepth[k]+LD->Lake[k].BedElev;
        DS->LakeGWElev[k] = LD->LakeGW[k] + LD->Lake[k].BaseElev;
    }
    
}
