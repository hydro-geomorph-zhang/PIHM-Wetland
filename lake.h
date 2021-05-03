

typedef struct lake_type_structure {
	double		SurfArea; /*surface area of lake*/
	double		SurfIC; /*Surface elevation when it is full of water*/
	double		BedElev; /*Lake Bed elevation*/
	double		BaseElev; /*Bedrock elevation*/
	int 		Soil; /*soil type beneath lake*/
	int 		Meteo; /*The Meteorological type of the lake*/
	int 		Macropore; /*1 or 0. 1 means macropore is considered*/
	double 		Kh; /*horizontal hydraulic conductivity*/
	double 		Kv; /*vertical hydraulic conductivity*/
	double 		ThetaS; /*Porosity*/
	double 		ThetaR; /*Residual porosity*/
	double 		macKh; /*Horizontal macropore hydraulic conductivity*/
	int			NumStreamEle; /*Number of stream element at boundary*/
	int			NumBankEle; /*Number of bank element at boundary*/
	int			NumEle; /*Number of total element at boundary*/
	int 		*StreamEle;
	int 		*BankEle;
    double      Precip; /* precipitation*/
	double      Temp; /* precipitation*/
	int			*StreamIndex;
	int			NumLakeNode;
    double      BankElevAdjust;
    double      OutElevDecline;
	double		Weir_length;
	double		Orifice_height;
} 	Lake;

typedef struct Lake_soil_structure {
	int			index;
	double 		Kh; /*horizontal hydraulic conductivity*/
	double 		Kv; /*vertical hydraulic conductivity*/
	double 		ThetaS; /*Porosity*/
	double 		ThetaR; /*Residual porosity*/
	double 		macKh; /*Horizontal macropore hydraulic conductivity*/
} 	Lake_soil;

typedef struct lake_data_structure {
	Lake 		*Lake; /*store lake information*/
	Lake_soil 	*LakeSoil;
	int			NumLake; /*Number of lakes*/
	int			NumLakeSoil; /*Number of lakes*/
    int			NumRiverInterest; /*Number of RiverInterest*/
    
	double      **BankSurf;
	double      **BankGW;
	double      *LakeSurfDepth;
	double      *LakeGW;
	double      *LakeDummyY;
	double      **FluxSurf; /*store lake information*/
	double      **FluxStream; /*store lake information*/
	double      **FluxSub; /*store lake information*/
	double      *ET; /*store lake information*/
	double      *Infil; /*store lake information*/
    int         **LakeBoundary;  //LakeBoundary[NumLake][NumStreamEle+NumBankEle]
    double      **LakeBoundaryNodes; //LakeBoundary[(NumStreamEle+NumBankEle)*2][2]
    int         *RiverInterest;
    int      ***RiverInterestInformation;
    double    *PrintVar[10];
} 	*LakeData;


