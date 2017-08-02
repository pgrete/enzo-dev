// LI implemented in this Enzo version the modifications performed originally by WS, Oct. 2011

#include "preincludes.h"
#include <iostream>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

#define MIN_EXPN -23.02585092994045684
#define TWO_PI 6.283185307179586

int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);

/* 
 * Increments the average momentum for the shear-improved technique
 *
 * @return returns SUCCESS or FAIL
 */
  
int grid::KalmanFiltering()
{
	if (ProcessorNumber != MyProcessorNumber)
	{
	    return SUCCESS;
	}
	
 	std::cout << "[" << MyProcessorNumber << "] Incrementing average momentum" << std::endl;

	int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;
	if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum)==FAIL)
	{
	    fprintf(stderr, "Error in IdentifyPhysicalQuantities.\n");
	    return FAIL;
	}
	
	int AveMomt1Num, AveMomt2Num, AveMomt3Num;
	
	if ((AveMomt1Num = FindField(AveMomtDensity1, FieldType, NumberOfBaryonFields)) < 0) 
	{
	    fprintf(stderr, "Cannot find average x-component of averaged momentum.\n");
	    return FAIL;
	}
	if ((AveMomt2Num = FindField(AveMomtDensity2, FieldType, NumberOfBaryonFields)) < 0) 
	{
	    fprintf(stderr, "Cannot find average y-component of averaged momentum.\n");
	    return FAIL;
	}
	if ((AveMomt3Num = FindField(AveMomtDensity3, FieldType, NumberOfBaryonFields)) < 0) 
	{
	    fprintf(stderr, "Cannot find average z-component of averaged momentum.\n");
	    return FAIL;
	}

	int size = 1;
	for (int dim = 0; dim < GridRank; dim++)
	    size *= GridDimension[dim];

	float DensityUnits = 1, LengthUnits = 1, TemperatureUnits = 1,
	    TimeUnits = 1, VelocityUnits = 1;
	if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits, &TimeUnits, &VelocityUnits, Time) == FAIL) {
	    ENZO_FAIL("Error in GetUnits.");
	}

        // convert KalmanFilterCorrlTime in Gyr to code units
        float corrl_time_code = 3.15e16*KalmanFilterCorrlTime/TimeUnits;
	float gamma = -TWO_PI*dtFixed/max(corrl_time_code, tiny_number);
	float coeff = max(1.0 + gamma, 0.0); 
        if (MyProcessorNumber == ROOT_PROCESSOR)
 	    std::cout << "KalmanFiltering: " << corrl_time_code << ", "<< dtFixed 
	              << ", " << gamma << ", " << coeff << std::endl; 

	for (int n = 0; n < size; n++) {
	    BaryonField[AveMomt1Num][n] = coeff*BaryonField[AveMomt1Num][n] + 
	                                  (1.0-coeff)*BaryonField[DensNum][n]*BaryonField[Vel1Num][n];
	    BaryonField[AveMomt2Num][n] = coeff*BaryonField[AveMomt2Num][n] + 
	                                  (1.0-coeff)*BaryonField[DensNum][n]*BaryonField[Vel2Num][n];
	    BaryonField[AveMomt3Num][n] = coeff*BaryonField[AveMomt3Num][n] + 
	                                  (1.0-coeff)*BaryonField[DensNum][n]*BaryonField[Vel3Num][n];
	}

	return SUCCESS;
}

