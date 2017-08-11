#include "preincludes.h"
#include <iostream>
#include <math.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "phys_constants.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

/* use this preprocessor symbol to switch between exponential and full Kalman filter */

#define NO_EXPONENTIAL_FILTER

int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);

/* 
 * Applies one iteration of the Kalman filtering algorithm for the current time step
 * (equations A28-A31 in Schmidt et al., MNRAS 440, 2014)
 *
 * @return returns SUCCESS or FAIL
 */
  
int grid::KalmanFiltering()
{
	if (ProcessorNumber != MyProcessorNumber) {
	  return SUCCESS;
	}

	if (NumberOfBaryonFields == 0) {
	  return SUCCESS;
	}
	
	int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;
	if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum)==FAIL)
	    ENZO_FAIL("Error in IdentifyPhysicalQuantities.\n");
	
	int AveVel1Num, AveVel2Num, AveVel3Num;
	
	if ((AveVel1Num = FindField(AveVelocity1, FieldType, NumberOfBaryonFields)) < 0) 
	    ENZO_FAIL("Cannot find x-component of averaged velocity.\n");
	if ((AveVel2Num = FindField(AveVelocity2, FieldType, NumberOfBaryonFields)) < 0) 
	    ENZO_FAIL("Cannot find y-component of averaged velocity.\n");
	if ((AveVel3Num = FindField(AveVelocity3, FieldType, NumberOfBaryonFields)) < 0) 
	    ENZO_FAIL("Cannot find z-component of averaged velocity.\n");

	int VarVel1Num, VarVel2Num, VarVel3Num;
	
	if ((VarVel1Num = FindField(VarVelocity1, FieldType, NumberOfBaryonFields)) < 0) 
	    ENZO_FAIL("Cannot find x-component of velocity variance.\n");
	if ((VarVel2Num = FindField(VarVelocity2, FieldType, NumberOfBaryonFields)) < 0) 
	    ENZO_FAIL("Cannot find y-component of velocity variance.\n");
	if ((VarVel3Num = FindField(VarVelocity3, FieldType, NumberOfBaryonFields)) < 0) 
	    ENZO_FAIL("Cannot find z-component of velocity variance.\n");

	int size = 1;
	for (int dim = 0; dim < GridRank; dim++)
	    size *= GridDimension[dim];

	float DensityUnits = 1, LengthUnits = 1, TemperatureUnits = 1,
	    TimeUnits = 1, VelocityUnits = 1;
	if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits, &TimeUnits, &VelocityUnits, Time) == FAIL) {
	    ENZO_FAIL("Error in GetUnits.");
	}

        // convert correlation time in Gyr to code units
        float corrl_time_code = 3.15e16*KalmanFilterCorrlTime/TimeUnits;
	float gamma = -2.0*pi*dtFixed/(sqrt(3.0)*max(corrl_time_code, tiny_number));

#ifdef EXPONENTIAL_FILTER
        
        /* simple exponential filter for testing purpose */

	float coeff = max(1.0 + gamma, 0.0); 
        if (debug && MyProcessorNumber == ROOT_PROCESSOR)
 	    std::cout << "Exponential filter: " << corrl_time_code << ", "<< dtFixed 
	              << ", " << gamma << ", " << coeff << std::endl;

	for (int n = 0; n < size; n++) {
	    BaryonField[AveVel1Num][n] = 
	        coeff*BaryonField[AveVel1Num][n] + (1.0-coeff)*BaryonField[Vel1Num][n];
	    BaryonField[AveVel2Num][n] = 
	        coeff*BaryonField[AveVel2Num][n] + (1.0-coeff)*BaryonField[Vel2Num][n];
	    BaryonField[AveVel3Num][n] = 
	        coeff*BaryonField[AveVel3Num][n] + (1.0-coeff)*BaryonField[Vel3Num][n];
	}
#else

        /* full Kalman filter (Chahuzac et al., Phys. Fluids 22, 2010; 
           see also Schmidt et al., MNRAS 440, 2014) */

        // convert velocity scale in km/s to code units
	float vel_scale_code = 1.0e5*KalmanFilterVelocityScale/VelocityUnits;
	float sigma = -gamma*vel_scale_code;
	float sigma_sqr = sigma*sigma;
        if (debug && MyProcessorNumber == ROOT_PROCESSOR)
 	    std::cout << "Kalman filter: " << corrl_time_code << ", "<< dtFixed 
	              << ", " << vel_scale_code << ", " << sigma << std::endl;

        float errvar, delta, gain, sigma_sqr_delta;

	for (int dim = 0; dim < GridRank; dim++)
	    for (int n = 0; n < size; n++) {
                // error variance prediction
	        errvar = BaryonField[VarVel1Num+dim][n] + sigma_sqr;
                // fluctuation of random component
                delta = BaryonField[Vel1Num+dim][n] - BaryonField[AveVel1Num+dim][n];
                // expression assumes that initial velocity fluctuation is non-zero
                // (applies to cosmological simulations with initial velocity increments)
                sigma_sqr_delta = fabs(delta)*vel_scale_code;
                // Kalman gain
                gain = errvar/(errvar + sigma_sqr_delta);
                // error variance correction
		BaryonField[VarVel1Num+dim][n] = (1.0 - gain)*errvar;
                // update smoothed component
                BaryonField[AveVel1Num+dim][n] += gain*delta;
            }
#endif

	return SUCCESS;
}

