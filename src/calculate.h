#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <omp.h>

#include "define.h"
#include "structs.h"
#include "loadData.h"
#include "writeInOut.h"
#include "water.h"
#include "utilities.h"

////////////////////////////////////////////////////////////////
/**
 * @brief This function computes the time integration of the equations.
 */
EXPORT_DLL int computeSimulation(
    t_parameters spar,
    t_mesh *mesh, 
    t_arrays *carrays,
    t_timers *timers,
	t_message *msg);
/*----------------------------*/

////////////////////////////////////////////////////////////////
/**
 * @brief This function computes the initial mass balance.
 */
EXPORT_DLL int computeInitialMassBalance(
    t_arrays *carrays,
    t_message *msg);
/*----------------------------*/

////////////////////////////////////////////////////////////////
/**
 * @brief This function updates the variables in the next the time step.
 */
EXPORT_DLL void generateTimeStep(
    double *t,
    t_arrays *carrays,
    t_timers *timers, 
    t_message *msg);
/*----------------------------*/

