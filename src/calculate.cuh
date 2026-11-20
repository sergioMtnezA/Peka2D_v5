#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <omp.h>

#include "cuda.h"
#include "cuda_runtime.h"
#include "cublas_v2.h"

#include "define.h"
#include "structs.h"
#include "loadData.h"
#include "writeInOut.cuh"
#include "memoryManage.cuh"
#include "water.cuh"
#include "boundary.cuh"
#include "cuTilities.cuh"

#if SET_SOLUTE
#include "solute.cuh"
#endif


////////////////////////////////////////////////////////////////
EXPORT_DLL int computeSimulation(
    t_parameters spar,
    t_mesh *mesh, 
    t_arrays *carrays,
    Peka2D_Rain *rainGroup, // Rainfall time series
    int nRainPolygons,      // Number of rainfall polygons
    t_timers *timers,
	t_message *msg);
/*----------------------------*/

////////////////////////////////////////////////////////////////
EXPORT_DLL int computeInitialBoundaryConditions(
    t_arrays *carrays,
    t_arrays *garrays,     
    t_cuPtr *cuPtr,
    t_message *msg);
/*----------------------------*/

////////////////////////////////////////////////////////////////
EXPORT_DLL int computeInitialMassBalance(
    t_arrays *carrays,
    t_arrays *garrays,     
    t_cuPtr *cuPtr,
    t_message *msg);
/*----------------------------*/

////////////////////////////////////////////////////////////////
EXPORT_DLL void generateTimeStep(
    double *t,
    t_arrays *carrays,
    t_arrays *garrays,     
    t_cuPtr *cuPtr,
    Peka2D_Rain *rainGroup, // Rainfall time series
    int nRainPolygons,      // Number of rainfall polygons
    double *hostRainIntensity, // Rainfall intensities on the CPU [m/s]
    t_timers *timers, 
    t_message *msg);
/*----------------------------*/

