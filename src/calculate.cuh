#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <omp.h>
#include <vector>

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
    t_timers *timers, 
    t_message *msg);
/*----------------------------*/

// --- Functions for preCICE
////////////////////////////////////////////////////////////////
EXPORT_DLL int getCellsNumbersAtBC(
    t_mesh *mesh
);
/*----------------------------*/

////////////////////////////////////////////////////////////////
EXPORT_DLL std::vector<double> getCellsCoordinatesAtBC(
    t_mesh *mesh
);

////////////////////////////////////////////////////////////////
EXPORT_DLL std::vector<int> getCellsIDsAtBC(
    t_mesh *mesh
);

////////////////////////////////////////////////////////////////
__global__ void g_get_coupling_values(
    int nTasks,
    t_arrays *arrays,
    int *cellIDs,
    double *hBC,
    double *huBC,
    double *hvBC);

////////////////////////////////////////////////////////////////
__global__ void g_write_coupling_values(
    int nTasks,
    t_arrays *arrays,
    int *cellIDs,
    double *hBC,
    double *huBC,
    double *hvBC);