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
#if SET_COUPLING
////////////////////////////////////////////////////////////////
__global__ void readWriteCouplingValues(
    int nTasks,
    char cond,
    t_arrays *arrays);
/*----------------------------*/

////////////////////////////////////////////////////////////////
EXPORT_DLL std::vector<double> getCellsCoordinatesAtBC(
    t_mesh *mesh,
    int k);
/*----------------------------*/

////////////////////////////////////////////////////////////////
// --- update the arrays in CPU <-> GPU
EXPORT_DLL void updateCouplingArrays(
    t_arrays *carrays,   
    t_cuPtr *cuPtr,
    char cond);
/*----------------------------*/

////////////////////////////////////////////////////////////////
// --- create a bounding box for preCICE mapping
EXPORT_DLL std::vector<double> CouplingBbox(
    t_mesh *mesh,
    int k,
    double margin);
/*----------------------------*/

////////////////////////////////////////////////////////////////
// --- Function to map the nearest PEKA coupling cell for each OF face in x-y plane
EXPORT_DLL void mapOFwithPEKA(
    t_mesh *mesh,
    int k,
    t_coupInterface &ci);
/*----------------------------*/

#endif
