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
#include "cuTilities.cuh"

////////////////////////////////////////////////////////////////
EXPORT_DLL int createArraysCudaMemory(
    t_arrays *carrays,
    t_arrays *garrays,
    t_cuPtr *cuPtr,
    t_message *msg);
/*----------------------------*/

////////////////////////////////////////////////////////////////
EXPORT_DLL int allocateArraysCudaMem(
    int NCwall, int ncells, int nWallCell, int nwc, int nwb,   
    t_cuPtr *cuPtr);
/*----------------------------*/

////////////////////////////////////////////////////////////////
EXPORT_DLL int copyControlArraysCudaMem(
    t_arrays *carrays,
    t_arrays *garrays,
    t_cuPtr *cuPtr);
/*----------------------------*/  

////////////////////////////////////////////////////////////////
int copyMeshArraysCudaMem(
    t_arrays *carrays,
    t_arrays *garrays,
    t_cuPtr *cuPtr);
/*----------------------------*/


////////////////////////////////////////////////////////////////
__global__ void assignGArraysToCudaMem(t_arrays *garrays,
	//------------------------cells
	int *activeC,
	int *actCells,
	int *cidx,
	int *nneig,
	double *z,
	double *h,
	double *hu,
	double *hv,
	double *u,
	double *v,
	double *modulou,
	double *sqrh,
	double *area,
	double *nman,
	double *SOX,
	double *SOY,
	double *mass,
	//--------------------- cells*NCwall
	double *dh,
	double *dhu,
	double *dhv,
	int *solidWallByCell,
	int *neighCell,
	int *neighWall,
	double *normalXbyCell,
	double *normalYbyCell,
	//---------------------- internal walls
	int *activeW,
	int *actWalls,
	int *widx,
	int *idx1,
	int *idx2,
	int *idw1,
	int *idw2,
	double *normalX,
	double *normalY,
	double *deltaX,
	double *length,
	double *distNormal,
	double *distCentX,
	double *distCentY,
	double *nman2wall,
	double *gp,
	int *typeOfBound,
	int *solidWall,
	double *qnormalL,
	double *localDt);
/*----------------------------*/

////////////////////////////////////////////////////////////////
EXPORT_DLL int freeCudaMemory(t_cuPtr *cuPtr);
/*----------------------------*/

