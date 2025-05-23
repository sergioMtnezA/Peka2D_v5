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
EXPORT_DLL int copyComputationControls(
    t_arrays *carrays,
    t_arrays *garrays,
    t_cuPtr *cuPtr);
/*----------------------------*/  

////////////////////////////////////////////////////////////////
EXPORT_DLL int allocateArraysCudaMem(
    int NCwall, int ncells, int nWallCell, int nwc, int nwb,   
    t_cuPtr *cuPtr);
/*----------------------------*/

////////////////////////////////////////////////////////////////
int copyMeshArraysCudaMem(
    t_arrays *carrays,
    t_arrays *garrays,
    t_cuPtr *cuPtr);
/*----------------------------*/

////////////////////////////////////////////////////////////////
__global__ void assignMeshArraysToCudaMem(t_arrays *garrays,
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
    int *typeWallByCell,
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

////////////////////////////////////////////////////////////////
EXPORT_DLL int allocateBoundArraysCudaMem(
    int nOBC, int nInlet, int nOutlet, 
    int nTotalBoundCells, int nTotalInnerCells,
    int nTotalPointSeries, 
    int nSolutes,
    t_cuPtr *cuPtr);
/*----------------------------*/

////////////////////////////////////////////////////////////////
EXPORT_DLL int copyBoundSetupArraysCudaMem(
    t_arrays *carrays,
    t_arrays *garrays,
    t_cuPtr *cuPtr);
/*----------------------------*/

////////////////////////////////////////////////////////////////
EXPORT_DLL int copyBoundMeshArraysCudaMem(
    t_arrays *carrays,
    t_arrays *garrays,
    t_cuPtr *cuPtr);
/*----------------------------*/

////////////////////////////////////////////////////////////////
__global__ void assignBoundArraysToCudaMem(t_arrays *garrays,
	//------------------------ bound geometry
	int *nCellsOBC,
    int *iniIndexOBC,
    int *idBoundOBC,
    int *typeOBC,
    int *flagInitializeOBC,
    double *blockSectionOBC,
    double *normalXOBC,
    double *normalYOBC,
    double *totalLengthOBC,
    double *totalAreaOBC,
    int *cellZminOBC,
    int *nInnerCellsOBC,
    int *iniInnerIndexOBC,
    //----------------------- bound cells
    int *cidxBound,
    double *zCellBound,
    double *areaCellBound,    
    double *nxWallBound,
    double *nyWallBound,
    double *lWallBound,
    //----------------------- inner cells
    int *cidxInner,
    //----------------------- time series
    int *nPointsSeriesOBC,
    int *iniIndexSeriesOBC,
    double *tSeriesOBC,
    double *qSeriesOBC,
    double *hzSeriesOBC,
    double *frSeriesOBC,
    double *phiSeriesOBC,
    //----------------------- mass balance pointers
    double *qBoundByCell,
    double *mBoundByCell,
    double *mInnerByCell,
    double *qInByInlet,
    double *qOutByOutlet,
    double *mInByInlet, 
    double *mOutByOutlet);
/*----------------------------*/

////////////////////////////////////////////////////////////////
EXPORT_DLL int freeBoundaCudaMemory(
    int nOBC, int nInlet, int nOutlet, 
    int nTotalBoundCells, int nTotalInnerCells,
    int nTotalPointSeries,
    t_cuPtr *cuPtr);
/*----------------------------*/



#if SET_SOLUTE
////////////////////////////////////////////////////////////////
EXPORT_DLL int allocateSoluteArraysCudaMem(
    int nSolutes, 
    int NCwall, int ncells, int nWallCell, int nwc, int nwb, 
    t_cuPtr *cuPtr);
/*----------------------------*/


////////////////////////////////////////////////////////////////
int copySoluteArraysCudaMem(
    int nSolutes,
    t_arrays *carrays,
    t_arrays *garrays,
    t_cuPtr *cuPtr);
/*----------------------------*/



////////////////////////////////////////////////////////////////
__global__ void assignSoluteArraysToCudaMem(int nSolutes, t_arrays *garrays,
	//------------------------solutes
	int *typeDiff,
	double *k_xx,
	double *k_yy,
	//------------------------solutes*cells
	double *hphi,
	double *phi,
    double *localDtd,
	double *BTcell,
	//------------------------solutes*cells*NCwalls
	double *dhphi,
    double *Bwall);
/*----------------------------*/

EXPORT_DLL int freeSoluteCudaMemory(
	int nSolutes, 
	t_cuPtr *cuPtr);
/*----------------------------*/
#endif

