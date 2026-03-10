#include "define.h"
#include "structs.h"
#include "loadData.h"
#include "cuTilities.cuh"

/**
 * @brief This function updates the open boundaries.
 * @param arrays This pointer variable passes the arrays structure.
 */
__device__ void d_initialize_inlet(t_arrays *arrays, 
	int i, int idb,
	double *qBoundByCell, double *mBoundByCell, double *mInnerByCell,
	double *qInByInlet, double *mInByInlet);


/**
 * @brief This function updates the open boundaries.
 * @param arrays This pointer variable passes the arrays structure.
 */
__global__ void g_update_open_boundary(int nTasks, t_arrays *arrays, 
	double *qBoundByCell, double *mBoundByCell, double *mInnerByCell,
	double *qInByInlet, double *mInByInlet,
	double *qOutByOutlet, double *mOutByOutlet);		   


/**
 * @brief This function updates the open boundaries.
 * @param arrays This pointer variable passes the arrays structure.
 */
__device__ int d_get_index(double valor, int n, int idx0, double *x);


/**
 * @brief This function updates the open boundaries.
 * @param arrays This pointer variable passes the arrays structure.
 */
__device__ double d_interpolate_vector(double valor, 
	int n, int idx0,
	int tidx, 
	double *x, 
	double *y);


/**
 * @brief This function updates the open boundaries.
 * @param arrays This pointer variable passes the arrays structure.
 */
__device__ double d_interpolate_matrix(double valor, 
	int line, int nPointsLine,
	int n, int idx0,
	int tidx, 
	double *x, 
	double *y);