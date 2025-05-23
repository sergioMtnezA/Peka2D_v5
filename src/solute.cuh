#include "define.h"
#include "structs.h"
#include "cuTilities.cuh"

#define NON_DIFF 0
#define CONSTANT_DIFF 1
#define ANISOTROPIC_DIFF 2
#define FRICTION_DIFF 3
#define _Dm_ 0.0 //Molecular diffusion

/**
 * @brief This function initializes the variation of conserved variables in the arrays structure to 0.0 and the array solidWallByCell to 0.
 * @param nTasks This integer variable passes the number of times the computation needs to be done.
 * @param arrays This pointer variable passes the arrays structure.
 */
__global__ void g_initialize_solute_delta(int nTasks, t_arrays *arrays);

/**
 * @brief This function initializes the variation of conserved variables in the arrays structure to 0.0 and the array solidWallByCell to 0.
 * @param nTasks This integer variable passes the number of times the computation needs to be done.
 * @param localDt This pointer variable gets the values of the local time steps.
 */
__global__ void g_wall_solute_calculus(int nTasks, t_arrays *arrays, double *dt);

/**
 * @brief This function initializes the variation of conserved variables in the arrays structure to 0.0 and the array solidWallByCell to 0.
 * @param nTasks This integer variable passes the number of times the computation needs to be done.
 */
__global__ void g_bound_solute_calculus(int nTasks, t_arrays *arrays);

/**
 * @brief This function computes the right-hand side of the equations to update each physical variable.
 * @param nTasks This integer variable passes the number of times the computation needs to be done.
 * @param arrays This pointer variable passes the arrays structure. arrays->dt is modifed.
 */
__global__ void g_update_solute_contributions(int nTasks, t_arrays *arrays);

/**
 * @brief This function updates the water depth and water discharges (in wet cells).
 * @param nTasks This integer variable passes the number of times the computation needs to be done.
 * @param arrays This pointer variable passes the arrays structure.
 */
__global__ void g_update_solute_cells(int nTasks, t_arrays *arrays);

/**
 * @brief This function initializes the variation of conserved variables in the arrays structure to 0.0 and the array solidWallByCell to 0.
 * @param nTasks This integer variable passes the number of times the computation needs to be done.
 * @param localDt This pointer variable gets the values of the local time steps.
 */
__global__ void g_initialize_solute_diffusion_delta(int nTasks, t_arrays *arrays);

/**
 * @brief This function initializes the variation of conserved variables in the arrays structure to 0.0 and the array solidWallByCell to 0.
 * @param nTasks This integer variable passes the number of times the computation needs to be done.
 * @param localDt This pointer variable gets the values of the local time steps.
 */
__global__ void g_wall_solute_diffusion_calculus(int nTasks, t_arrays *arrays);

/**
 * @brief This function computes the right-hand side of the equations to update each physical variable.
 * @param nTasks This integer variable passes the number of times the computation needs to be done.
 * @param arrays This pointer variable passes the arrays structure. arrays->dt is modifed.
 */
__global__ void g_update_solute_diffusion_contributions(int nTasks, t_arrays *arrays, double *localDtd);

/**
 * @brief This function updates the water depth and water discharges (in wet cells).
 * @param nTasks This integer variable passes the number of times the computation needs to be done.
 * @param arrays This pointer variable passes the arrays structure.
 */
__global__ void g_get_solute_diffusion_dtmin(t_arrays *arrays, double *localDtd, int *idmin);

/**
 * @brief This function updates the water depth and water discharges (in wet cells).
 * @param nTasks This integer variable passes the number of times the computation needs to be done.
 * @param arrays This pointer variable passes the arrays structure.
 */
__global__ void g_update_solute_diffusion_cells(int nTasks, t_arrays *arrays, double *Dtd);



